"""Run the five-algorithm benchmark after installing its pinned dependencies."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def verify_package(root):
    manifest = json.loads((root / 'package_manifest.json').read_text())
    require(manifest['format_version'] == 2, 'Unsupported package manifest')
    for relative, record in manifest['files'].items():
        path = root / relative
        require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(root),
                f'Missing or redirected package file: {relative}')
        require(sha(path) == record['sha256'], f'Changed package file: {relative}')
    return manifest


def trace_identity(element):
    case = next((c.get('value') for c in element if c.tag.rsplit('}', 1)[-1] == 'string'
                 and c.get('key') == 'concept:name'), '')
    activities = []
    for event in element:
        if event.tag.rsplit('}', 1)[-1] == 'event':
            labels = [c.get('value') for c in event if c.get('key') == 'concept:name']
            require(len(labels) == 1, 'Missing/ambiguous event concept:name')
            activities.append(labels[0])
    return case, activities


def prepare_samples(work, datasets, manifest):
    metadata = json.loads((work / 'reference_results/metadata.json').read_text())
    report_path = work / 'sample_reconstruction.json'
    previous = json.loads(report_path.read_text()) if report_path.exists() else {}
    report = dict(status='running', datasets={})
    for year in ['2012', '2017']:
        source = datasets[year]
        expected = metadata['samples'][year]['source_xes_sha256']
        print(f'{year}: checking the supplied full XES log...', flush=True)
        require(sha(source) == expected,
                f'{source.name} differs from the original benchmark dataset. Expected SHA-256 {expected}.')
        output = work / f'inputs/bpi{year}/sampled_1000_unique_traces.xes'
        saved = previous.get('datasets', {}).get(year, {})
        if saved.get('source_sha256') == expected and output.exists() and sha(output) == manifest['original_sample_hashes'][year]:
            report['datasets'][year] = saved
            print(f'{year}: saved 1,000-trace selection already verified.', flush=True)
            continue
        records = json.loads((output.parent / 'sample_manifest.json').read_text())
        require(len(records) == len({tuple(r['activities']) for r in records}) == 1000, 'Invalid saved selection')
        selection = {int(r['source_trace_index']): i for i, r in enumerate(records)}
        require(len(selection) == 1000, 'Duplicate source indices')
        selected = [None] * 1000
        parser = ET.iterparse(source, events=('start', 'end'))
        depth, source_index, out = 0, 0, None
        for event, element in parser:
            if event == 'start':
                depth += 1
                if depth == 1:
                    original_root = element
                    out = ET.Element(element.tag, element.attrib)
                    if element.tag.startswith('{'):
                        ET.register_namespace('', element.tag[1:].split('}')[0])
                continue
            if depth == 2:
                if element.tag.rsplit('}', 1)[-1] == 'trace':
                    source_index += 1
                    if source_index in selection:
                        index = selection[source_index]
                        record = records[index]
                        case, activities = trace_identity(element)
                        encoded = json.dumps(activities, ensure_ascii=False, separators=(',', ':')).encode()
                        require(case == record['case_id'] and activities == record['activities']
                                and len(activities) == record['trace_length']
                                and hashlib.sha256(encoded).hexdigest() == record['trace_sha256'],
                                f'{year}: selected trace {index + 1} does not match the saved manifest')
                        selected[index] = copy.deepcopy(element)
                elif element.tag.rsplit('}', 1)[-1] in ['extension', 'global', 'classifier']:
                    out.append(copy.deepcopy(element))
                original_root.remove(element)
                element.clear()
            depth -= 1
        require(all(t is not None for t in selected), f'{year}: some selected source traces were not found')
        out.extend(selected)  # Saved sample order, not source-log order.
        temporary = output.with_suffix('.xes.tmp')
        ET.ElementTree(out).write(temporary, encoding='utf-8', xml_declaration=True)
        temporary.replace(output)
        actual_hash = sha(output)
        require(actual_hash == manifest['original_sample_hashes'][year],
                f'{year}: reconstructed XES differs from the saved benchmark sample')
        report['datasets'][year] = dict(source_file=str(source), source_sha256=expected,
                                       source_trace_count=source_index, selected_unique_traces=1000,
                                       sampled_xes_sha256=actual_hash,
                                       byte_identical_to_original_sample=actual_hash == manifest['original_sample_hashes'][year],
                                       all_selected_case_ids_activities_indices_and_order_match=True)
        print(f'{year}: reconstructed 1,000 unique traces in their original sample order.', flush=True)
    report.update(status='passed', checked_utc=datetime.now(timezone.utc).isoformat())
    save(report_path, report)


def child_environment():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='42', RAYON_NUM_THREADS='1')
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    return env


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description='Five algorithms, saved 1,000 traces per dataset, eight models, minimum of three runs.')
    parser.add_argument('--data-dir', type=Path, default=root, help='Folder containing BPI_Challenge_2012.xes and BPI_Challenge_2017.xes')
    parser.add_argument('--smoke', action='store_true', help='Functional test: first three selected traces, all models/algorithms, three runs')
    parser.add_argument('--prepare-only', action='store_true', help='Prepare inputs and a run without alignments; dependencies must already be installed')
    parser.add_argument('--resume', action='store_true', help='Resume a prepared run between completed jobs')
    parser.add_argument('--output', help='New run name; default reproduction_01 or smoke_01')
    args = parser.parse_args()
    require(platform.python_implementation() == 'CPython' and platform.python_version() == '3.14.7'
            and platform.system() == 'Linux' and platform.machine() == 'x86_64',
            'This preserved build requires CPython 3.14.7 on compatible Linux x86_64. Run it using python3.14.')
    name = args.output or ('smoke_01' if args.smoke else 'reproduction_01')
    require(name not in ['', '.', '..'] and Path(name).name == name, '--output must be a single directory name')
    work = root
    result = work / 'runs' / name
    require(args.resume or not result.exists(), f'Run already exists: {result}. Choose a new --output name or --resume.')
    datasets = {year: (args.data_dir / f'BPI_Challenge_{year}.xes').resolve() for year in ['2012', '2017']}
    missing = [str(path) for path in datasets.values() if not path.is_file()]
    require(not missing, 'Add both original BPI datasets. Missing: ' + ', '.join(missing))
    manifest = verify_package(root)
    runtime = root / 'environment/runtime.json'
    require(runtime.is_file() and json.loads(runtime.read_text()).get('status') == 'ready',
            'Install dependencies first: python3.14 install_dependencies.py')
    subprocess.run([sys.executable, '-B', str(root / 'install_dependencies.py'), '--verify-only'],
                   cwd=root, env=child_environment(), check=True)
    prepare_samples(root, datasets, manifest)
    command = [sys.executable, '-B', str(work / 'scripts/run_benchmark.py'), '--output', f'runs/{name}']
    for flag, enabled in [('--smoke', args.smoke), ('--prepare-only', args.prepare_only), ('--resume', args.resume)]:
        if enabled:
            command.append(flag)
    subprocess.run(command, cwd=work, env=child_environment(), check=True)
    if args.prepare_only:
        print(f'Preparation complete. Start using the same command with --resume instead of --prepare-only. Run: {result}', flush=True)
    else:
        subprocess.run([sys.executable, '-B', str(work / 'scripts/verify_results.py'), '--results', f'runs/{name}'],
                       cwd=work, env=child_environment(), check=True)
        print(f'Completed and independently verified. Comparison CSV: {result / "all_programs_trace_timings.csv"}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error
