#!/usr/bin/env python3
"""Run the preserved five-algorithm benchmark from a relocatable package."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import xml.etree.ElementTree as ET

BUNDLE = Path(__file__).resolve().parent.parent


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def link(target, path):
    path.symlink_to(os.path.relpath(target, path.parent), target_is_directory=True)


def prepare(output, smoke, native_binaries=None):
    runtime_path = BUNDLE / 'environment/runtime.json'
    if not runtime_path.exists():
        raise SystemExit('Run python3 install_dependencies.py first.')
    runtime = json.loads(runtime_path.read_text())
    if runtime.get('status') != 'ready':
        raise SystemExit('Python environments are not ready; run install_dependencies.py.')
    output.mkdir(parents=True, exist_ok=False)
    runtime_snapshot = output / 'environment_runtime.json'
    shutil.copy2(runtime_path, runtime_snapshot)
    original = json.loads((BUNDLE / 'reference_results/metadata.json').read_text())
    count = 3 if smoke else 1000
    for name in ['sources', 'support']:
        link(BUNDLE / name, output / name)
    native = (BUNDLE / native_binaries).resolve() if native_binaries else BUNDLE / 'bin'
    assert all((native / name).is_file() for name in ['ebi-native-dynamic-smoke' if smoke else 'ebi-native-dynamic-benchmark', 'standard-alignments-benchmark'])
    link(native, output / 'bin')
    if not smoke:
        link(BUNDLE / 'inputs', output / 'inputs')
    else:
        # Only the functional smoke test shortens inputs. The full benchmark
        # always consumes the exact untouched archived XES and manifests.
        inputs = output / 'inputs'
        inputs.mkdir()
        link(BUNDLE / 'inputs/models', inputs / 'models')
        (inputs / 'reference').mkdir()
        import csv
        for year in ['2012', '2017']:
            target = inputs / f'bpi{year}'
            target.mkdir()
            source = BUNDLE / f'inputs/bpi{year}'
            records = json.loads((source / 'sample_manifest.json').read_text())[:count]
            save(target / 'sample_manifest.json', records)
            tree = ET.parse(source / 'sampled_1000_unique_traces.xes')
            root = tree.getroot()
            if root.tag.startswith('{'):
                ET.register_namespace('', root.tag[1:].split('}')[0])
            traces = [c for c in root if c.tag.rsplit('}', 1)[-1] == 'trace']
            assert len(traces) == 1000
            for trace in traces[count:]:
                root.remove(trace)
            tree.write(target / 'sampled_1000_unique_traces.xes', encoding='utf-8', xml_declaration=True)
            for level in ['00', '10', '25', '50']:
                name = f'bpi{year}_pt{level}.csv'
                with (BUNDLE / 'inputs/reference' / name).open(newline='') as stream:
                    reader = csv.DictReader(stream)
                    fields = reader.fieldnames
                    rows = list(reader)[:count]
                with (inputs / 'reference' / name).open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(rows)
                cost_name = f'bpi{year}_pt{level}_costs.json'
                costs = json.loads((BUNDLE / 'inputs/reference' / cost_name).read_text())[:count]
                save(inputs / 'reference' / cost_name, costs)
    meta = copy.deepcopy(original)
    for key in ['started_utc', 'finished_utc', 'verification', 'current_job', 'current_job_started_utc',
                'failed_job', 'returncode']:
        meta.pop(key, None)
    # Preserve the venv executable path: resolving its symlink would execute
    # the base interpreter directly and bypass the archived environment.
    executables = {a: str(BUNDLE / name) for a, name in runtime['python_executables'].items()}
    assert all(Path(path).is_file() for path in executables.values())
    protected = {str(runtime_snapshot): sha(runtime_snapshot)}
    for directory in ['inputs', 'bin', 'sources', 'support']:
        for path in (output / directory).rglob('*'):
            if path.is_file():
                protected[str(path.resolve())] = sha(path)
    for path in (BUNDLE / 'scripts').glob('*.py'):
        protected[str(path)] = sha(path)
    cpp_module = (BUNDLE / runtime['cpp_module']).resolve()
    assert sha(cpp_module) == runtime['cpp_module_sha256'], 'Selected C++ module changed after setup'
    protected[str(cpp_module)] = sha(cpp_module)
    for key in ['reference_result_kind', 'source_campaigns', 'assembled_utc']:
        meta.pop(key, None)
    meta.update(status='prepared', created_utc=datetime.now(timezone.utc).isoformat(),
                completed_jobs=[], samples_per_dataset=count,
                expected_timed_alignments=len(meta['pairs']) * count * 5 * 3,
                python_executables=executables, protected_files_sha256=protected,
                smoke_test=smoke, reproduction_package=str(BUNDLE),
                historical_sources=['best_of_three_run_01 (four comparators)', 'native_dynamic_best_of_three_run_03 (native DP)'],
                reference_metadata_sha256=sha(BUNDLE / 'reference_results/metadata.json'),
                runtime_configuration_sha256=sha(runtime_snapshot),
                runtime_configuration='environment_runtime.json',
                cpp_module=str((BUNDLE / runtime['cpp_module']).resolve()),
                cpp_module_sha256=runtime['cpp_module_sha256'],
                cpp_original_binary_reused=runtime['cpp_original_binary_reused'],
                native_binary_selection=str(native),
                historical_ebi_binaries_reused=native == BUNDLE / 'bin',
                logical_cpu_count=os.cpu_count(), platform=platform.platform(),
                cpu_affinity=sorted(os.sched_getaffinity(0)),
                cpu_models=sorted({s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                                   if s.startswith('model name')}),
                controller_python=sys.version,
                experiment_kind='functional smoke test; not a performance result' if smoke else 'full five-algorithm protocol with native DP replacement')
    save(output / 'metadata.json', meta)
    print(f'Prepared {count} traces per dataset, 8 models, 5 algorithms, 3 repetitions: '
          f'{meta["expected_timed_alignments"]:,} timed calls. Output: {output}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', help='New output directory inside this package; default runs/reproduction_01')
    parser.add_argument('--smoke', action='store_true', help='Test first three traces per dataset across all eight models and five algorithms')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--resume', action='store_true', help='Resume a prepared/interrupted run between completed jobs')
    parser.add_argument('--report-only', action='store_true', help='Rebuild tables from an existing run')
    parser.add_argument('--native-binaries', help='Optional rebuilt Rust runner directory; default uses preserved historical executables')
    args = parser.parse_args()
    name = args.output or ('runs/native_smoke_01' if args.smoke else 'runs/reproduction_01')
    output = Path(name)
    if not output.is_absolute():
        output = BUNDLE / output
    output = output.resolve()
    if not output.is_relative_to(BUNDLE / 'runs') or output == BUNDLE / 'runs':
        parser.error('Use a new subdirectory of the package runs/ directory.')
    if args.report_only or args.resume:
        if not (output / 'metadata.json').exists():
            parser.error('No prepared run at this output path.')
    else:
        if output.exists():
            parser.error('Output already exists. Use --resume, --report-only, or a new output name.')
        prepare(output, args.smoke, args.native_binaries)
    if args.prepare_only:
        return
    os.environ['EBI_REPRO_RUN_DIR'] = str(output)
    path = BUNDLE / 'scripts/benchmark_engine.py'
    spec = importlib.util.spec_from_file_location('portable_benchmark_engine', path)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    if args.report_only:
        engine.report()
    else:
        engine.run_all()


if __name__ == '__main__':
    main()
