"""Run the unchanged Python dyn_align on the existing seed-42 benchmark samples."""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import gc
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import statistics
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

BATCH = Path(__file__).resolve().parent
PROJECT = BATCH.parent.parent
SOURCE = PROJECT.parent / 'process-tree-alignments-refactoring'
SAMPLES = PROJECT / 'benchmark_results/rerun_02'
REFERENCE = PROJECT / 'benchmark_results/optimized_special_alignments_run_03'
IDS = ['sample_index', 'source_trace_index', 'case_id', 'variant_frequency', 'trace_length', 'trace_sha256']
FIELDS = IDS + ['elapsed_ns', 'elapsed_ms', 'alignment_cost', 'ebi_alignment_cost',
                'cost_matches_ebi', 'status', 'attempt_elapsed_ns', 'error']
TIMEOUT = 60


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fields, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def tag(element):
    return element.tag.rsplit('}', 1)[-1]


def read_sample(path, manifest):
    traces = []
    for _, element in ET.iterparse(path, events=('end',)):
        if tag(element) != 'trace':
            continue
        case = next((c.get('value') for c in element if tag(c) == 'string' and c.get('key') == 'concept:name'), '')
        activities = []
        for event in element:
            if tag(event) == 'event':
                label = next((c.get('value') for c in event if tag(c) == 'string' and c.get('key') == 'concept:name'), None)
                assert label is not None, 'Event missing concept:name'
                activities.append(label)
        record = manifest[len(traces)]
        assert record['sample_index'] == len(traces) + 1
        assert case == record['case_id'] and activities == record['activities']
        assert len(activities) == record['trace_length']
        digest = hashlib.sha256(json.dumps(activities, ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()
        assert digest == record['trace_sha256']
        traces.append(tuple(activities))
        element.clear()
    assert len(traces) == len(manifest) == len(set(traces))
    return traces


def canonical(operator, label, children):
    if operator is None:
        return ('leaf', label)
    if operator in ('X', '+'):
        children = sorted(children, key=repr)
    return (operator, tuple(children))


def check_model_import(path, tree):
    """Allow commutative child order and removal of an explicit silent loop exit."""
    source = ET.parse(path).getroot()[0]
    nodes = {node.get('id'): node for node in source if node.get('name') is not None}
    children = {key: [] for key in nodes}
    for edge in source:
        if edge.get('sourceId') is not None:
            children[edge.get('sourceId')].append(edge.get('targetId'))
    operators = {'sequence': '->', 'and': '+', 'xor': 'X', 'xorLoop': '*'}

    def raw(key):
        node = nodes[key]
        if tag(node) in ('manualTask', 'automaticTask'):
            return canonical(None, node.get('name') if tag(node) == 'manualTask' else None, [])
        descendants = [raw(child) for child in children[key]]
        if tag(node) == 'xorLoop' and len(descendants) == 3:
            assert descendants[2] == ('leaf', None), 'Unexpected non-silent loop exit'
            descendants.pop()
        return canonical(operators[tag(node)], None, descendants)

    def imported(node):
        return canonical(str(node.operator) if node.operator else None, node.label, [imported(c) for c in node.children])

    assert raw(source.get('root')) == imported(tree), 'Imported model differs from source PTML'


class AlignmentTimeout(Exception):
    pass


def alarm_handler(signum, frame):
    raise AlignmentTimeout(f'Exceeded {TIMEOUT}-second per-trace limit')


def measure(align, trace, tree, expected):
    # Alarm setup/reset and validation are outside the measured call.
    signal.setitimer(signal.ITIMER_REAL, TIMEOUT)
    started = time.perf_counter_ns()
    try:
        cost = align(trace, tree)
        elapsed = time.perf_counter_ns() - started
    except Exception as exception:
        elapsed = time.perf_counter_ns() - started
        return {'elapsed_ns': '', 'elapsed_ms': '', 'alignment_cost': '',
                'ebi_alignment_cost': expected, 'cost_matches_ebi': '',
                'status': 'timeout' if isinstance(exception, AlignmentTimeout) else 'error',
                'attempt_elapsed_ns': elapsed, 'error': f'{type(exception).__name__}: {exception}'}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    assert isinstance(cost, int) and cost >= 0, f'Invalid alignment cost: {cost!r}'
    return {'elapsed_ns': elapsed, 'elapsed_ms': str(Decimal(elapsed) / 1_000_000),
            'alignment_cost': cost, 'ebi_alignment_cost': expected,
            'cost_matches_ebi': int(cost == expected), 'status': 'ok' if cost == expected else 'cost_mismatch',
            'attempt_elapsed_ns': elapsed, 'error': ''}


def summarize(rows):
    times = [int(row['elapsed_ns']) / 1_000_000 for row in rows if row['elapsed_ns'] != '']
    return {'attempted_traces': len(rows), 'completed_traces': len(times),
            'cost_matches_ebi': sum(row['status'] == 'ok' for row in rows),
            'cost_mismatches': sum(row['status'] == 'cost_mismatch' for row in rows),
            'timeouts': sum(row['status'] == 'timeout' for row in rows),
            'errors': sum(row['status'] == 'error' for row in rows),
            'mean_ms': statistics.mean(times) if times else '',
            'median_ms': statistics.median(times) if times else '',
            'p95_ms': statistics.quantiles(times, n=100, method='inclusive')[94] if len(times) > 1 else '',
            'max_ms': max(times) if times else '',
            'sum_alignment_seconds': sum(int(row['elapsed_ns']) for row in rows if row['elapsed_ns'] != '') / 1e9,
            'sum_attempt_seconds': sum(int(row['attempt_elapsed_ns']) for row in rows) / 1e9}


def worker(name):
    import pm4py
    directory = BATCH / name
    meta = json.loads((directory / 'metadata.json').read_text())
    root_meta = json.loads((BATCH / 'metadata.json').read_text())
    assert meta['status'] == 'prepared'
    for path, digest in meta['inputs_sha256'].items():
        assert sha256(PROJECT / path) == digest, path
    snapshot = BATCH / 'source_snapshot/dynamic_alignment.py'
    assert sha256(snapshot) == root_meta['algorithm_source_sha256']
    spec = importlib.util.spec_from_file_location('benchmarked_dynamic_alignment', snapshot)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    manifest = json.loads((PROJECT / meta['sample_manifest']).read_text())
    traces = read_sample(PROJECT / meta['sample_xes'], manifest)
    reference = read_csv(PROJECT / meta['reference_csv'])
    assert len(reference) == len(manifest)
    for record, row in zip(manifest, reference):
        assert all(str(record[field] if record[field] is not None else '') == row[field] for field in IDS)
    tree = pm4py.read_ptml(str(PROJECT / meta['model_ptml']))
    check_model_import(PROJECT / meta['model_ptml'], tree)
    (directory / 'model_tree_string.txt').write_text(str(tree) + '\n')
    signal.signal(signal.SIGALRM, alarm_handler)
    meta.update(status='running', started_utc=now(), python=sys.version,
                pm4py_version=pm4py.__version__, platform=platform.platform(),
                garbage_collection_enabled=gc.isenabled(), gc_threshold=gc.get_threshold(),
                source_ptml_structure_verified=True, sample_xes_manifest_verified=True,
                logical_cpu_count=os.cpu_count(), cpu_affinity=sorted(os.sched_getaffinity(0)),
                cpu_models=sorted({line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')}))
    save(directory / 'metadata.json', meta)
    dependencies = sorted(f'{dist.metadata["Name"]}=={dist.version}' for dist in importlib.metadata.distributions())
    (directory / 'environment.txt').write_text('\n'.join(dependencies) + '\n')
    meta['warmup_results'] = []
    for index in range(3):
        outcome = measure(module.dyn_align, traces[index], tree, int(reference[index]['alignment_cost']))
        meta['warmup_results'].append({'sample_index': index + 1, **outcome})
    save(directory / 'metadata.json', meta)
    print(f'{name}: 1,000 verified unique traces; three warm-ups complete', flush=True)
    timings = directory / 'python_alignment_timings.csv'
    with timings.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for index, (trace, record, baseline) in enumerate(zip(traces, manifest, reference), 1):
            row = {field: record[field] for field in IDS}
            row.update(measure(module.dyn_align, trace, tree, int(baseline['alignment_cost'])))
            writer.writerow(row)
            stream.flush()
            if index % 25 == 0 or row['status'] != 'ok':
                print(f'{name}: {index}/1000; latest status={row["status"]}, ms={row["elapsed_ms"]}', flush=True)
    rows = read_csv(timings)
    meta.update(status='completed', finished_utc=now(), csv_sha256=sha256(timings), summary=summarize(rows))
    save(directory / 'metadata.json', meta)
    print(json.dumps(meta['summary']), flush=True)


def report():
    meta = json.loads((BATCH / 'metadata.json').read_text())
    combined, summaries = [], []
    for name in meta['run_order']:
        directory = BATCH / name
        child = json.loads((directory / 'metadata.json').read_text())
        assert child['status'] == 'completed'
        path = directory / 'python_alignment_timings.csv'
        assert sha256(path) == child['csv_sha256']
        rows = read_csv(path)
        manifest = json.loads((PROJECT / child['sample_manifest']).read_text())
        reference = read_csv(PROJECT / child['reference_csv'])
        assert len(rows) == len(manifest) == len(reference) == 1000
        assert len({r['trace_sha256'] for r in rows}) == 1000
        for row, record, baseline in zip(rows, manifest, reference):
            assert all(row[k] == str(record[k]) == baseline[k] for k in IDS)
            assert row['ebi_alignment_cost'] == baseline['alignment_cost']
            assert row['status'] in ('ok', 'cost_mismatch', 'timeout', 'error')
            if row['elapsed_ns']:
                assert int(row['elapsed_ns']) > 0
                assert Decimal(row['elapsed_ms']) * 1_000_000 == int(row['elapsed_ns'])
                assert (row['alignment_cost'] == baseline['alignment_cost']) == (row['status'] == 'ok')
            else:
                assert row['status'] in ('timeout', 'error') and row['alignment_cost'] == ''
            combined.append({'dataset': child['dataset'], 'model': child['model'], **row})
        summaries.append({'dataset': child['dataset'], 'model': child['model'], **summarize(rows)})
    assert len(combined) == 8000
    for path, digest in meta['protected_files_sha256'].items():
        assert sha256(Path(path)) == digest, f'Input/source changed during run: {path}'
    write_csv(BATCH / 'all_trace_timings.csv', ['dataset', 'model'] + FIELDS, combined)
    write_csv(BATCH / 'all_datasets_summary.csv', list(summaries[0]), summaries)
    write_csv(BATCH / 'failures_and_cost_mismatches.csv', ['dataset', 'model'] + FIELDS,
              [row for row in combined if row['status'] != 'ok'])
    assert read_csv(BATCH / 'all_trace_timings.csv') == [{k: str(v) for k, v in row.items()} for row in combined]
    meta.update(status='completed', finished_utc=now(), summary=summarize(combined),
                input_and_source_hashes_verified_after_run=True,
                outputs_sha256={p.name: sha256(p) for p in BATCH.glob('*.csv')})
    save(BATCH / 'metadata.json', meta)
    lines = [
        '# Python process-tree-alignments-refactoring benchmark', '',
        'This is a new run of the unchanged `dyn_align` implementation from the sibling '
        '`process-tree-alignments-refactoring` project. It measures all 1,000 existing seed-42 '
        'unique traces per year against pt00, pt10, pt25 and pt50 (8,000 trace/model pairs).', '',
        '## Method', '',
        '- Saved XES samples and manifests come from `../rerun_02`; no traces were resampled or reordered. '
        'Every case ID, activity sequence, trace length and SHA-256 was checked. The same 1,000 traces '
        'are used for all four models within each year; all lifecycle events are retained.',
        '- Original PTML files in `DATAFILES/ptml_python` are imported once with PM4Py. Imported tree '
        'structure is checked against source PTML, allowing commutative child order and silent loop exits.',
        '- Eight models run sequentially in fresh Python processes, with three warm-ups on samples 1–3 '
        'followed by one timed call for each of the 1,000 samples. No benchmark workers run concurrently.',
        '- `time.perf_counter_ns()` measures `dyn_align(tuple_trace, tree)`. XES/PTML parsing, tuple '
        'creation, cost validation, timeout setup/reset and CSV output are outside the timer. Normal '
        'Python garbage collection remains enabled. Label sets remain cached after warm-up; unique-label '
        'validation runs inside each call, and DP cost caches are new for each call.',
        '- Python returns a cost only, whereas the Ebi implementations also reconstruct alignments. '
        'Costs use one per visible deviation and zero per silent move. Every returned cost is compared '
        'with `optimized_special_alignments_run_03`; this does not validate an alignment path.',
        f'- A {TIMEOUT}-second per-trace wall-time limit matches the previous C++/default benchmark limit. '
        'Timeouts/errors, if any, have blank completed times/costs and a separate `attempt_elapsed_ns`. '
        'Timing statistics use completed calls only, including any completed cost mismatches.',
        '- Python 3.14.7 and PM4Py 2.7.23.8 from the project venv; `PYTHONHASHSEED=42`. '
        'Exact installed packages and machine details are recorded per model.',
        '- `source_snapshot/dynamic_alignment.py` is a byte-for-byte copy of the Python source and is '
        'loaded directly to avoid importing the unrelated Gurobi/graph alignment implementations. '
        'The project evaluation driver is not used because it selects variants independently, runs '
        'workers concurrently and takes the best of multiple repetitions.',
        '- SHA-256 hashes of input logs, PTMLs, saved samples, manifests, reference CSVs, the algorithm '
        'and runner are checked before and after execution. Earlier benchmark results are preserved.', '',
        '## Results', '',
        '| Dataset | Model | Completed | Cost matches | Timeouts | Errors | Mean (ms) | Median (ms) | Max (ms) |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for row in summaries:
        lines.append(f'| {row["dataset"]} | {row["model"]} | {row["completed_traces"]} | '
                     f'{row["cost_matches_ebi"]} | {row["timeouts"]} | {row["errors"]} | '
                     f'{row["mean_ms"]:.3f} | {row["median_ms"]:.3f} | {row["max_ms"]:.3f} |')
    summary = meta['summary']
    lines += ['', f'Total completed-call time: **{summary["sum_alignment_seconds"]:.3f} seconds**. '
              f'Cost mismatches: **{summary["cost_mismatches"]}**.', '', '## Files', '',
              '- [All 8,000 per-trace timings](all_trace_timings.csv), including dataset/model and original trace IDs.',
              '- [Per-model summary](all_datasets_summary.csv).',
              '- [Failures and cost mismatches](failures_and_cost_mismatches.csv) (header only if none).',
              '- Each model directory contains `python_alignment_timings.csv`, `run.log`, `metadata.json`, '
              '`environment.txt` and the imported `model_tree_string.txt`.',
              '- [Run metadata](metadata.json), [runner](run_benchmark.py) and [source snapshot](source_snapshot/dynamic_alignment.py).', '',
              'Regenerate the combined files and report without rerunning alignments: '
              '`python3 -B benchmark_results/python_refactoring_run_01/run_benchmark.py --report`.', '']
    (BATCH / 'README.md').write_text('\n'.join(lines))
    print(json.dumps(meta['summary']), flush=True)


def run_all():
    assert not (BATCH / 'metadata.json').exists(), 'Refusing to overwrite an existing run'
    snapshot = BATCH / 'source_snapshot'
    snapshot.mkdir()
    protected = {}

    def protect(path, expected=None):
        digest = sha256(path)
        assert expected is None or digest == expected, path
        protected[str(path)] = digest
        return digest

    source_file = SOURCE / 'process_tree_alignments/dynamic_alignment.py'
    for path in (source_file, SOURCE / 'evaluation.py', SOURCE / 'pyproject.toml', SOURCE / 'uv.lock'):
        protect(path)
        shutil.copy2(path, snapshot / path.name)
    protect(Path(__file__).resolve())
    algorithm_hash = protect(snapshot / 'dynamic_alignment.py')
    assert algorithm_hash == protected[str(source_file)]
    reference_meta = json.loads((REFERENCE / 'metadata.json').read_text())
    protect(REFERENCE / 'metadata.json')
    by_model = {(row['dataset'], row['model']): row for row in reference_meta['runs']}
    names, sample_pairs = [], {}
    for year in (2012, 2017):
        for level in ('00', '10', '25', '50'):
            name = f'bpi{year}_pt{level}_seed42_1000'
            names.append(name)
            original = SAMPLES / name
            old = json.loads((original / 'metadata.json').read_text())
            paths = [(PROJECT / old['source_xes'], old['source_xes_sha256']),
                     (PROJECT / old['model_ptml'], old['model_ptml_sha256']),
                     (original / 'sampled_1000_unique_traces.xes', old['sample_xes_sha256']),
                     (original / 'sample_manifest.json', old['sample_manifest_sha256']),
                     (REFERENCE / name / 'special_alignment_timings.csv', by_model[(year, f'pt{level}')]['csv_sha256'])]
            for path, expected in paths:
                if str(path) in protected:
                    assert protected[str(path)] == expected
                else:
                    protect(path, expected)
            pair = (paths[2][1], paths[3][1])
            assert year not in sample_pairs or sample_pairs[year] == pair
            sample_pairs[year] = pair
            child = {key: old[key] for key in ('source_xes', 'model_ptml', 'seed', 'sample_size', 'sampling',
                                              'representative_case', 'activity_attribute', 'event_filter', 'unique_trace_definition')}
            child.update(status='prepared', dataset=year, model=f'pt{level}',
                         sample_xes=str(paths[2][0].relative_to(PROJECT)),
                         sample_manifest=str(paths[3][0].relative_to(PROJECT)),
                         reference_csv=str(paths[4][0].relative_to(PROJECT)),
                         inputs_sha256={str(path.relative_to(PROJECT)): digest for path, digest in paths[1:]})
            directory = BATCH / name
            directory.mkdir()
            save(directory / 'metadata.json', child)
    python = str(SOURCE / 'venv/bin/python')
    meta = dict(status='prepared', created_utc=now(), source_directory=str(SOURCE),
                algorithm_source_sha256=algorithm_hash, python_executable=python,
                interface='dynamic_alignment.dyn_align(tuple[str, ...], pm4py.ProcessTree) -> int',
                warmup_traces=3, measurements_per_trace=1, benchmark_parallelism=1,
                timeout_seconds=TIMEOUT, cross_trace_dp_cache=False, label_cache='retained on model after warm-up',
                python_hash_seed=42, timer='time.perf_counter_ns; wall time around dyn_align only',
                run_order=names, completed_runs=[], protected_files_sha256=protected)
    save(BATCH / 'metadata.json', meta)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='42',
               MPLCONFIGDIR='/tmp/ebi_python_benchmark_matplotlib')
    meta.update(status='running', started_utc=now())
    save(BATCH / 'metadata.json', meta)
    for name in names:
        command = [python, '-B', str(Path(__file__).resolve()), '--worker', name]
        print(f'Starting {name}', flush=True)
        with (BATCH / name / 'run.log').open('x') as stream:
            result = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            meta.update(status='failed', failed_model=name, returncode=result.returncode)
            save(BATCH / 'metadata.json', meta)
            raise RuntimeError(f'{name} failed; see its run.log')
        meta['completed_runs'].append({'name': name, 'command': command, 'finished_utc': now()})
        save(BATCH / 'metadata.json', meta)
        child = json.loads((BATCH / name / 'metadata.json').read_text())
        print(f'Finished {name}: {json.dumps(child["summary"])}', flush=True)
    report()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--worker')
    group.add_argument('--report', action='store_true')
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments.worker)
    elif arguments.report:
        report()
    else:
        run_all()
