"""Benchmark locally corrected PM4Py process-tree alignments on the saved Ebi samples."""
import argparse
import csv
from decimal import Decimal
import gc
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

import benchmark_helpers as helpers

BATCH = Path(__file__).resolve().parent
PROJECT = BATCH.parent.parent
PREVIOUS = PROJECT / 'benchmark_results/python_refactoring_run_01'
IDS = helpers.IDS
DETAILS = ['alignment_moves', 'log_moves', 'visible_model_moves', 'silent_model_moves',
           'synchronous_moves', 'alignment_checks_passed', 'reported_optimal']
FIELDS = helpers.FIELDS + DETAILS
ALL_FIELDS = ['dataset', 'model'] + FIELDS
TIMEOUT = 60
RETRY_TIMEOUT = 600


class AlignmentTimeout(Exception):
    pass


def alarm_handler(signum, frame):
    raise AlignmentTimeout('Per-trace wall-time limit exceeded')


def measure(align, trace, tree, expected, limit):
    """Time the unchanged public API; inspect the returned alignment after timing."""
    row = {field: '' for field in FIELDS if field not in IDS}
    row['ebi_alignment_cost'] = expected
    result = None
    signal.setitimer(signal.ITIMER_REAL, limit)
    started = time.perf_counter_ns()
    try:
        result = align(trace, tree)
        elapsed = time.perf_counter_ns() - started
    except Exception as exception:
        elapsed = time.perf_counter_ns() - started
        row.update(status='timeout' if isinstance(exception, AlignmentTimeout) else 'error',
                   attempt_elapsed_ns=elapsed, error=f'{type(exception).__name__}: {exception}')
        return row
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    row.update(elapsed_ns=elapsed, elapsed_ms=str(Decimal(elapsed) / 1_000_000), attempt_elapsed_ns=elapsed)
    try:
        assert isinstance(result, dict), f'Unexpected result: {type(result).__name__}'
        cost = result['cost']
        assert isinstance(cost, (int, float)) and int(cost) == cost and cost >= 0, f'Invalid cost: {cost}'
        row.update(alignment_cost=int(cost), cost_matches_ebi=int(cost == expected),
                   reported_optimal=int(result.get('optimal') is True))
        counts = dict(log_moves=0, visible_model_moves=0, silent_model_moves=0, synchronous_moves=0)
        consumed = []
        leaves = set()
        stack = [tree]
        while stack:
            node = stack.pop()
            if node.operator is None:
                leaves.add(id(node))
            stack.extend(node.children)
        for log, model in result['alignment']:
            if log != '>>':
                consumed.append(log)
            if isinstance(model, str):
                assert model == '>>' and log != '>>', 'Invalid log move'
                counts['log_moves'] += 1
            else:
                assert id(model) in leaves, 'Alignment references an unknown model leaf'
                if log == '>>':
                    counts['silent_model_moves' if model.label is None else 'visible_model_moves'] += 1
                else:
                    assert model.label == log, 'Invalid synchronous move'
                    counts['synchronous_moves'] += 1
        assert consumed == [event['concept:name'] for event in trace], 'Alignment does not consume the input trace'
        assert counts['log_moves'] + counts['visible_model_moves'] == cost, 'Move cost differs from reported cost'
        state = result['state']
        assert state.index == len(trace) and state.state[(id(tree), tree)] == tree.OperatorState.CLOSED
        assert result.get('optimal') is True, 'PM4Py did not report an optimal alignment'
        row.update(**counts, alignment_moves=len(result['alignment']), alignment_checks_passed=1,
                   status='ok' if cost == expected else 'cost_mismatch')
    except Exception as exception:
        row.update(status='validation_error', alignment_checks_passed=0,
                   error=f'{type(exception).__name__}: {exception}')
    return row


def summarize(rows):
    summary = helpers.summarize(rows)
    summary['cost_matches_ebi'] = sum(str(row['cost_matches_ebi']) == '1' for row in rows)
    summary['cost_mismatches'] = sum(str(row['cost_matches_ebi']) == '0' for row in rows)
    summary['fully_validated_alignments'] = sum(row['status'] == 'ok' for row in rows)
    summary['validation_errors'] = sum(row['status'] == 'validation_error' for row in rows)
    summary['alignment_checks_passed'] = sum(str(row['alignment_checks_passed']) == '1' for row in rows)
    return summary


def source_modules(algorithm):
    modules = [algorithm, algorithm.DEFAULT_VARIANT.value]
    for name in (
        'pm4py.algo.conformance.alignments.process_tree.util.search_graph_pt_replay_semantics',
        'pm4py.objects.process_tree.obj', 'pm4py.objects.process_tree.utils.generic',
        'pm4py.objects.petri_net.utils.align_utils', 'pm4py.objects.log.obj',
        'pm4py.util.exec_utils', 'pm4py.objects.process_tree.importer.importer',
        'pm4py.objects.process_tree.importer.variants.ptml',
    ):
        modules.append(importlib.import_module(name))
    return {module.__name__: Path(module.__file__).resolve() for module in modules}


def prepare():
    import pm4py
    from pm4py.algo.conformance.alignments.process_tree import algorithm
    assert not (BATCH / 'metadata.json').exists(), 'Refusing to overwrite an existing run'
    assert algorithm.DEFAULT_VARIANT == algorithm.Variants.SEARCH_GRAPH_PT
    previous = json.loads((PREVIOUS / 'metadata.json').read_text())
    patch = json.loads((BATCH / 'patch_metadata.json').read_text())
    assert helpers.sha256(Path(patch['installed_path'])) == patch['patched_sha256']
    protected = {str(BATCH / 'patch_metadata.json'): helpers.sha256(BATCH / 'patch_metadata.json')}
    for patch_file in (BATCH / 'patch').iterdir():
        protected[str(patch_file)] = helpers.sha256(patch_file)
    snapshot = BATCH / 'source_snapshot'
    snapshot.mkdir()
    module_files = source_modules(algorithm)
    module_hashes = {}
    for name, path in module_files.items():
        digest = helpers.sha256(path)
        protected[str(path)] = digest
        module_hashes[name] = {'path': str(path), 'sha256': digest}
        target = snapshot / (name.replace('.', '/') + '.py')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        protected[str(target)] = digest
    for path in (Path(__file__).resolve(), BATCH / 'benchmark_helpers.py'):
        protected[str(path)] = helpers.sha256(path)
    assert helpers.sha256(BATCH / 'benchmark_helpers.py') == helpers.sha256(PREVIOUS / 'run_benchmark.py')
    run_order = previous['run_order']
    sample_pairs = {}
    for name in run_order:
        original_path = PREVIOUS / name / 'metadata.json'
        old = json.loads(original_path.read_text())
        assert old['status'] == 'completed'
        protected[str(original_path)] = helpers.sha256(original_path)
        for path, digest in old['inputs_sha256'].items():
            full = PROJECT / path
            assert helpers.sha256(full) == digest, full
            protected[str(full)] = digest
        full_log = PROJECT / old['source_xes']
        digest = previous['protected_files_sha256'][str(full_log)]
        if str(full_log) not in protected:
            assert helpers.sha256(full_log) == digest, full_log
            protected[str(full_log)] = digest
        pair = (old['inputs_sha256'][old['sample_xes']], old['inputs_sha256'][old['sample_manifest']])
        assert old['dataset'] not in sample_pairs or sample_pairs[old['dataset']] == pair
        sample_pairs[old['dataset']] = pair
        keys = ['dataset', 'model', 'source_xes', 'model_ptml', 'sample_xes', 'sample_manifest', 'reference_csv',
                'inputs_sha256', 'seed', 'sample_size', 'sampling', 'representative_case',
                'activity_attribute', 'event_filter', 'unique_trace_definition']
        child = {key: old[key] for key in keys}
        child.update(status='prepared', reference_metadata_sha256=helpers.sha256(original_path))
        directory = BATCH / name
        directory.mkdir()
        helpers.save(directory / 'metadata.json', child)
    dependencies = sorted(f'{dist.metadata["Name"]}=={dist.version}' for dist in importlib.metadata.distributions())
    (BATCH / 'environment.txt').write_text('\n'.join(dependencies) + '\n')
    meta = dict(status='prepared', created_utc=helpers.now(), python=sys.version, python_executable=sys.executable,
                pm4py_version=pm4py.__version__, platform=platform.platform(), local_patch=patch,
                implementation_label=f'PM4Py {pm4py.__version__} with local correctness fixes',
                selected_method='Default process-tree alignment (explicitly selected by user)',
                interface='pm4py.algo.conformance.alignments.process_tree.algorithm.apply(Trace, ProcessTree)',
                default_variant=algorithm.DEFAULT_VARIANT.name, algorithm_modules=module_hashes,
                warmup_traces=3, measurements_per_trace=1, benchmark_parallelism=1,
                python_hash_seed=42, timeout_seconds=TIMEOUT, retry_timeout_seconds=RETRY_TIMEOUT,
                cross_trace_alignment_cache=False, timer='time.perf_counter_ns; wall time around public apply call',
                timing_includes=['per-call leaf extraction and activity projection', 'search', 'alignment reconstruction', 'internal cleanup'],
                timing_excludes=['imports', 'XES/PTML parsing', 'input Trace construction', 'timeout setup/reset',
                                 'returned-result validation and destruction', 'CSV output'],
                costs={'visible_log_move': 1, 'visible_model_move': 1, 'silent_model_move': 0, 'synchronous_move': 0},
                run_order=run_order, completed_runs=[], protected_files_sha256=protected)
    helpers.save(BATCH / 'metadata.json', meta)


def worker(name, retry=False):
    import pm4py
    from pm4py.algo.conformance.alignments.process_tree import algorithm
    from pm4py.objects.log.obj import Event, Trace
    directory = BATCH / name
    original_meta = json.loads((directory / 'metadata.json').read_text())
    root = json.loads((BATCH / 'metadata.json').read_text())
    assert original_meta['status'] == ('completed' if retry else 'prepared')
    assert algorithm.DEFAULT_VARIANT.name == root['default_variant']
    for name_, item in root['algorithm_modules'].items():
        assert helpers.sha256(Path(item['path'])) == item['sha256'], name_
    for path, digest in original_meta['inputs_sha256'].items():
        assert helpers.sha256(PROJECT / path) == digest, path
    manifest = json.loads((PROJECT / original_meta['sample_manifest']).read_text())
    sequences = helpers.read_sample(PROJECT / original_meta['sample_xes'], manifest)
    reference = helpers.read_csv(PROJECT / original_meta['reference_csv'])
    assert len(reference) == 1000
    for record, baseline in zip(manifest, reference):
        assert all(str(record[key]) == baseline[key] for key in IDS)
    traces = [Trace([Event({'concept:name': activity}) for activity in word],
                    attributes={'concept:name': record['case_id']}) for word, record in zip(sequences, manifest)]
    assert all('>>' not in sequence for sequence in sequences)
    tree = pm4py.read_ptml(str(PROJECT / original_meta['model_ptml']))
    helpers.check_model_import(PROJECT / original_meta['model_ptml'], tree)
    signal.signal(signal.SIGALRM, alarm_handler)
    meta = dict(original_meta)
    meta.update(status='running', started_utc=helpers.now(), python=sys.version, pm4py_version=pm4py.__version__,
                garbage_collection_enabled=gc.isenabled(), gc_threshold=gc.get_threshold(),
                source_ptml_structure_verified=True, sample_xes_manifest_verified=True,
                logical_cpu_count=os.cpu_count(), cpu_affinity=sorted(os.sched_getaffinity(0)),
                cpu_models=sorted({line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')}))
    meta_path = directory / ('retry_metadata.json' if retry else 'metadata.json')
    filename = 'retry_timings.csv' if retry else 'pm4py_alignment_timings.csv'
    limit = RETRY_TIMEOUT if retry else TIMEOUT
    meta['timeout_seconds'] = limit
    helpers.save(meta_path, meta)
    if not retry:
        (directory / 'model_tree_string.txt').write_text(str(tree) + '\n')
        indices = list(range(1000))
    else:
        original = helpers.read_csv(directory / 'pm4py_alignment_timings.csv')
        indices = [int(row['sample_index']) - 1 for row in original if row['status'] == 'timeout']
        assert indices and helpers.sha256(directory / 'pm4py_alignment_timings.csv') == original_meta['csv_sha256']
        meta['protocol'] = 'Fresh process/model, three warm-ups, then retry only original timeouts in original order'
    meta['warmup_results'] = []
    for index in range(3):
        outcome = measure(algorithm.apply, traces[index], tree, int(reference[index]['alignment_cost']), limit)
        meta['warmup_results'].append({'sample_index': index + 1, **outcome})
    helpers.save(meta_path, meta)
    print(f'{name}: three warm-ups complete; measuring {len(indices)} trace(s)', flush=True)
    with (directory / filename).open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for count, index in enumerate(indices, 1):
            row = {field: manifest[index][field] for field in IDS}
            row.update(measure(algorithm.apply, traces[index], tree, int(reference[index]['alignment_cost']), limit))
            writer.writerow(row)
            stream.flush()
            if count % 25 == 0 or row['status'] != 'ok' or retry:
                print(f'{name}: {count}/{len(indices)}, sample {index + 1}, status={row["status"]}, '
                      f'ms={row["elapsed_ms"]}, cost={row["alignment_cost"]}, error={row["error"]}', flush=True)
    rows = helpers.read_csv(directory / filename)
    meta.update(status='completed', finished_utc=helpers.now(), csv_sha256=helpers.sha256(directory / filename), summary=summarize(rows))
    helpers.save(meta_path, meta)
    print(json.dumps(meta['summary']), flush=True)


def report():
    root = json.loads((BATCH / 'metadata.json').read_text())
    originals, combined, summaries, retry_summaries, retry_rows = [], [], [], [], []
    for name in root['run_order']:
        directory = BATCH / name
        meta = json.loads((directory / 'metadata.json').read_text())
        assert meta['status'] == 'completed'
        assert helpers.sha256(directory / 'pm4py_alignment_timings.csv') == meta['csv_sha256']
        original = helpers.read_csv(directory / 'pm4py_alignment_timings.csv')
        reference = helpers.read_csv(PROJECT / meta['reference_csv'])
        manifest = json.loads((PROJECT / meta['sample_manifest']).read_text())
        assert len(original) == len(reference) == len(manifest) == 1000
        retries = {}
        if (directory / 'retry_timings.csv').exists():
            retry_meta = json.loads((directory / 'retry_metadata.json').read_text())
            assert retry_meta['status'] == 'completed'
            assert helpers.sha256(directory / 'retry_timings.csv') == retry_meta['csv_sha256']
            retries = {row['sample_index']: row for row in helpers.read_csv(directory / 'retry_timings.csv')}
            assert set(retries) == {row['sample_index'] for row in original if row['status'] == 'timeout'}
        selected = []
        for row, baseline, record in zip(original, reference, manifest):
            assert all(row[key] == baseline[key] == str(record[key]) for key in IDS)
            assert row['ebi_alignment_cost'] == baseline['alignment_cost']
            original_row = {'dataset': meta['dataset'], 'model': meta['model'], **row}
            originals.append(original_row)
            chosen = retries.get(row['sample_index'], row)
            assert all(chosen[key] == row[key] for key in IDS)
            assert chosen['ebi_alignment_cost'] == baseline['alignment_cost']
            for attempt in (row, chosen):
                if attempt['elapsed_ns']:
                    assert Decimal(attempt['elapsed_ms']) * 1_000_000 == int(attempt['elapsed_ns'])
                if attempt['status'] == 'ok':
                    assert attempt['alignment_cost'] == baseline['alignment_cost']
                    assert attempt['alignment_checks_passed'] == '1'
                    assert int(attempt['log_moves']) + int(attempt['visible_model_moves']) == int(attempt['alignment_cost'])
            final_row = {'dataset': meta['dataset'], 'model': meta['model'], **chosen,
                         'measurement_attempt': 2 if row['sample_index'] in retries else 1,
                         'timeout_seconds': RETRY_TIMEOUT if row['sample_index'] in retries else TIMEOUT,
                         'original_status': row['status'], 'original_attempt_elapsed_ns': row['attempt_elapsed_ns']}
            combined.append(final_row)
            selected.append(final_row)
            if row['sample_index'] in retries:
                retry_rows.append({'dataset': meta['dataset'], 'model': meta['model'], **chosen})
        summaries.append({'dataset': meta['dataset'], 'model': meta['model'], **summarize(original)})
        retry_summaries.append({'dataset': meta['dataset'], 'model': meta['model'],
                                'retried_traces': len(retries), **summarize(selected)})
    assert len(combined) == len({(r['dataset'], r['model'], r['trace_sha256']) for r in combined}) == 8000
    for path, digest in root['protected_files_sha256'].items():
        assert helpers.sha256(Path(path)) == digest, f'Source/input changed during benchmark: {path}'
    helpers.write_csv(BATCH / 'all_trace_timings.csv', ALL_FIELDS, originals)
    helpers.write_csv(BATCH / 'all_datasets_summary.csv', list(summaries[0]), summaries)
    extra = ['measurement_attempt', 'timeout_seconds', 'original_status', 'original_attempt_elapsed_ns']
    helpers.write_csv(BATCH / 'all_trace_timings_with_retries.csv', ALL_FIELDS + extra, combined)
    helpers.write_csv(BATCH / 'all_datasets_summary_with_retries.csv', list(retry_summaries[0]), retry_summaries)
    helpers.write_csv(BATCH / 'retry_timings.csv', ALL_FIELDS, retry_rows)
    helpers.write_csv(BATCH / 'initial_failures_and_cost_mismatches.csv', ALL_FIELDS, [r for r in originals if r['status'] != 'ok'])
    helpers.write_csv(BATCH / 'final_failures_and_cost_mismatches.csv', ALL_FIELDS + extra, [r for r in combined if r['status'] != 'ok'])
    assert helpers.read_csv(BATCH / 'all_trace_timings_with_retries.csv') == [{k: str(v) for k, v in r.items()} for r in combined]
    root.update(status='completed', finished_utc=helpers.now(), original_summary=summarize(originals),
                final_summary=summarize(combined), retried_traces=len(retry_rows), input_and_source_hashes_verified_after_run=True,
                outputs_sha256={p.name: helpers.sha256(p) for p in BATCH.glob('*.csv')})
    helpers.save(BATCH / 'metadata.json', root)
    lines = ['# Corrected PM4Py process-tree alignment benchmark', '',
             f'PM4Py **{root["pm4py_version"]}**, variant **{root["default_variant"]}**, using '
             '`pm4py.algo.conformance.alignments.process_tree.algorithm.apply(Trace, ProcessTree)` '
             'with its default parameters and the local correctness fixes described in [PATCH.md](PATCH.md).', '',
             '## Conditions', '',
             '- The exact saved seed-42 samples from `rerun_02`: 1,000 unique ordered activity sequences '
             'for each year, reused in the same order for pt00, pt10, pt25 and pt50. '
             'All lifecycle events are retained; activity attribute is `concept:name`.',
             '- Original XES, saved samples, manifests, PTML and reference timing files are verified by SHA-256. '
             'Sample case IDs, activity sequences, lengths and sequence hashes are checked. '
             'PM4Py tree imports are checked against source PTML structure.',
             '- One fresh Python process per model; models run sequentially. Three warm-ups on samples '
             '1–3, then one timed call per selected trace. The interpreter, dependencies and '
             '`PYTHONHASHSEED=42` match the previous Python benchmark; normal garbage collection stays enabled.',
             '- `time.perf_counter_ns` measures the public alignment call, including leaf extraction, '
             'activity projection, search, reconstruction and internal cleanup. Imports, XES/PTML parsing, '
             'input Trace construction, timer alarm setup/reset, result validation/destruction and CSV output '
             'are outside the timer. Each call creates a new search; no cross-trace alignment cache.',
             '- The input is a PM4Py `Trace` containing the saved ordered activity names and representative '
             'case ID. Other event attributes are not needed by this algorithm. The original tree is used '
             'directly, with corrections to silent-move reconstruction and priority-queue maintenance.',
             '- Costs are 1 per visible log/model move and 0 per silent/synchronous move. '
             'Each cost is compared with Ebi special run 03. Returned moves are checked for exact trace '
             'consumption, matching synchronous labels, valid model leaf references and consistent cost. '
             'The returned search state must consume the trace and close the root. These checks do not '
             'independently replay every tree-state transition.',
             '- PM4Py returns a full alignment; the previous Python `dyn_align` returns only a cost. '
             'This difference in work must be considered when comparing timings.',
             f'- Initial per-trace limit: {TIMEOUT} seconds. After the eight original runs, timed-out '
             f'traces are retried once with {RETRY_TIMEOUT} seconds, in fresh processes with three warm-ups. '
             'Original outcomes remain preserved. Retry timings have a different preceding call history '
             'and are marked with `measurement_attempt=2`.',
             '- Failed/timeout calls have blank completed timings and costs; `attempt_elapsed_ns` records '
             'the attempted duration. Completed timing statistics include any completed validation errors '
             'or cost mismatches, with separate status counts. Percentiles use inclusive interpolation.', '',
             '## Results including separate retries', '',
             '| Dataset | Model | Completed | Fully checked | Retries | Mean (ms) | Median (ms) | Max (ms) |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    def fmt(value):
        return f'{value:.3f}' if value != '' else '—'
    for row in retry_summaries:
        lines.append(f'| {row["dataset"]} | {row["model"]} | {row["completed_traces"]} | '
                     f'{row["fully_validated_alignments"]} | {row["retried_traces"]} | {fmt(row["mean_ms"])} | '
                     f'{fmt(row["median_ms"])} | {fmt(row["max_ms"])} |')
    final = root['final_summary']
    lines += ['', f'Initial outcomes: {root["original_summary"]["completed_traces"]}/8000 completed, '
              f'{root["original_summary"]["timeouts"]} timeouts. Final outcomes: '
              f'{final["completed_traces"]}/8000 completed, {final["cost_mismatches"]} cost mismatches, '
              f'{final["validation_errors"]} validation errors, {final["errors"]} execution errors, '
              f'{final["timeouts"]} remaining timeouts. Total completed-call time: '
              f'**{final["sum_alignment_seconds"]:.3f} seconds**.', '', '## Files', '',
              '- [All per-trace timings including retries](all_trace_timings_with_retries.csv).',
              '- [Summary including retries](all_datasets_summary_with_retries.csv).',
              '- [Original timings](all_trace_timings.csv) and [original summary](all_datasets_summary.csv).',
              '- [Original issues](initial_failures_and_cost_mismatches.csv) and [remaining issues](final_failures_and_cost_mismatches.csv).',
              '- [Retry timings](retry_timings.csv), [metadata](metadata.json), [environment](environment.txt).',
              '- Each model directory contains original timings, metadata, imported tree and progress log; '
              'retry files are added when needed.',
              '- [Runner](run_benchmark.py); `benchmark_helpers.py` is an unchanged copy of the previous '
              'Python runner, reused only for input checks, CSV/JSON helpers and summary statistics. '
              '`source_snapshot` contains the installed PM4Py alignment, semantics and supporting source files.', '']
    (BATCH / 'README.md').write_text('\n'.join(lines))
    print(json.dumps(root['final_summary']), flush=True)


def run_all():
    assert not (BATCH / 'metadata.json').exists(), 'Refusing to overwrite an existing run'
    previous = json.loads((PREVIOUS / 'metadata.json').read_text())
    python = previous['python_executable']
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='42',
               MPLCONFIGDIR='/tmp/ebi_pm4py_benchmark_matplotlib')
    subprocess.run([python, '-B', str(Path(__file__).resolve()), '--prepare'], env=env, check=True)
    root = json.loads((BATCH / 'metadata.json').read_text())
    root.update(status='running', started_utc=helpers.now())
    helpers.save(BATCH / 'metadata.json', root)

    def execute(name, retry=False):
        command = [python, '-B', str(Path(__file__).resolve()), '--worker', name]
        if retry:
            command.append('--retry')
        print(f'{"Retrying" if retry else "Starting"} {name}', flush=True)
        with (BATCH / name / ('retry.log' if retry else 'run.log')).open('x') as stream:
            result = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            root.update(status='failed', failed_model=name, returncode=result.returncode)
            helpers.save(BATCH / 'metadata.json', root)
            raise RuntimeError(f'{name} worker failed; see its log')
        child = json.loads((BATCH / name / ('retry_metadata.json' if retry else 'metadata.json')).read_text())
        root.setdefault('retry_runs' if retry else 'completed_runs', []).append({'name': name, 'command': command, 'finished_utc': helpers.now()})
        helpers.save(BATCH / 'metadata.json', root)
        print(f'Finished {name}: {json.dumps(child["summary"])}', flush=True)

    for name in root['run_order']:
        execute(name)
    for name in root['run_order']:
        rows = helpers.read_csv(BATCH / name / 'pm4py_alignment_timings.csv')
        if any(row['status'] == 'timeout' for row in rows):
            execute(name, retry=True)
    report()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--worker')
    parser.add_argument('--retry', action='store_true')
    parser.add_argument('--report', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    elif args.worker:
        worker(args.worker, args.retry)
    elif args.report:
        report()
    else:
        run_all()
