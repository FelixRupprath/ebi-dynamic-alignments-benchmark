# Portable adapter of the original five-way protocol. The native DP runner
# from native_dynamic_best_of_three_run_03 replaces historical Special Alignments.
# Other timed interfaces, warmups and minimum selection are preserved.
"""Repeat the five preserved alignment implementations three times per input pair."""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import gc
import hashlib
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

ROOT = Path(os.environ['EBI_REPRO_RUN_DIR']).resolve()
PROJECT = Path(__file__).resolve().parent.parent
ALGORITHMS = ['ebi_native_dynamic', 'cpp', 'ebi_default', 'python_dynamic', 'pm4py']
LABELS = dict(ebi_native_dynamic='Ebi Native Dynamic Alignments', cpp='C++ Dynamic Alignments',
              ebi_default='Ebi Default Alignments', python_dynamic='Python Dynamic Alignments',
              pm4py='PM4Py process-tree alignments (corrected)')
IDS = ['sample_index', 'source_trace_index', 'case_id', 'variant_frequency', 'trace_length', 'trace_sha256']
DETAILS = ['elapsed_ns', 'elapsed_ms', 'alignment_cost', 'reference_alignment_cost',
           'cost_matches_reference', 'alignment_checks_passed', 'status', 'attempt_elapsed_ns', 'error']
FIELDS = IDS + DETAILS
KEYS = ['dataset', 'model_dataset', 'model', 'algorithm', 'repetition']
LEVELS = ['00', '10', '25', '50']


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fields, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def pair_name(year, model_year, level):
    return f'bpi{year}_model{model_year}_pt{level}'


def reference_rows(meta, pair):
    if pair['dataset'] == pair['model_dataset']:
        path = ROOT / f'inputs/reference/bpi{pair["dataset"]}_{pair["model"]}.csv'
    else:
        path = ROOT / 'raw/repetition_01/ebi_native_dynamic' / pair['name'] / 'timings.csv'
    return read_csv(path) if path.exists() else None


def worker(algorithm, year, model_year, level, repetition):
    sys.path.insert(0, str(ROOT / 'support'))
    import benchmark_helpers as helpers
    meta = json.loads((ROOT / 'metadata.json').read_text())
    pair = next(p for p in meta['pairs'] if p['name'] == pair_name(year, model_year, level))
    directory = ROOT / f'raw/repetition_{repetition:02d}' / algorithm / pair['name']
    directory.mkdir(parents=True, exist_ok=False)
    sample = ROOT / meta['samples'][str(year)]['sample_xes']
    manifest_path = ROOT / meta['samples'][str(year)]['manifest']
    model_path = ROOT / meta['models'][f'{model_year}_pt{level}']
    manifest = json.loads(manifest_path.read_text())
    traces = helpers.read_sample(sample, manifest)
    reference = reference_rows(meta, pair)
    if reference is not None:
        assert len(reference) == meta['samples_per_dataset']
        for record, row in zip(manifest, reference):
            assert all(str(record[k]) == row[k] for k in IDS)
    run_meta = dict(status='running', started_utc=now(), algorithm=algorithm, repetition=repetition, **pair,
                    python=sys.version, executable=sys.executable, cpu_affinity=sorted(os.sched_getaffinity(0)),
                    garbage_collection_enabled=gc.isenabled(), gc_threshold=gc.get_threshold(), warmups=[])
    save(directory / 'metadata.json', run_meta)
    if algorithm == 'ebi_native_dynamic':
        native_csv = directory / 'native_timings.csv'
        runner = 'ebi-native-dynamic-smoke' if meta.get('smoke_test') else 'ebi-native-dynamic-benchmark'
        costs = ROOT / f'inputs/reference/bpi{year}_pt{level}_costs.json'
        command = [str(ROOT / 'bin' / runner), str(sample), str(model_path), str(manifest_path), str(costs), str(native_csv)]
        subprocess.run(command, check=True)
        native = read_csv(native_csv)
        assert len(native) == meta['samples_per_dataset']
        rows = []
        for record, result in zip(manifest, native):
            assert all(str(record[k]) == result[k] for k in IDS)
            expected = int(reference[len(rows)]['alignment_cost']) if reference else int(result['alignment_cost'])
            cost = int(result['alignment_cost'])
            rows.append(dict(**{k: record[k] for k in IDS}, elapsed_ns=result['elapsed_ns'], elapsed_ms=result['elapsed_ms'],
                             alignment_cost=cost, reference_alignment_cost=expected, cost_matches_reference=int(cost == expected),
                             alignment_checks_passed=int(result['replay_valid'] == 'true'),
                             status=result['status'] if cost == expected and result['cost_matches_reference'] == 'true' and result['replay_valid'] == 'true' else 'validation_failed',
                             attempt_elapsed_ns=result['elapsed_ns'], error=result['error']))
        write_csv(directory / 'timings.csv', FIELDS, rows)
        run_meta['warmups'] = 'Samples 1, 2, 3, performed and replay-validated by preserved native runner'
    else:
        native_worker = None
        if algorithm == 'ebi_default':
            control = module('default_controller', ROOT / 'support/default_controller.py')
            control.BATCH, control.PROJECT, control.TIMEOUT = ROOT / 'bin', PROJECT, 600
            for name, path in [('sampled_1000_unique_traces.xes', sample), ('sample_manifest.json', manifest_path)]:
                (directory / name).symlink_to(os.path.relpath(path, directory))
            native_worker = control.Worker(directory, {'model_ptml': str(model_path)})

            def measure(index):
                result = native_worker.run(index + 1)
                result['alignment_checks_passed'] = int(result.get('replay_valid', False))
                result['attempt_elapsed_ns'] = result.get('elapsed_ns') or result.get('observed_wall_ns', '')
                return result
        else:
            import pm4py
            tree = pm4py.read_ptml(str(model_path))
            helpers.check_model_import(model_path, tree)
            (directory / 'model_tree_string.txt').write_text(str(tree) + '\n')
            run_meta['pm4py_version'] = pm4py.__version__
            if algorithm == 'cpp':
                alignment = module('alignment', Path(meta['cpp_module']))
                aligner = alignment.AlignmentWrapper()
                aligner.loadTree(str(tree))
                traces = [list(trace) for trace in traces]

                def measure(index):
                    started = time.perf_counter_ns()
                    try:
                        cost = aligner.align(traces[index])
                        elapsed = time.perf_counter_ns() - started
                    except Exception as error:
                        return dict(status='error', attempt_elapsed_ns=time.perf_counter_ns() - started, error=str(error))
                    return dict(status='ok' if cost >= 0 else 'timeout', elapsed_ns=elapsed if cost >= 0 else '',
                                alignment_cost=cost if cost >= 0 else '', attempt_elapsed_ns=elapsed)
            elif algorithm == 'python_dynamic':
                implementation = module('benchmarked_dynamic_alignment', ROOT / 'sources/dynamic_alignment.py')
                helpers.TIMEOUT = 600
                signal.signal(signal.SIGALRM, helpers.alarm_handler)

                def measure(index):
                    return helpers.measure(implementation.dyn_align, traces[index], tree, int(reference[index]['alignment_cost']))
            else:
                from pm4py.algo.conformance.alignments.process_tree import algorithm as implementation
                from pm4py.objects.log.obj import Event, Trace
                pm_runner = module('pm4py_runner', ROOT / 'support/pm4py_runner.py')
                assert implementation.DEFAULT_VARIANT.name == 'SEARCH_GRAPH_PT'
                assert digest(implementation.DEFAULT_VARIANT.value.__file__) == meta['pm4py_local_patch']['patched_sha256']
                traces = [Trace([Event({'concept:name': activity}) for activity in trace],
                                attributes={'concept:name': record['case_id']}) for trace, record in zip(traces, manifest)]
                signal.signal(signal.SIGALRM, pm_runner.alarm_handler)

                def measure(index):
                    return pm_runner.measure(implementation.apply, traces[index], tree, int(reference[index]['alignment_cost']), 600)
        try:
            for index in range(3):
                result = measure(index)
                assert result['status'] == 'ok', ('Warm-up failed', algorithm, pair, index, result)
                run_meta['warmups'].append({'sample_index': index + 1, **result})
            save(directory / 'metadata.json', run_meta)
            with (directory / 'timings.csv').open('x', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS)
                writer.writeheader()
                for index, record in enumerate(manifest):
                    result = measure(index)
                    expected = int(reference[index]['alignment_cost'])
                    cost = result.get('alignment_cost', '')
                    status = result['status']
                    if status == 'ok' and cost != expected:
                        status = 'cost_mismatch'
                    elapsed = result.get('elapsed_ns') or ''
                    row = {k: result.get(k, '') for k in DETAILS}
                    row.update({k: record[k] for k in IDS})
                    row.update(elapsed_ns=elapsed, elapsed_ms=str(Decimal(elapsed) / 1_000_000) if elapsed else '',
                               alignment_cost=cost, reference_alignment_cost=expected,
                               cost_matches_reference=int(cost == expected) if cost != '' else '', status=status)
                    writer.writerow(row)
                    stream.flush()
                    if (index + 1) % 25 == 0 or status != 'ok':
                        print(f'{algorithm} {pair["name"]} repetition {repetition}: {index + 1}/{meta["samples_per_dataset"]}; {status}; {row["elapsed_ms"]} ms', flush=True)
        finally:
            if native_worker is not None:
                native_worker.stop(graceful=True)
    rows = read_csv(directory / 'timings.csv')
    assert len(rows) == meta['samples_per_dataset']
    run_meta.update(status='completed', finished_utc=now(), csv_sha256=digest(directory / 'timings.csv'),
                    completed=sum(r['status'] == 'ok' for r in rows), failures=sum(r['status'] != 'ok' for r in rows),
                    sum_alignment_seconds=sum(int(r['elapsed_ns']) for r in rows if r['elapsed_ns']) / 1e9)
    save(directory / 'metadata.json', run_meta)
    print(json.dumps({k: run_meta[k] for k in ['algorithm', 'name', 'repetition', 'completed', 'failures', 'sum_alignment_seconds']}), flush=True)


def report():
    meta = json.loads((ROOT / 'metadata.json').read_text())
    repetitions, minima, wide, summaries = [], [], [], []
    failures = []
    for pair in meta['pairs']:
        manifest = json.loads((ROOT / meta['samples'][str(pair['dataset'])]['manifest']).read_text())
        collected = {}
        for algorithm in ALGORITHMS:
            runs = []
            for repetition in [1, 2, 3]:
                directory = ROOT / f'raw/repetition_{repetition:02d}' / algorithm / pair['name']
                run_meta = json.loads((directory / 'metadata.json').read_text())
                assert run_meta['status'] == 'completed' and digest(directory / 'timings.csv') == run_meta['csv_sha256']
                rows = read_csv(directory / 'timings.csv')
                assert len(rows) == meta['samples_per_dataset']
                for record, row in zip(manifest, rows):
                    assert all(str(record[k]) == row[k] for k in IDS)
                    if row['elapsed_ns']:
                        assert int(row['elapsed_ns']) > 0 and Decimal(row['elapsed_ms']) * 1_000_000 == int(row['elapsed_ns'])
                    enriched = dict(dataset=pair['dataset'], model_dataset=pair['model_dataset'], model=pair['model'],
                                    algorithm=algorithm, repetition=repetition, **row)
                    repetitions.append(enriched)
                    if row['status'] != 'ok':
                        failures.append(enriched)
                runs.append(rows)
            selected = []
            for index, record in enumerate(manifest):
                candidates = [(rep + 1, rows[index]) for rep, rows in enumerate(runs) if rows[index]['status'] == 'ok']
                assert len(candidates) == 3, ('Three successful repetitions required', pair, algorithm, index + 1)
                assert len({r['alignment_cost'] for _, r in candidates}) == 1
                repetition, best = min(candidates, key=lambda item: (int(item[1]['elapsed_ns']), item[0]))
                chosen = dict(dataset=pair['dataset'], model_dataset=pair['model_dataset'], model=pair['model'],
                              algorithm=algorithm, repetition=repetition, **best)
                for rep, result in candidates:
                    chosen[f'repetition_{rep}_elapsed_ns'] = result['elapsed_ns']
                selected.append(chosen)
                minima.append(chosen)
            collected[algorithm] = selected
            times = [int(r['elapsed_ns']) for r in selected]
            summaries.append(dict(dataset=pair['dataset'], model_dataset=pair['model_dataset'], model=pair['model'], algorithm=algorithm,
                                  traces=meta['samples_per_dataset'], repetitions_per_trace=3, minimum_ms=min(times) / 1e6,
                                  mean_ms=statistics.mean(times) / 1e6, median_ms=statistics.median(times) / 1e6,
                                  p95_ms=statistics.quantiles(times, n=100, method='inclusive')[94] / 1e6,
                                  maximum_ms=max(times) / 1e6, sum_minimum_seconds=sum(times) / 1e9,
                                  sum_all_repetitions_seconds=sum(int(r['elapsed_ns']) for rows in runs for r in rows) / 1e9))
        for index, record in enumerate(manifest):
            row = dict(dataset=pair['dataset'], model_dataset=pair['model_dataset'], model=pair['model'], **{k: record[k] for k in IDS})
            costs = {collected[a][index]['alignment_cost'] for a in ALGORITHMS}
            assert len(costs) == 1, ('Algorithms disagree', pair, index + 1, costs)
            row['alignment_cost'] = next(iter(costs))
            for algorithm in ALGORITHMS:
                chosen = collected[algorithm][index]
                for field in ['elapsed_ns', 'elapsed_ms', 'repetition_1_elapsed_ns', 'repetition_2_elapsed_ns', 'repetition_3_elapsed_ns']:
                    row[f'{algorithm}_{field}'] = chosen[field]
                row[f'{algorithm}_selected_repetition'] = chosen['repetition']
            wide.append(row)
    assert len(repetitions) == meta['expected_timed_alignments']
    assert len({tuple(str(r[k]) for k in KEYS + ['sample_index']) for r in repetitions}) == len(repetitions)
    for path, expected in meta['protected_files_sha256'].items():
        protected = Path(path) if Path(path).is_absolute() else PROJECT / path
        assert digest(protected) == expected, f'Protected file changed: {path}'
    overall = []
    for algorithm in ALGORITHMS:
        times = [int(r['elapsed_ns']) for r in minima if r['algorithm'] == algorithm]
        overall.append(dict(algorithm=algorithm, traces=len(times), repetitions_per_trace=3,
                            mean_ms=statistics.mean(times) / 1e6, median_ms=statistics.median(times) / 1e6,
                            p95_ms=statistics.quantiles(times, n=100, method='inclusive')[94] / 1e6,
                            maximum_ms=max(times) / 1e6, sum_minimum_seconds=sum(times) / 1e9))
    for filename, rows, fields in [
        ('all_trace_repetitions.csv', repetitions, KEYS + FIELDS),
        ('minimum_trace_timings.csv', minima, list(minima[0])),
        ('all_programs_trace_timings.csv', wide, list(wide[0])),
        ('all_datasets_summary.csv', summaries, list(summaries[0])),
        ('overall_summary.csv', overall, list(overall[0])),
        ('failures.csv', failures, KEYS + FIELDS),
    ]:
        write_csv(ROOT / filename, fields, rows)
    verification = dict(verified_utc=now(), timed_alignments=len(repetitions), minimum_rows=len(minima),
                        comparison_rows=len(wide), dataset_algorithm_summary_rows=len(summaries),
                        all_three_repetitions_successful=True, all_costs_match=True, all_selected_times_are_minima=True,
                        identifiers_match_original_samples=True, protected_inputs_sources_and_previous_final_unchanged=True,
                        failures=len(failures), outputs_sha256={p.name: digest(p) for p in ROOT.glob('*.csv')})
    save(ROOT / 'verification.json', verification)
    meta.update(status='completed', finished_utc=now(), verification=verification)
    save(ROOT / 'metadata.json', meta)
    print(json.dumps(overall, indent=2), flush=True)


def run_all():
    import fcntl
    with (ROOT / 'run.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta = json.loads((ROOT / 'metadata.json').read_text())
        assert meta['status'] in ['prepared', 'running'], meta['status']
        meta.update(status='running', started_utc=meta.get('started_utc', now()))
        save(ROOT / 'metadata.json', meta)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='42',
                   MPLCONFIGDIR=str(ROOT / 'mpl_config'), RAYON_NUM_THREADS='1')
        for repetition in [1, 2, 3]:
            for algorithm in ALGORITHMS:
                for pair in meta['pairs']:
                    job = f'repetition_{repetition:02d}/{algorithm}/{pair["name"]}'
                    if job in meta['completed_jobs']:
                        continue
                    output_dir = ROOT / 'logs' / f'repetition_{repetition:02d}' / algorithm
                    output_dir.mkdir(parents=True, exist_ok=True)
                    executable = meta['python_executables'].get(algorithm, sys.executable)
                    command = [executable, '-B', str(Path(__file__).resolve()), '--worker', algorithm,
                               str(pair['dataset']), str(pair['model_dataset']), pair['model'][2:], str(repetition)]
                    meta.update(current_job=job, current_job_started_utc=now())
                    save(ROOT / 'metadata.json', meta)
                    print(f'{now()} START {job}', flush=True)
                    with (output_dir / f'{pair["name"]}.log').open('x') as log:
                        outcome = subprocess.run(command, env=env, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT)
                    if outcome.returncode:
                        meta.update(status='failed', failed_job=job, returncode=outcome.returncode)
                        save(ROOT / 'metadata.json', meta)
                        raise RuntimeError(f'{job} failed; see log')
                    child = json.loads((ROOT / 'raw' / job / 'metadata.json').read_text())
                    assert child['status'] == 'completed'
                    meta['completed_jobs'].append(job)
                    save(ROOT / 'metadata.json', meta)
                    print(f'{now()} DONE {job}: {child["completed"]}/{meta["samples_per_dataset"]} valid, {child["sum_alignment_seconds"]:.3f} s measured; {len(meta["completed_jobs"])}/{3 * 5 * len(meta["pairs"])} jobs', flush=True)
        report()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', nargs=5, metavar=('ALGORITHM', 'YEAR', 'MODEL_YEAR', 'LEVEL', 'REPETITION'))
    parser.add_argument('--report', action='store_true')
    args = parser.parse_args()
    if args.worker:
        algorithm, year, model_year, level, repetition = args.worker
        worker(algorithm, int(year), int(model_year), level, int(repetition))
    elif args.report:
        report()
    else:
        run_all()
