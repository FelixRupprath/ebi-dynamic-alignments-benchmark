"""Run the standard command routine sequentially, preserving censored attempts."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import queue
import resource
import statistics
import subprocess
import threading
import time

BATCH = Path(__file__).resolve().parent
PROJECT = BATCH.parent.parent
TIMEOUT = 60
ADDRESS_SPACE_BYTES = 8 * 1024**3
IDS = ['sample_index', 'source_trace_index', 'case_id', 'variant_frequency', 'trace_length', 'trace_sha256']
FIELDS = IDS + ['elapsed_ns', 'elapsed_ms', 'alignment_cost', 'weighted_alignment_cost',
                'silent_moves', 'synchronous_moves', 'alignment_moves', 'reference_alignment_cost',
                'cost_matches_reference', 'replay_valid', 'status', 'timeout_limit_seconds',
                'observed_wall_ns', 'error']


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def set_limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_AS, (ADDRESS_SPACE_BYTES, ADDRESS_SPACE_BYTES))


class Worker:
    def __init__(self, directory, meta):
        self.directory = directory
        self.meta = meta
        self.process = None
        self.starts = 0

    def start(self):
        self.frames = queue.Queue()
        self.diagnostics = (self.directory / 'diagnostics.log').open('a')
        self.diagnostic_offset = self.diagnostics.tell()
        command = [str(BATCH / 'standard-alignments-benchmark'),
                   str(self.directory / 'sampled_1000_unique_traces.xes'),
                   str(PROJECT / self.meta['model_ptml']), str(self.directory / 'sample_manifest.json')]
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=self.diagnostics, text=True, bufsize=1,
                                        cwd=PROJECT, preexec_fn=set_limits,
                                        env=dict(os.environ, RAYON_NUM_THREADS='1'))
        self.starts += 1
        process, frames = self.process, self.frames

        def reader():
            try:
                for line in process.stdout:
                    frames.put(json.loads(line))
            except Exception as error:
                frames.put({'phase': 'protocol_error', 'error': str(error)})
            finally:
                frames.put({'phase': 'exited'})

        self.reader = threading.Thread(target=reader, daemon=True)
        self.reader.start()
        frame = self.frames.get(timeout=30)
        if frame.get('phase') != 'ready' or frame.get('traces') != len(json.loads((self.directory / 'sample_manifest.json').read_text())):
            self.stop()
            raise RuntimeError(f'Worker import failed: {frame}; see diagnostics.log')

    def stop(self, graceful=False):
        if self.process is None:
            return
        if self.process.poll() is None:
            if graceful:
                self.process.stdin.close()
                try:
                    self.process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    self.process.kill()
            else:
                self.process.kill()
        self.process.wait()
        self.reader.join(timeout=5)
        self.diagnostics.close()
        self.process = None

    def run(self, index):
        if self.process is None:
            self.start()
        self.process.stdin.write(json.dumps({'sample_index': index}) + '\n')
        self.process.stdin.flush()
        frame = self.frames.get(timeout=30)
        if frame.get('phase') != 'started':
            self.stop()
            raise RuntimeError(f'Worker did not start trace {index}: {frame}')
        observed_start = time.monotonic_ns()
        elapsed_ns = None
        try:
            frame = self.frames.get(timeout=TIMEOUT)
            observed_ns = time.monotonic_ns() - observed_start
            if frame.get('phase') == 'computed':
                elapsed_ns = frame['elapsed_ns']
                frame = self.frames.get(timeout=30)
            if frame.get('phase') == 'result':
                assert frame['sample_index'] == index
                return {**frame, 'observed_wall_ns': observed_ns}
            self.stop()
            diagnostic = (self.directory / 'diagnostics.log').read_text()[self.diagnostic_offset:][-2000:]
            status = 'memory_limit' if 'memory allocation of' in diagnostic else 'worker_error'
            return {'status': status, 'elapsed_ns': elapsed_ns, 'observed_wall_ns': observed_ns, 'error': diagnostic}
        except queue.Empty:
            observed_ns = time.monotonic_ns() - observed_start
            self.stop()
            return {'status': 'timeout' if elapsed_ns is None else 'validation_timeout',
                    'elapsed_ns': elapsed_ns, 'observed_wall_ns': observed_ns,
                    'error': 'Alignment exceeded 60 seconds' if elapsed_ns is None else 'Replay validation did not finish within 30 seconds'}


def summary(rows):
    completed = [row for row in rows if row['status'] in ('ok', 'cost_mismatch')]
    times = [int(row['elapsed_ns']) / 1_000_000 for row in completed]
    return {
        'attempted_traces': len(rows), 'completed_traces': len(completed),
        'mean_ms': statistics.mean(times) if times else None,
        'median_ms': statistics.median(times) if times else None,
        'sum_alignment_seconds': sum(int(row['elapsed_ns']) for row in completed) / 1_000_000_000,
        'min_ms': min(times) if times else None, 'max_ms': max(times) if times else None,
        'p95_ms': statistics.quantiles(times, n=100, method='inclusive')[94] if len(times) > 1 else (times[0] if times else None),
        'zero_cost_traces': sum(int(row['alignment_cost']) == 0 for row in completed),
        'cost_matches_reference': sum(row['status'] == 'ok' for row in rows),
        'cost_mismatches': sum(row['status'] == 'cost_mismatch' for row in rows),
        'timeouts': sum(row['status'] == 'timeout' for row in rows),
        'memory_limit_failures': sum(row['status'] == 'memory_limit' for row in rows),
        'other_failures': sum(row['status'] not in ('ok', 'cost_mismatch', 'timeout', 'memory_limit') for row in rows),
        'observed_attempt_seconds': sum(int(row['observed_wall_ns']) for row in rows) / 1_000_000_000,
    }


def main():
    meta_path = BATCH / 'metadata.json'
    meta = json.loads(meta_path.read_text())
    assert meta['status'] == 'preparing', 'Batch already started; existing results are protected'
    for path, expected in meta['source_sha256'].items():
        assert digest(PROJECT / path) == expected
    for path, expected in meta['original_results_sha256'].items():
        assert digest(PROJECT / path) == expected
    names = [f'bpi{year}_pt{level}_seed42_1000' for year in (2012, 2017) for level in ('00', '10', '25', '50')]
    meta.update(status='running', started_utc=now(), run_order=names, completed_runs=[],
                timeout_seconds=TIMEOUT, address_space_limit_bytes=ADDRESS_SPACE_BYTES,
                timeout_note='Same 60-second per-trace cutoff as the preceding C++ test; enforced externally because the standard routine has no cancellation API.',
                memory_note='8 GiB virtual-address-space limit per worker; limit failures are explicit censored rows, not completed alignment times.',
                restart_note='A failed worker is replaced before the next trace. The first three sample traces are attempted as warm-ups once per model; no additional alignment warm-ups follow a worker restart.',
                benchmark_runner_sha256=digest(BATCH / 'benchmark_runner.rs'),
                controller_sha256=digest(Path(__file__)), executable_sha256=digest(BATCH / 'standard-alignments-benchmark'),
                cpu_models=sorted({line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')}),
                logical_cpu_count=os.cpu_count(),
                cost_model={'log_move': 10000, 'visible_model_move': 10000, 'silent_move': 1, 'synchronous_move': 0},
                cost_columns='alignment_cost counts visible deviations; weighted_alignment_cost is 10000 * alignment_cost + silent_moves, computed from the replayed alignment.',
                timing_includes='One align_stochastic_language call on a one-trace language: activity translation, A* search, move reconstruction, probability/result assembly and internal cleanup.',
                timing_excludes='Startup, XES/PTML import, PTML-to-Petri-net conversion, one-trace language creation, replay validation, serialization/output and dropping the returned result.',
                summary_policy='Timing statistics include only completed replay-valid alignments. Timed-out/failed measurements are censored and have separate counts; missing times are never replaced by zero or the timeout threshold.')
    save(meta_path, meta)
    try:
        for name in names:
            directory = BATCH / name
            child_path = directory / 'metadata.json'
            child = json.loads(child_path.read_text())
            assert child['status'] == 'prepared'
            for path, expected in [(PROJECT / child['source_xes'], child['source_xes_sha256']),
                                   (PROJECT / child['model_ptml'], child['model_ptml_sha256']),
                                   (directory / 'sampled_1000_unique_traces.xes', child['sample_xes_sha256']),
                                   (directory / 'sample_manifest.json', child['sample_manifest_sha256'])]:
                assert digest(path) == expected
            manifest = json.loads((directory / 'sample_manifest.json').read_text())
            reference = read_csv(PROJECT / child['reference_directory'] / 'special_alignment_timings.csv')
            assert len(manifest) == len(reference) == 1000
            assert len({record['trace_sha256'] for record in manifest}) == 1000
            for record, row in zip(manifest, reference):
                assert all(str(record[field] if record[field] is not None else '') == row[field] for field in IDS)
            print(f'\nStarting {name}', flush=True)
            child.update(status='running', started_utc=now(), warmup_results=[])
            save(child_path, child)
            worker = Worker(directory, child)
            try:
                for index in (1, 2, 3):
                    result = worker.run(index)
                    child['warmup_results'].append(result)
                    save(child_path, child)
                    print(f'{name}: warm-up {index}/3 {result["status"]}', flush=True)
                with (directory / 'alignment_timings.csv').open('x', newline='') as stream, (directory / 'run.log').open('w') as log:
                    writer = csv.DictWriter(stream, fieldnames=FIELDS)
                    writer.writeheader()
                    failures = mismatches = 0
                    for index, (record, old) in enumerate(zip(manifest, reference), 1):
                        result = worker.run(index)
                        expected = int(old['alignment_cost'])
                        status = result['status']
                        if status == 'ok' and result['alignment_cost'] != expected:
                            status = 'cost_mismatch'
                            mismatches += 1
                        failures += int(status not in ('ok', 'cost_mismatch'))
                        elapsed = result.get('elapsed_ns')
                        row = {**{field: record[field] for field in IDS},
                               **{field: result.get(field, '') for field in FIELDS if field not in IDS},
                               'elapsed_ns': elapsed if elapsed is not None else '',
                               'elapsed_ms': f'{elapsed / 1_000_000:.6f}' if elapsed is not None else '',
                               'reference_alignment_cost': expected, 'status': status,
                               'cost_matches_reference': int(status == 'ok'),
                               'replay_valid': int(result.get('replay_valid', False)), 'timeout_limit_seconds': TIMEOUT}
                        writer.writerow(row)
                        stream.flush()
                        if index % 25 == 0 or status not in ('ok', 'cost_mismatch'):
                            line = f'{name}: {index}/1000 attempted; {failures} incomplete, {mismatches} cost differences; last {status}'
                            print(line, flush=True)
                            log.write(line + '\n')
                            log.flush()
            finally:
                worker.stop(graceful=True)
            path = directory / 'alignment_timings.csv'
            rows = read_csv(path)
            assert len(rows) == 1000
            child.update(status='completed', finished_utc=now(), worker_starts=worker.starts, csv_sha256=digest(path), summary=summary(rows))
            save(child_path, child)
            meta['completed_runs'].append(name)
            save(meta_path, meta)
            print(json.dumps({'model': name, **child['summary']}), flush=True)
    except Exception as error:
        meta.update(status='failed', error=str(error))
        save(meta_path, meta)
        raise
    meta.update(status='benchmarks_completed', finished_utc=now())
    save(meta_path, meta)
    print('All eight standard-alignments benchmarks attempted.', flush=True)


if __name__ == '__main__':
    main()
