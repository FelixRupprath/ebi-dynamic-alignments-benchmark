# Five-algorithm benchmark with downloaded dependencies

This folder contains ordinary scripts, model/input manifests, and the preserved
alignment binaries. It does **not** contain `benchmark.pyz`, Python/C++ environment
archives, or pre-created virtual environments. The installer downloads the pinned
Python packages from PyPI into two local virtual environments.

## Setup and execution

Use **CPython 3.14.7 on compatible Linux x86_64**, including `venv` and `ensurepip`.
The preserved native binaries were built with glibc 2.44 and require compatible
system libraries. The C++ alignment extension is already compiled; installing its
Python environment does not require rebuilding the benchmark's C++ implementation.
Internet access is needed for the initial dependency installation.

Run these commands from this folder:

```bash
# Download dependencies, apply PM4Py corrections, and verify both environments.
python3.14 install_dependencies.py

# Run the full benchmark using the two original XES datasets.
python3.14 run_benchmark.py --data-dir /path/to/DATAFILES
```

The data directory must contain:

```text
BPI_Challenge_2012.xes
BPI_Challenge_2017.xes
```

Alternatively, put those files directly in this folder and run
`python3.14 run_benchmark.py` without `--data-dir`.

The launcher verifies the source-log checksums and reconstructs exactly the
previously selected 1,000 unique traces per dataset. The resulting sampled XES
files must match the original samples byte for byte. No new sampling or model
discovery occurs. All eight original PTML models and saved reference costs are
already included.

For a quick test before starting the full experiment:

```bash
python3.14 run_benchmark.py --data-dir /path/to/DATAFILES --smoke
```

This performs 360 measured calls across all five algorithms and eight models,
using the first three selected traces and three repetitions, plus excluded warmups.
The full run performs **120,000 measured calls**.

## Files and folders

| Path | Purpose |
|---|---|
| `install_dependencies.py` | Downloads exact package versions, creates both environments, applies/checks PM4Py fixes, and tests the interfaces |
| `run_benchmark.py` | Checks inputs, rebuilds sampled XES logs, runs all algorithms, and verifies results |
| `environment/requirements_cpp.txt` | Exact pins for the 22-package C++ interface/import environment |
| `environment/requirements_python.txt` | Exact pins for the 46-package Python Dynamic/PM4Py environment |
| `environment/manifest.json` | Original platform/version inventories, requirement hashes and PM4Py correction hashes |
| `bin/` | Preserved native Dynamic, Default and C++ executables/module, plus the native smoke runner |
| `inputs/` | Two selection manifests, eight PTML models and reference costs; sampled XES logs are generated here |
| `sources/` | Python Dynamic implementation and the ten preserved PM4Py replacement modules |
| `support/`, `scripts/` | Measurement adapters, controller and independent result verifier |
| `reference_results/metadata.json` | Benchmark configuration template, required by the controller |
| `package_manifest.json` | Checksums of the distributed runtime files and expected sampled XES checksums |
| `licenses/` | Retained Ebi and pybind11 license texts |

Generated files are **not needed when distributing this folder**:

- `.venvs/`: installed dependencies, created by the installer.
- `.cache/` and `__pycache__/`: generated caches.
- `runs/`: new benchmark measurements and logs.
- `environment/runtime.json` and `environment/install_*.json`: readiness and pip
  artifact-download reports; these are recreated by installation/verification.
- `sample_reconstruction.json` and `inputs/bpi*/sampled_1000_unique_traces.xes`:
  generated from the supplied full datasets.

## Installation details

The environments retain the original complete inventories, including PM4Py
**2.7.23.8**. Scientific dependencies are installed from compatible wheels using
exact version pins, without resolving additional dependencies. The C++ environment
also historically contained the unrelated pure-Python package `alignment==1.0.10`;
PyPI supplies it as a source distribution, so pip builds its small Python wheel
in an isolated build environment. This package is not the benchmark's C++ solver:
the controller loads the preserved `bin/alignment*.so` explicitly.

The ten saved PM4Py modules are copied into the **Python environment only**, retaining
silent-move reconstruction and priority-queue corrections. The C++ environment's
PM4Py importer remains unpatched. Setup checks the exact installed versions,
`pip check`, binary/module hashes, and known alignment examples before marking
both environments ready. A partially installed environment is never reported ready.

Downloaded packages retain the historical version numbers, but downloaded wheel
artifacts can differ from the original installed environment bytes. Pip reports
in `environment/install_*.json` record artifact URLs and hashes. The alignment
binaries, Python Dynamic implementation, and PM4Py corrections remain the
preserved versions. Absolute timings can vary with machine, load, and builds.

Useful installer options:

```bash
# Check installed environments without accessing the network.
python3.14 install_dependencies.py --verify-only

# Install/recreate one environment at a time.
python3.14 install_dependencies.py --kind cpp
python3.14 install_dependencies.py --kind python

# Repair incomplete generated environments before starting a new experiment.
python3.14 install_dependencies.py --recreate

# Use an alternative Python package index, if needed.
python3.14 install_dependencies.py --index-url https://your-index.example/simple
```

Use the installer instead of installing the requirement files alone: it also
applies the required PM4Py changes and creates the runtime configuration.
Keep the two environments separate and do not request PM4Py extras, which can
impose different dependency versions. Allow around **1.5 GiB of free space** for
installation and temporary downloads, plus space for the input logs and results.

## Benchmark protocol and output

The five approaches are **Ebi Native Dynamic**, **Ebi Default**, **C++ Dynamic**,
**Python Dynamic**, and **corrected PM4Py process-tree alignment**. Native Dynamic
uses the preserved `native_dynamic_best_of_three_run_03` implementation; the other
four retain the original `best_of_three_run_01` interfaces. Special Alignments
is not executed.

Each year's 1,000 selected unique traces is aligned to its own four models
(`pt00`, `pt10`, `pt25`, `pt50`). Three complete passes are performed, choosing the
lowest validated time independently for each trace/model/algorithm combination.
Ties choose the earliest repetition. Jobs run sequentially; each algorithm/model/pass
uses a fresh process/model and three excluded warmups, followed by all selected
traces in the saved order. Each pass visits Native, C++, Default, Python Dynamic,
then PM4Py; each algorithm visits 2012 then 2017 with ascending model thresholds.
Rayon pools use one thread, `PYTHONHASHSEED=42`, and normal Python garbage collection.

Rust timers use `std::time::Instant`; Python-facing calls use `perf_counter_ns`.
File parsing, input construction, benchmark logging, validation and CSV output
are excluded. Native includes per-call tree preparation and standard alignment
assembly; Default excludes PTML-to-Petri-net conversion. C++/Python Dynamic return
costs while Ebi and PM4Py reconstruct alignments. Default retains its original
weighted A* objective and reports visible-deviation cost for comparison.

New results appear under:

```text
runs/reproduction_01/all_programs_trace_timings.csv
```

The run folder also contains every individual timing, raw jobs/logs, selected
minima, per-model and overall summaries, and independent verification. Smoke
results use `runs/smoke_01`. Smoke times are functional checks, not full-workload
performance results.

Additional run options:

```bash
python3.14 run_benchmark.py --data-dir /path/to/DATAFILES --prepare-only
python3.14 run_benchmark.py --data-dir /path/to/DATAFILES --resume
python3.14 run_benchmark.py --data-dir /path/to/DATAFILES --output reproduction_02
```

Use the same output name when resuming. Resume is supported between completed
jobs; partially written or failed jobs need inspection and are not silently
discarded. Keep the original package files unchanged while a prepared run is
active, because their hashes are checked. For a new machine, copy the distributed
files and run the installer again instead of copying `.venvs`. Saved run metadata
contains the paths of the machine where it was prepared.

The original full benchmark took roughly three hours on its original machine.
The complete historical reproduction package remains the place for archived
measurements and source rebuilds; this folder is the smaller execution package.

## Validation of this package

The download installer was tested on 21 September 2026. Both environments matched
their pinned inventories and passed dependency and interface checks. Reconstructing
the samples reproduced both original 1,000-trace XES files byte for byte. A smoke
benchmark completed all 360 measured calls and 360 warmups across the five algorithms,
eight models and three passes. Independent verification confirmed the saved costs
and best-of-three selection. The full benchmark was not rerun for this packaging change.

Validation evidence is saved separately in
`../download_dependency_validation_01/`. Generated environments and smoke outputs
were removed from this distribution after validation; its size is approximately
39 MiB, compared with the previous 341 MiB archive.
