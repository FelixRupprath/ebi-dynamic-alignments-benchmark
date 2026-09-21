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
