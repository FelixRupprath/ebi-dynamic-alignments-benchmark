#!/usr/bin/env python3
"""Download pinned Python packages into the benchmark's isolated environments.

Reuses the preserved C++ extension and Python/PM4Py corrections. Environment
archives are not required. Installation downloads packages from the configured
Python package index; --verify-only performs no downloads or installations.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
ENVIRONMENT = ROOT / "environment"
ORIGINAL_CPP_SHA256 = "8633957fa162189e013c7e6134ed0b771dbc0b9863cccb7e3a6cc6da1de228c7"
DYNAMIC_SOURCE_SHA256 = "4676bedcd6997a5a7d2764b5cbbac3bccc1678e47b26a6f9df2a3c2cd341c9c5"
MARKER = "BENCHMARK_ENVIRONMENT_JSON="
PROBE = r'''
import importlib.metadata as m, json, platform, sys, sysconfig
print("BENCHMARK_ENVIRONMENT_JSON=" + json.dumps(dict(
    version=platform.python_version(), implementation=platform.python_implementation(),
    machine=platform.machine(), system=platform.system(), libc=platform.libc_ver(),
    soabi=sysconfig.get_config_var("SOABI"), extension_suffix=sysconfig.get_config_var("EXT_SUFFIX"),
    abiflags=sys.abiflags, prefix=sys.prefix, purelib=sysconfig.get_paths()["purelib"],
    distributions={d.metadata["Name"]: d.version for d in m.distributions()})))
'''


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def normalize(name):
    return name.lower().replace("_", "-").replace(".", "-")


def child_environment():
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="42",
                       MPLCONFIGDIR=str(ROOT / ".cache" / "matplotlib"),
                       PIP_DISABLE_PIP_VERSION_CHECK="1")
    environment.pop("PYTHONPATH", None)
    environment.pop("PYTHONHOME", None)
    return environment


def run_json(python, code, *arguments):
    result = subprocess.run([str(python), "-B", "-c", code, *map(str, arguments)],
                            cwd=ROOT, env=child_environment(), text=True, capture_output=True)
    require(result.returncode == 0,
            f"Python environment check failed ({python}):\n{result.stdout}\n{result.stderr}")
    encoded = [line[len(MARKER):] for line in result.stdout.splitlines() if line.startswith(MARKER)]
    require(len(encoded) == 1, f"Environment probe did not return one result: {result.stdout}")
    return json.loads(encoded[0])


def validate_platform(actual, expected):
    fields = ("version", "implementation", "machine", "system", "soabi", "abiflags")
    differences = {field: {"required": expected[field], "actual": actual[field]}
                   for field in fields if actual[field] != expected[field]}
    require(not differences,
            "This benchmark requires the original CPython 3.14.7 Linux x86_64 ABI. "
            "Select that interpreter with --python. Differences: " + json.dumps(differences))


def validate_distributions(actual, expected):
    present = {normalize(name): version for name, version in actual.items()}
    required = {normalize(name): version for name, version in expected.items()}
    differences = dict(missing={name: version for name, version in required.items() if name not in present},
                       changed={name: {"required": version, "actual": present[name]}
                                for name, version in required.items()
                                if name in present and present[name] != version},
                       extras={name: version for name, version in present.items() if name not in required})
    require(not any(differences.values()), "Installed distribution inventory differs: " + json.dumps(differences))
    return differences


def read_requirements(expected):
    path = ENVIRONMENT / expected["requirements_file"]
    require(digest(path) == expected["requirements_sha256"], f"Requirements file changed: {path}")
    pins = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        require(re.fullmatch(r"[A-Za-z0-9._-]+==[A-Za-z0-9.!+_-]+", line),
                f"Expected a single exact package pin, found {line!r} in {path}")
        name, version = line.split("==", 1)
        key = normalize(name)
        require(key not in pins, f"Duplicate package pin in {path}: {name}")
        pins[key] = (name, version)
    require({key: version for key, (_, version) in pins.items()} ==
            {normalize(name): version for name, version in expected["distributions"].items()},
            f"Requirements do not match the historical distribution inventory: {path}")
    return pins


def pip_command(python):
    return [str(python), "-B", "-m", "pip", "--isolated", "--disable-pip-version-check", "--no-input"]


def install_pins(kind, python, pins, args):
    """Install the complete inventory without resolving or adding dependencies.

    Every compiled dependency must have a compatible wheel. The historical
    C++ environment also contains the unrelated pure-Python package alignment;
    it may be built from its exact pinned source release in pip's temporary
    isolated build environment if the index does not provide a wheel. This is
    not the benchmark's preserved alignment extension.
    """
    common = ["install", "--no-cache-dir", "--no-compile", "--no-deps",
              "--index-url", args.index_url]
    # ensurepip only bootstraps the installer; install the recorded pip version
    # before installing the other exact pins. Build dependencies, if needed
    # for the pure-Python alignment package, stay outside the target venv.
    subprocess.run([str(python), "-B", "-m", "ensurepip"], cwd=ROOT,
                   env=child_environment(), check=True)
    subprocess.run([*pip_command(python), *common, "--only-binary=:all:",
                    "--report", str(ENVIRONMENT / f"install_{kind}_pip.json"),
                    "pip==" + pins["pip"][1]],
                   cwd=ROOT, env=child_environment(), check=True)
    wheel_pins = [f"{name}=={version}" for key, (name, version) in pins.items()
                  if key not in {"pip", "alignment"}]
    if wheel_pins:
        subprocess.run([*pip_command(python), *common, "--only-binary=:all:",
                        "--report", str(ENVIRONMENT / f"install_{kind}_wheels.json"), *wheel_pins],
                       cwd=ROOT, env=child_environment(), check=True)
    if "alignment" in pins:
        name, version = pins["alignment"]
        subprocess.run([*pip_command(python), *common, "--prefer-binary",
                        "--report", str(ENVIRONMENT / f"install_{kind}_alignment.json"),
                        f"{name}=={version}"], cwd=ROOT, env=child_environment(), check=True)


def check_environment(kind, expected):
    directory = ROOT / ".venvs" / kind
    python = directory / "bin" / "python"
    require(directory.resolve().is_relative_to(ROOT.resolve()) and not directory.is_symlink(),
            "Generated environment must remain inside this benchmark directory and not be a symlink")
    require(directory.is_dir() and python.is_file() and (directory / "pyvenv.cfg").is_file(),
            f"Incomplete environment: {directory}; run setup with --kind {kind} --recreate")
    try:
        info = run_json(python, PROBE)
        require(Path(info["prefix"]).resolve() == directory.resolve(),
                f"Interpreter does not belong to the expected environment: {directory}")
        require(Path(info["purelib"]).resolve().is_relative_to(directory.resolve()),
                "Environment site-packages escaped its directory")
        validate_platform(info, expected)
        info["distribution_differences"] = validate_distributions(info["distributions"], expected["distributions"])
        require(not list(Path(info["purelib"]).glob("*.pth")), "Unexpected .pth import dependency")
        result = subprocess.run([*pip_command(python), "check"], cwd=ROOT,
                                env=child_environment(), capture_output=True, text=True)
        require(result.returncode == 0, f"Package dependency check failed: {result.stdout}\n{result.stderr}")
        info["pip_check"] = result.stdout.strip()
    except RuntimeError as error:
        raise RuntimeError(f"{error}\nIf this environment is incomplete, run setup with --kind {kind} --recreate.") from error
    return python, info


def create_environment(kind, base_python, expected, pins, args):
    directory = ROOT / ".venvs" / kind
    require(directory.resolve().is_relative_to(ROOT.resolve()) and not directory.is_symlink(),
            "Generated environment must remain inside this benchmark directory and not be a symlink")
    if args.recreate and directory.exists():
        require(directory.is_dir(), f"Environment path is not a directory: {directory}")
        shutil.rmtree(directory)
    if not directory.exists():
        require(not args.verify_only, f"Environment missing: {directory}; run setup first")
        print(f"Downloading pinned packages into {directory.relative_to(ROOT)}", flush=True)
        subprocess.run([str(base_python), "-m", "venv", "--without-pip", str(directory)],
                       cwd=ROOT, env=child_environment(), check=True)
        python = directory / "bin" / "python"
        info = run_json(python, PROBE)
        require(Path(info["prefix"]).resolve() == directory.resolve(), "Unexpected venv interpreter prefix")
        require(Path(info["purelib"]).resolve().is_relative_to(directory.resolve()), "Unexpected venv site-packages location")
        install_pins(kind, python, pins, args)
    return check_environment(kind, expected)


def apply_overlay(python_info, manifest, verify_only):
    site = Path(python_info["purelib"])
    files = manifest["pm4py_overlay"]["files_sha256"]
    for relative, expected in files.items():
        source = ROOT / relative
        require(digest(source) == expected, f"Changed PM4Py overlay: {relative}")
        target = site / source.relative_to(ROOT / "sources")
        require(target.is_file(), f"Expected installed PM4Py module missing: {target}")
        if not verify_only:
            shutil.copy2(source, target)
            # Never use bytecode compiled before an overlay replacement.
            for cache in target.parent.glob(f"__pycache__/{target.stem}.*.pyc"):
                cache.unlink()
        require(digest(target) == expected, f"Installed PM4Py module differs from the preserved overlay: {target}")
    return files


CPP_SMOKE = r'''
import importlib.util, json, sys
spec=importlib.util.spec_from_file_location("alignment", sys.argv[1])
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
solver=module.AlignmentWrapper()
solver.loadTree("->('a','b')")
costs=[solver.align(trace) for trace in (["a","b"], ["a"], ["x"])]
assert costs == [0,1,3], costs
print("BENCHMARK_ENVIRONMENT_JSON="+json.dumps(dict(module_file=module.__file__,costs=costs)))
'''


PYTHON_SMOKE = r'''
import hashlib, importlib.util, json, sys
import pm4py
from pm4py.objects.process_tree.utils.generic import parse
from pm4py.objects.log.obj import Trace, Event
from pm4py.algo.conformance.alignments.process_tree import algorithm
spec=importlib.util.spec_from_file_location("preserved_dynamic_alignment",sys.argv[1])
dynamic=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=dynamic
spec.loader.exec_module(dynamic)
assert pm4py.__version__ == "2.7.23.8"
assert algorithm.DEFAULT_VARIANT.name == "SEARCH_GRAPH_PT"
tree=parse("->('a','b')")
dynamic_costs=[dynamic.dyn_align(tuple(trace),tree) for trace in (["a","b"],["a"],["x"])]
assert dynamic_costs == [0,1,3],dynamic_costs
pm4py_costs=[]
for activities,expected in [(["a","b"],0),(["a"],1),(["x"],3)]:
    trace=Trace([Event({"concept:name":a}) for a in activities])
    result=algorithm.apply(trace,tree)
    assert result["cost"] == expected,result
    assert result["optimal"] is True,result
    pm4py_costs.append(result["cost"])
print("BENCHMARK_ENVIRONMENT_JSON="+json.dumps(dict(pm4py_version=pm4py.__version__,
    default_variant=algorithm.DEFAULT_VARIANT.name,dynamic_costs=dynamic_costs,
    pm4py_costs=pm4py_costs,search_graph_pt_file=algorithm.DEFAULT_VARIANT.value.__file__)))
'''


def save_runtime(data):
    path = ENVIRONMENT / "runtime.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable, help="CPython 3.14.7 base interpreter")
    parser.add_argument("--index-url", default="https://pypi.org/simple",
                        help="Python package index used only when creating an environment")
    parser.add_argument("--kind", choices=("all", "cpp", "python"), default="all",
                        help="Create or recreate one environment at a time, or both (default)")
    parser.add_argument("--verify-only", action="store_true",
                        help="Check installed environments without downloading, installing or replacing packages")
    parser.add_argument("--recreate", action="store_true",
                        help="Replace only the selected generated .venvs/cpp and/or .venvs/python")
    args = parser.parse_args()
    require(not (args.verify_only and args.recreate), "--verify-only cannot recreate environments")
    require(ENVIRONMENT.is_dir(), f"Environment configuration directory is missing: {ENVIRONMENT}")
    save_runtime(dict(status="checking", setup_mode="pip", checked_utc=datetime.now(timezone.utc).isoformat()))
    try:
        manifest = json.loads((ENVIRONMENT / "manifest.json").read_text())
        pins = {kind: read_requirements(manifest["environments"][kind]) for kind in ("cpp", "python")}
        base_python = Path(shutil.which(args.python) or args.python).absolute()
        base = run_json(base_python, PROBE)
        validate_platform(base, manifest["environments"]["python"])
        selected = ("cpp", "python") if args.kind == "all" else (args.kind,)
        interpreters, environments = {}, {}
        for kind in selected:
            interpreters[kind], environments[kind] = create_environment(
                kind, base_python, manifest["environments"][kind], pins[kind], args)
        # A single-environment invocation may establish global readiness only
        # after checking the other installed environment too; never trust a
        # previous runtime.json or mark a half-created installation ready.
        for kind in ("cpp", "python"):
            if kind not in environments and (ROOT / ".venvs" / kind).exists():
                interpreters[kind], environments[kind] = check_environment(kind, manifest["environments"][kind])
        overlay = {}
        smoke_checks = {}
        module = ROOT / "bin" / ("alignment" + manifest["environments"]["cpp"]["extension_suffix"])
        require(module.is_file(), f"Preserved C++ extension missing: {module}")
        require(digest(module) == ORIGINAL_CPP_SHA256, "Preserved C++ extension hash mismatch")
        dynamic = ROOT / "sources" / "dynamic_alignment.py"
        require(digest(dynamic) == DYNAMIC_SOURCE_SHA256, "Preserved Python Dynamic source hash mismatch")
        if "python" in environments:
            overlay = apply_overlay(environments["python"], manifest,
                                    args.verify_only or "python" not in selected)
            smoke_checks["python_and_pm4py"] = run_json(interpreters["python"], PYTHON_SMOKE, dynamic)
        if "cpp" in environments:
            cpp_variant = Path(environments["cpp"]["purelib"]) / "pm4py/algo/conformance/alignments/process_tree/variants/search_graph_pt.py"
            require(digest(cpp_variant) == manifest["environments"]["cpp"]["search_graph_pt_sha256"],
                    "The C++ environment's original, unpatched PM4Py module has changed")
            try:
                smoke_checks["cpp"] = run_json(interpreters["cpp"], CPP_SMOKE, module)
            except RuntimeError as error:
                raise RuntimeError(f"{error}\nThe preserved extension requires compatible Linux x86_64 system libraries "
                                   "and was built with glibc 2.44. This package does not rebuild binaries.") from error
        ready = set(environments) == {"cpp", "python"}
        executable_paths = {}
        if "cpp" in interpreters:
            executable_paths["cpp"] = str(interpreters["cpp"].relative_to(ROOT))
        if "python" in interpreters:
            executable_paths["python_dynamic"] = str(interpreters["python"].relative_to(ROOT))
            executable_paths["pm4py"] = executable_paths["python_dynamic"]
        runtime = dict(status="ready" if ready else "partial", validated_utc=datetime.now(timezone.utc).isoformat(),
                       setup_mode="pip", python_executables=executable_paths,
                       cpp_module=str(module.relative_to(ROOT)), cpp_module_sha256=digest(module),
                       cpp_original_binary_reused=True,
                       environment_manifest_sha256=digest(ENVIRONMENT / "manifest.json"),
                       environments=environments, pm4py_overlay_sha256=overlay,
                       smoke_checks=smoke_checks,
                       requirements_sha256={kind: digest(ENVIRONMENT / manifest["environments"][kind]["requirements_file"])
                                            for kind in ("cpp", "python")},
                       pip_install_reports_sha256={str(path.relative_to(ROOT)): digest(path)
                                                  for path in sorted(ENVIRONMENT.glob("install_*.json"))},
                       package_provenance_note="Exact historical versions; downloaded wheels can differ from the historical installed bytes. "
                                               "Pip reports record downloaded artifact URLs and hashes. The preserved alignment binaries and PM4Py corrections are hash checked.",
                       relocation_note="Run setup or --verify-only after moving the benchmark directory.")
        save_runtime(runtime)
    except Exception as error:
        save_runtime(dict(status="failed", setup_mode="pip", checked_utc=datetime.now(timezone.utc).isoformat(),
                          error=str(error)))
        raise
    if ready:
        print("Environment ready: exact package inventories, dependencies, PM4Py corrections and three alignment interfaces verified.")
    else:
        missing = sorted({"cpp", "python"} - set(environments))
        print("Selected environment verified. Full benchmark setup is incomplete; run setup with --kind " + missing[0] + ".")
    print(ENVIRONMENT / "runtime.json")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError, OSError) as error:
        raise SystemExit(str(error)) from error
