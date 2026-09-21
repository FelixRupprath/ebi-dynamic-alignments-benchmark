#!/usr/bin/env python3
"""Independent, portable audit of the archived or reproduced five-method benchmark.

No benchmark/controller helpers are imported. Raw jobs are independently joined
to aggregate tables, minima, original sample identities and reference costs.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import statistics
import xml.etree.ElementTree as ET


BUNDLE = Path(__file__).resolve().parents[1]
ALGORITHMS = ("ebi_native_dynamic", "cpp", "ebi_default", "python_dynamic", "pm4py")
IDS = ("sample_index", "source_trace_index", "case_id", "variant_frequency", "trace_length", "trace_sha256")
DETAILS = ("elapsed_ns", "elapsed_ms", "alignment_cost", "reference_alignment_cost",
           "cost_matches_reference", "alignment_checks_passed", "status", "attempt_elapsed_ns", "error")
KEYS = ("dataset", "model_dataset", "model", "algorithm", "repetition")
OUTPUTS = {"all_trace_repetitions.csv", "minimum_trace_timings.csv", "all_programs_trace_timings.csv",
           "all_datasets_summary.csv", "overall_summary.csv", "failures.csv"}
NATIVE_FIELDS = IDS + ("elapsed_ns", "elapsed_ms", "alignment_cost", "reference_alignment_cost",
                      "alignment_moves", "replay_valid", "cost_matches_reference", "status", "error")
NATIVE_WARMUPS = "Samples 1, 2, 3, performed and replay-validated by preserved native runner"


def ensure(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def read(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        ensure(reader.fieldnames is not None and len(set(reader.fieldnames)) == len(reader.fieldnames),
               f"Missing/duplicate columns: {path}")
        ensure(all(None not in r and all(v is not None for v in r.values()) for r in rows), f"Malformed CSV: {path}")
        return rows


def key(row):
    return tuple(str(row[field]) for field in ("dataset", "model_dataset", "model", "sample_index"))


def float_equal(actual, expected, context):
    actual = float(actual)
    ensure(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-9),
           f"Incorrect {context}: got {actual}, expected {expected}")


def check_samples(root, metadata, count):
    manifests = {}
    for year, item in metadata["samples"].items():
        manifest = load(root / item["manifest"])
        ensure(len(manifest) == count, f"Manifest count differs: {year}")
        ensure(len({tuple(r["activities"]) for r in manifest}) == count, f"Repeated variant in {year}")
        log = ET.parse(root / item["sample_xes"]).getroot()
        traces = log.findall("{*}trace")
        ensure(len(traces) == count, f"XES count differs: {year}")
        for index, (record, trace) in enumerate(zip(manifest, traces), start=1):
            ensure(record["sample_index"] == index, "Nonsequential sample identity")
            ensure(record["trace_length"] == len(record["activities"]), "Manifest length differs")
            encoded = json.dumps(record["activities"], ensure_ascii=False, separators=(",", ":")).encode()
            ensure(hashlib.sha256(encoded).hexdigest() == record["trace_sha256"], "Manifest activity digest differs")
            case = trace.find("{*}string[@key='concept:name']")
            ensure(case is not None and case.get("value") == record["case_id"], "XES case identity differs")
            labels = []
            for event in trace.findall("{*}event"):
                activity = event.find("{*}string[@key='concept:name']")
                ensure(activity is not None, "XES event lacks concept:name")
                labels.append(activity.get("value"))
            ensure(labels == record["activities"], "XES activity sequence differs")
        manifests[year] = manifest
    return manifests


def row_valid(row, record, expected_cost, algorithm, context):
    ensure(set(row) == set(IDS + DETAILS), f"Unexpected raw schema: {context}")
    ensure(all(row[f] == str(record[f]) for f in IDS), f"Wrong identity: {context}")
    ensure(row["status"] == "ok" and row["cost_matches_reference"] == "1" and not row["error"],
           f"Unsuccessful row: {context}")
    ensure(row["alignment_cost"] == row["reference_alignment_cost"] == str(expected_cost),
           f"Reference cost mismatch: {context}")
    ensure(int(row["elapsed_ns"]) > 0 and int(row["attempt_elapsed_ns"]) == int(row["elapsed_ns"]),
           f"Invalid elapsed/attempt time: {context}")
    ensure(Decimal(row["elapsed_ms"]) * 1_000_000 == int(row["elapsed_ns"]), f"Elapsed-unit mismatch: {context}")
    expected_flag = "1" if algorithm in ("ebi_native_dynamic", "ebi_default", "pm4py") else ""
    ensure(row["alignment_checks_passed"] == expected_flag, f"Incorrect validation flag: {context}")


def computed_stats(minima, all_times):
    return dict(traces=len(minima), repetitions_per_trace=3,
                minimum_ms=min(minima) / 1e6, mean_ms=statistics.mean(minima) / 1e6,
                median_ms=statistics.median(minima) / 1e6,
                p95_ms=statistics.quantiles(minima, n=100, method="inclusive")[94] / 1e6,
                maximum_ms=max(minima) / 1e6, sum_minimum_seconds=sum(minima) / 1e9,
                sum_all_repetitions_seconds=sum(all_times) / 1e9)


def verify(root, reference):
    metadata = load(root / "metadata.json")
    ensure(metadata["status"] == "completed", "Results must be complete before verification")
    if reference:
        ensure(metadata.get("reference_result_kind") == "combined_historical_campaigns",
               "Reference results must identify their combined historical provenance")
    ensure(set(metadata["algorithms"]) == set(ALGORITHMS) and metadata["repetitions"] == 3,
           "Expected five algorithms and three repetitions")
    count = metadata["samples_per_dataset"]
    ensure(count == (3 if metadata.get("smoke_test", False) else 1000),
           "Expected exactly three smoke samples or all 1,000 preserved samples")
    pairs = metadata["pairs"]
    ensure(len(pairs) == 8 and len({p["name"] for p in pairs}) == 8, "Expected eight distinct models")
    ensure(all(p["dataset"] == p["model_dataset"] for p in pairs), "Expected matched-year models")
    ensure(metadata["warmup_traces_per_process"] == 3 and metadata["benchmark_parallelism"] == 1,
           "Unexpected warmup/parallelism protocol")
    manifests = check_samples(root, metadata, count)
    raw = read(root / "all_trace_repetitions.csv")
    ensure(len(raw) == metadata["expected_timed_alignments"] == count * 8 * 5 * 3, "Total measurement coverage differs")
    indexed = {key(r) + (r["algorithm"], r["repetition"]): r for r in raw}
    ensure(len(indexed) == len(raw), "Duplicate aggregate measurement keys")
    grouped = defaultdict(list)
    actual_jobs, job_hashes = set(), {}
    explicit_warmups = inferred_warmups = alignment_validations = 0
    for pair in pairs:
        manifest = manifests[str(pair["dataset"])]
        references = read(root / f"inputs/reference/bpi{pair['dataset']}_{pair['model']}.csv")
        ensure(len(references) == count, "Reference cost coverage differs")
        for record, ref in zip(manifest, references):
            ensure(all(ref[field] == str(record[field]) for field in IDS), "Reference trace identity differs")
        for algorithm in ALGORITHMS:
            for repetition in (1, 2, 3):
                relative = f"repetition_{repetition:02d}/{algorithm}/{pair['name']}"
                actual_jobs.add(relative)
                directory = root / "raw" / relative
                job = load(directory / "metadata.json")
                ensure(job["status"] == "completed" and job["completed"] == count and job["failures"] == 0,
                       f"Unsuccessful raw job: {relative}")
                expected_identity = dict(algorithm=algorithm, repetition=repetition, **pair)
                ensure(all(job[field] == value for field, value in expected_identity.items()), "Raw job identity differs")
                ensure(job["csv_sha256"] == sha(directory / "timings.csv"), f"Raw CSV digest differs: {relative}")
                job_hashes[relative] = {"csv_sha256": sha(directory / "timings.csv"),
                                        "metadata_sha256": sha(directory / "metadata.json")}
                rows = read(directory / "timings.csv")
                ensure(len(rows) == count, "Raw job trace coverage differs")
                for record, ref, row in zip(manifest, references, rows):
                    row_valid(row, record, int(ref["alignment_cost"]), algorithm, f"{relative}/{record['sample_index']}")
                    enriched = {"dataset": str(pair["dataset"]), "model_dataset": str(pair["model_dataset"]),
                                "model": pair["model"], "algorithm": algorithm, "repetition": str(repetition), **row}
                    identity = key(enriched) + (algorithm, str(repetition))
                    ensure(indexed.get(identity) == enriched, f"Aggregate differs from original raw row: {identity}")
                    grouped[key(enriched) + (algorithm,)].append(enriched)
                    alignment_validations += row["alignment_checks_passed"] == "1"
                float_equal(job["sum_alignment_seconds"], sum(int(r["elapsed_ns"]) for r in rows) / 1e9,
                            f"raw job time sum: {relative}")
                if algorithm == "ebi_native_dynamic":
                    ensure(job["warmups"] == NATIVE_WARMUPS,
                           "Native warmup/replay-validation record differs")
                    native = read(directory / "native_timings.csv")
                    ensure(len(native) == count, "Native Dynamic trace coverage differs")
                    for original, normalized in zip(native, rows):
                        ensure(set(original) == set(NATIVE_FIELDS), "Unexpected original native schema")
                        ensure(original["status"] == "ok" and original["replay_valid"] == "true"
                               and original["cost_matches_reference"] == "true" and not original["error"],
                               "Original native result did not pass cost/replay checks")
                        for field in IDS + ("elapsed_ns", "elapsed_ms", "alignment_cost", "reference_alignment_cost",
                                            "status", "error"):
                            ensure(original[field] == normalized[field], f"Native normalization differs: {relative}/{field}")
                        ensure(normalized["alignment_checks_passed"] == str(int(original["replay_valid"] == "true"))
                               and normalized["cost_matches_reference"] == str(int(original["cost_matches_reference"] == "true"))
                               and normalized["attempt_elapsed_ns"] == original["elapsed_ns"],
                               "Native boolean/attempt normalization differs")
                        ensure(int(original["alignment_moves"]) >= int(original["alignment_cost"]), "Invalid native move count")
                    if reference:
                        source = BUNDLE / "original_native_results" / "raw" / f"repetition_{repetition:02d}" / (
                            f"bpi{pair['dataset']}_{pair['model']}") / "timings.csv"
                        ensure(sha(directory / "native_timings.csv") == sha(source),
                               f"Native reference timings differ from the original run03 artifact: {relative}")
                    inferred_warmups += 3
                else:
                    ensure(len(job["warmups"]) == 3, "Explicit warmup count differs")
                    for index, warmup in enumerate(job["warmups"], start=1):
                        ensure(warmup["sample_index"] == index and warmup["status"] == "ok", "Warmup failed")
                        ensure(warmup["alignment_cost"] == int(references[index - 1]["alignment_cost"]), "Warmup reference cost differs")
                        ensure(int(warmup["elapsed_ns"]) > 0, "Warmup duration invalid")
                        if algorithm in ("ebi_default", "pm4py"):
                            ensure(warmup["alignment_checks_passed"] == 1, "Warmup alignment validation failed")
                        if algorithm == "ebi_default":
                            ensure(warmup["replay_valid"] is True, "Warmup replay failed")
                            ensure(warmup["weighted_alignment_cost"] == warmup["alignment_cost"] * 10000 + warmup["silent_moves"],
                                   "Warmup weighted cost differs")
                        explicit_warmups += 1
    ensure(len(actual_jobs) == 120 and set(metadata["completed_jobs"]) == actual_jobs
           and len(metadata["completed_jobs"]) == 120, "Completed-job coverage differs")
    actual_paths = {str(p.parent.relative_to(root / "raw")) for p in (root / "raw").glob("*/*/*/timings.csv")}
    ensure(actual_paths == actual_jobs, "Raw directory coverage differs")
    ensure(sum(len(g) for g in grouped.values()) == len(raw), "Raw rows omitted from independent join")

    minima = {}
    for group, candidates in grouped.items():
        ensure(sorted(r["repetition"] for r in candidates) == ["1", "2", "3"], "Missing/duplicate repetition")
        ensure(len({r["alignment_cost"] for r in candidates}) == 1, "Cost changes across repetitions")
        minima[group] = min(candidates, key=lambda r: (int(r["elapsed_ns"]), int(r["repetition"])))
    minimum = read(root / "minimum_trace_timings.csv")
    ensure(len(minimum) == len(minima) == count * 40, "Minimum table coverage differs")
    ensure(len({key(r) + (r["algorithm"],) for r in minimum}) == len(minimum), "Duplicate selected minima")
    for row in minimum:
        group = key(row) + (row["algorithm"],)
        selected = minima[group]
        expected = dict(selected, **{f"repetition_{trial['repetition']}_elapsed_ns": trial["elapsed_ns"]
                                     for trial in grouped[group]})
        ensure(row == expected, f"Selected minimum/earliest tie differs: {group}")
    wide = read(root / "all_programs_trace_timings.csv")
    ensure(len(wide) == count * 8 and len({key(r) for r in wide}) == len(wide), "Wide table coverage differs")
    for row in wide:
        record = manifests[row["dataset"]][int(row["sample_index"]) - 1]
        expected = {field: row[field] for field in ("dataset", "model_dataset", "model")}
        expected.update({field: str(record[field]) for field in IDS})
        costs = set()
        for algorithm in ALGORITHMS:
            group = key(row) + (algorithm,)
            best = minima[group]
            costs.add(best["alignment_cost"])
            expected.update({f"{algorithm}_elapsed_ns": best["elapsed_ns"],
                             f"{algorithm}_elapsed_ms": best["elapsed_ms"],
                             f"{algorithm}_selected_repetition": best["repetition"]})
            expected.update({f"{algorithm}_repetition_{trial['repetition']}_elapsed_ns": trial["elapsed_ns"]
                             for trial in grouped[group]})
        ensure(len(costs) == 1, "Five algorithms disagree on visible cost")
        expected["alignment_cost"] = next(iter(costs))
        ensure(row == expected, f"Wide table differs: {key(row)}")

    model_summaries = read(root / "all_datasets_summary.csv")
    model_keys = {(str(p["dataset"]), str(p["model_dataset"]), p["model"], a) for p in pairs for a in ALGORITHMS}
    observed_keys = [(r["dataset"], r["model_dataset"], r["model"], r["algorithm"]) for r in model_summaries]
    ensure(len(observed_keys) == 40 and set(observed_keys) == model_keys, "Model summary coverage differs")
    for row in model_summaries:
        prefix = (row["dataset"], row["model_dataset"], row["model"])
        selected = [int(value["elapsed_ns"]) for group, value in minima.items() if group[:3] == prefix and group[4] == row["algorithm"]]
        attempts = [int(trial["elapsed_ns"]) for group, trials in grouped.items()
                    if group[:3] == prefix and group[4] == row["algorithm"] for trial in trials]
        expected = computed_stats(selected, attempts)
        ensure(len(selected) == count, "Per-model selected coverage differs")
        ensure(set(row) == {"dataset", "model_dataset", "model", "algorithm"} | set(expected), "Model summary schema differs")
        for field, value in expected.items():
            float_equal(row[field], value, f"model summary/{prefix}/{row['algorithm']}/{field}")
    overall = read(root / "overall_summary.csv")
    ensure(len(overall) == 5 and {r["algorithm"] for r in overall} == set(ALGORITHMS), "Overall summary coverage differs")
    complete_overall = {}
    for row in overall:
        algorithm = row["algorithm"]
        selected = [int(v["elapsed_ns"]) for group, v in minima.items() if group[4] == algorithm]
        attempts = [int(trial["elapsed_ns"]) for group, trials in grouped.items() if group[4] == algorithm for trial in trials]
        expected = computed_stats(selected, attempts)
        complete_overall[algorithm] = expected
        ensure(set(row) == {"algorithm"} | (set(expected) - {"minimum_ms", "sum_all_repetitions_seconds"}),
               "Historical overall summary schema differs")
        for field in set(row) - {"algorithm"}:
            float_equal(row[field], expected[field], f"overall summary/{algorithm}/{field}")
    ensure(not read(root / "failures.csv"), "Failure rows are present")
    verification = load(root / "verification.json")
    ensure(verification == metadata["verification"], "Embedded/separate verification differs")
    expected_verification = dict(timed_alignments=len(raw), minimum_rows=len(minimum), comparison_rows=len(wide),
                                 dataset_algorithm_summary_rows=40, failures=0)
    ensure(all(verification[k] == value for k, value in expected_verification.items()), "Verification counts differ")
    for flag in ("all_three_repetitions_successful", "all_costs_match", "all_selected_times_are_minima",
                 "identifiers_match_original_samples", "protected_inputs_sources_and_previous_final_unchanged"):
        ensure(verification[flag] is True, f"Verification flag false: {flag}")
    ensure(set(verification["outputs_sha256"]) == OUTPUTS, "Output hash coverage differs")
    output_hashes = {name: sha(root / name) for name in sorted(OUTPUTS)}
    ensure(output_hashes == verification["outputs_sha256"], "Output artifact hashes differ")
    protected_count = 0
    old_bundle = Path(metadata.get("reproduction_package", str(BUNDLE)))
    for path, expected in metadata["protected_files_sha256"].items():
        # Combined references use package-relative paths. New measurements may
        # record absolute paths under the package's then-current root; rebase
        # those paths after relocation without inspecting external originals.
        saved = Path(path)
        if saved.is_absolute():
            ensure(saved.is_relative_to(old_bundle), f"Protected path lies outside reproduction package: {path}")
            current = BUNDLE / saved.relative_to(old_bundle)
        else:
            current = BUNDLE / saved
        current = current.resolve()
        ensure(current.is_relative_to(BUNDLE), f"Protected path escapes reproduction package: {path}")
        ensure(sha(current) == expected, f"Protected reproduction input/source changed: {current}")
        protected_count += 1
    return dict(status="passed", verified_utc=datetime.now(timezone.utc).isoformat(),
                reference_results=reference, result_directory=str(root.relative_to(BUNDLE)),
                reference_result_kind=metadata.get("reference_result_kind"),
                smoke_test=metadata.get("smoke_test", False), samples_per_dataset=count,
                raw_jobs=len(actual_jobs), raw_repetitions=len(raw), selected_minima=len(minimum), comparison_rows=len(wide),
                explicit_successful_warmup_records=explicit_warmups,
                native_dynamic_warmups_inferred_from_completed_runner=inferred_warmups,
                total_warmup_calls=explicit_warmups + inferred_warmups,
                timed_alignment_check_records=alignment_validations,
                exactly_three_successful_repetitions_per_combination=True,
                independently_recomputed_minima_and_earliest_ties_match=True,
                raw_job_csvs_match_aggregate=True, wide_csv_matches_raw_repetitions=True,
                all_identifiers_match_manifests_and_xes=True, all_five_algorithms_agree_on_visible_cost=True,
                all_costs_match_saved_references=True, model_summary_rows=40, overall_summary_rows=5,
                fully_recomputed_overall_statistics=complete_overall,
                protected_reproduction_files_checked=protected_count,
                protection_scope="All recorded package-relative or reproduction-root absolute paths checked, rebased within the package if relocated; external historical paths are never accessed.",
                warmup_note="Native Dynamic validates the cost and native replay of its three warmups before measuring traces; successful runner completion establishes those checks, but individual warmup durations are not recorded. Other methods retain explicit warmup results.",
                raw_artifact_hashes=job_hashes, outputs_sha256=output_hashes,
                metadata_sha256=sha(root / "metadata.json"), validator_sha256=sha(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--reference", action="store_true", help="Audit preserved reference results (default)")
    group.add_argument("--results", help="Completed result directory within this package, e.g. runs/smoke_01")
    parser.add_argument("--output", help="Audit JSON path within the package; historical results cannot be overwritten")
    args = parser.parse_args()
    root = Path(args.results) if args.results else BUNDLE / "reference_results"
    if not root.is_absolute():
        root = BUNDLE / root
    root = root.resolve()
    ensure(root.is_relative_to(BUNDLE), "Result directory must be inside the package")
    reference = root == BUNDLE / "reference_results"
    output = Path(args.output) if args.output else (BUNDLE / "validation/reference_verification.json" if reference
                                                   else root / "reproduction_verification.json")
    if not output.is_absolute():
        output = BUNDLE / output
    output = output.resolve()
    ensure(output.is_relative_to(BUNDLE) and not output.is_relative_to(BUNDLE / "reference_results"),
           "Audit output must stay inside package and outside preserved reference results")
    ensure(output.suffix == ".json" and output.name not in {"metadata.json", "verification.json", "bundle_manifest.json"},
           "Choose a separate JSON audit output")
    result = verify(root, reference)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("raw_artifact_hashes", "outputs_sha256", "fully_recomputed_overall_statistics")}, indent=2))
    print(f"Audit saved to {output}")


if __name__ == "__main__":
    main()
