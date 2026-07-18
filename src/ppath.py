#!/usr/bin/env python3
"""Deterministic gene mapping, pathway analysis, and provenance reports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

VERSION = "0.1.0"
UNUSABLE_MAPPING_STATUSES = {"ambiguous", "unmapped"}


class PPathError(ValueError):
    pass


def load_study(path: str | Path) -> dict:
    with Path(path).open(encoding="utf-8") as handle:
        study = json.load(handle)
    required = {
        "study_id",
        "title",
        "organism",
        "gene_id_column",
        "control_columns",
        "case_columns",
        "minimum_absolute_log2_fold_change",
    }
    missing = sorted(required - study.keys())
    if missing:
        raise PPathError(f"study is missing fields: {', '.join(missing)}")
    if not study["control_columns"] or not study["case_columns"]:
        raise PPathError("control_columns and case_columns must not be empty")
    return study


def read_tsv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise PPathError(f"{path} has no header")
        return list(reader.fieldnames), list(reader)


def write_tsv(rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fieldnames})
    return output.getvalue()


def parse_tsv(text: str) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def require_columns(columns: Iterable[str], required: Iterable[str], label: str) -> None:
    missing = sorted(set(required) - set(columns))
    if missing:
        raise PPathError(f"{label} is missing columns: {', '.join(missing)}")


def parse_number(value: str, field: str, gene_id: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise PPathError(f"{gene_id}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(number) or number < 0:
        raise PPathError(f"{gene_id}: {field} must be a finite non-negative number")
    return number


def map_expression(expression_path: str | Path, mapping_path: str | Path, study: Mapping) -> str:
    expression_columns, expression_rows = read_tsv(expression_path)
    mapping_columns, mapping_rows = read_tsv(mapping_path)
    gene_column = study["gene_id_column"]
    sample_columns = [*study["control_columns"], *study["case_columns"]]
    require_columns(expression_columns, [gene_column, *sample_columns], "expression table")
    require_columns(
        mapping_columns,
        ["organism", "input_id", "entrez_id", "symbol", "mapping_status", "source"],
        "mapping table",
    )

    mapping_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mapping_rows:
        if row["organism"] == study["organism"]:
            mapping_index[row["input_id"]].append(row)

    seen: set[str] = set()
    mapped_rows: list[dict[str, object]] = []
    for expression in expression_rows:
        input_id = expression[gene_column].strip()
        if not input_id:
            raise PPathError("expression table contains an empty gene identifier")
        if input_id in seen:
            raise PPathError(f"expression table contains duplicate identifier: {input_id}")
        seen.add(input_id)
        values = {column: parse_number(expression[column], column, input_id) for column in sample_columns}
        candidates = mapping_index.get(input_id, [])
        unique_candidates = {
            (candidate["entrez_id"], candidate["symbol"], candidate["mapping_status"], candidate["source"])
            for candidate in candidates
        }

        if not unique_candidates:
            receipt = {
                "entrez_id": "",
                "symbol": "",
                "mapping_status": "unmapped",
                "mapping_source": "",
                "candidate_entrez_ids": "",
                "used_in_analysis": "no",
            }
        elif len({candidate[0] for candidate in unique_candidates}) > 1:
            receipt = {
                "entrez_id": "",
                "symbol": "",
                "mapping_status": "ambiguous",
                "mapping_source": ";".join(sorted({candidate[3] for candidate in unique_candidates})),
                "candidate_entrez_ids": ";".join(sorted({candidate[0] for candidate in unique_candidates})),
                "used_in_analysis": "no",
            }
        else:
            entrez_id, symbol, status, source = sorted(unique_candidates)[0]
            receipt = {
                "entrez_id": entrez_id,
                "symbol": symbol,
                "mapping_status": status,
                "mapping_source": source,
                "candidate_entrez_ids": entrez_id,
                "used_in_analysis": "yes",
            }
        mapped_rows.append({"input_id": input_id, **receipt, **values})

    fields = [
        "input_id",
        "entrez_id",
        "symbol",
        "mapping_status",
        "mapping_source",
        "candidate_entrez_ids",
        "used_in_analysis",
        *sample_columns,
    ]
    return write_tsv(mapped_rows, fields)


def mapping_summary(mapped_text: str) -> dict:
    rows = parse_tsv(mapped_text)
    statuses = Counter(row["mapping_status"] for row in rows)
    usable = [row for row in rows if row["used_in_analysis"] == "yes"]
    unique_entrez = {row["entrez_id"] for row in usable}
    return {
        "input_genes": len(rows),
        "usable_mapping_rows": len(usable),
        "unique_entrez_ids": len(unique_entrez),
        "duplicate_rows_collapsed": len(usable) - len(unique_entrez),
        "status_counts": dict(sorted(statuses.items())),
    }


def analyze_expression(mapped_text: str, study: Mapping) -> str:
    rows = parse_tsv(mapped_text)
    control_columns = study["control_columns"]
    case_columns = study["case_columns"]
    pseudocount = float(study.get("pseudocount", 0.0))
    threshold = float(study["minimum_absolute_log2_fold_change"])
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["used_in_analysis"] == "yes":
            grouped[row["entrez_id"]].append(row)

    results: list[dict[str, object]] = []
    for entrez_id in sorted(grouped, key=lambda value: (int(value) if value.isdigit() else math.inf, value)):
        gene_rows = grouped[entrez_id]
        control_values = [float(row[column]) for row in gene_rows for column in control_columns]
        case_values = [float(row[column]) for row in gene_rows for column in case_columns]
        control_mean = sum(control_values) / len(control_values)
        case_mean = sum(case_values) / len(case_values)
        if control_mean + pseudocount <= 0 or case_mean + pseudocount <= 0:
            raise PPathError("pseudocount must make fold-change inputs positive")
        log2_fold_change = math.log2((case_mean + pseudocount) / (control_mean + pseudocount))
        results.append(
            {
                "entrez_id": entrez_id,
                "symbol": sorted({row["symbol"] for row in gene_rows})[0],
                "input_ids": ";".join(sorted(row["input_id"] for row in gene_rows)),
                "collapsed_input_count": len(gene_rows),
                "control_mean": format_float(control_mean),
                "case_mean": format_float(case_mean),
                "log2_fold_change": format_float(log2_fold_change),
                "direction": "up" if log2_fold_change > 0 else "down" if log2_fold_change < 0 else "unchanged",
                "changed": "yes" if abs(log2_fold_change) >= threshold else "no",
            }
        )
    fields = [
        "entrez_id",
        "symbol",
        "input_ids",
        "collapsed_input_count",
        "control_mean",
        "case_mean",
        "log2_fold_change",
        "direction",
        "changed",
    ]
    return write_tsv(results, fields)


def format_float(value: float) -> str:
    return f"{value:.8g}"


def hypergeometric_tail(population: int, successes: int, draws: int, observed: int) -> float:
    if population == 0 or draws == 0 or successes == 0:
        return 1.0
    denominator = math.comb(population, draws)
    maximum = min(successes, draws)
    minimum = max(observed, 0)
    probability = 0.0
    for overlap in range(minimum, maximum + 1):
        failures_drawn = draws - overlap
        if failures_drawn <= population - successes:
            probability += math.comb(successes, overlap) * math.comb(population - successes, failures_drawn) / denominator
    return min(probability, 1.0)


def adjust_bh(rows: list[dict[str, object]]) -> None:
    ordered = sorted(enumerate(rows), key=lambda item: (float(item[1]["p_value"]), item[1]["pathway_id"]))
    adjusted = [1.0] * len(rows)
    running = 1.0
    total = len(rows)
    for reverse_rank in range(total - 1, -1, -1):
        original_index, row = ordered[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, float(row["p_value"]) * total / rank)
        adjusted[original_index] = min(running, 1.0)
    for row, value in zip(rows, adjusted):
        row["fdr"] = format_float(value)


def analyze_pathways(stats_text: str, pathway_path: str | Path, study: Mapping) -> str:
    stats = parse_tsv(stats_text)
    pathway_columns, pathway_rows = read_tsv(pathway_path)
    require_columns(pathway_columns, ["organism", "pathway_id", "pathway_name", "entrez_ids"], "pathway table")
    genes = {row["entrez_id"]: row for row in stats}
    universe = set(genes)
    changed = {gene_id for gene_id, row in genes.items() if row["changed"] == "yes"}
    results: list[dict[str, object]] = []
    for pathway in pathway_rows:
        if pathway["organism"] != study["organism"]:
            continue
        pathway_genes = {value for value in pathway["entrez_ids"].split(";") if value}
        measured = pathway_genes & universe
        changed_in_pathway = measured & changed
        fold_changes = [float(genes[gene_id]["log2_fold_change"]) for gene_id in measured]
        p_value = hypergeometric_tail(len(universe), len(changed), len(measured), len(changed_in_pathway))
        results.append(
            {
                "pathway_id": pathway["pathway_id"],
                "pathway_name": pathway["pathway_name"],
                "pathway_gene_count": len(pathway_genes),
                "measured_gene_count": len(measured),
                "changed_gene_count": len(changed_in_pathway),
                "up_gene_count": sum(genes[gene_id]["direction"] == "up" for gene_id in changed_in_pathway),
                "down_gene_count": sum(genes[gene_id]["direction"] == "down" for gene_id in changed_in_pathway),
                "coverage": format_float(len(measured) / len(pathway_genes) if pathway_genes else 0.0),
                "mean_log2_fold_change": format_float(sum(fold_changes) / len(fold_changes) if fold_changes else 0.0),
                "changed_symbols": ";".join(sorted(genes[gene_id]["symbol"] for gene_id in changed_in_pathway)),
                "p_value": format_float(p_value),
                "fdr": "",
            }
        )
    adjust_bh(results)
    results.sort(key=lambda row: (float(row["fdr"]), float(row["p_value"]), row["pathway_id"]))
    fields = [
        "pathway_id",
        "pathway_name",
        "pathway_gene_count",
        "measured_gene_count",
        "changed_gene_count",
        "up_gene_count",
        "down_gene_count",
        "coverage",
        "mean_log2_fold_change",
        "changed_symbols",
        "p_value",
        "fdr",
    ]
    return write_tsv(results, fields)


def render_table(rows: Sequence[Mapping[str, str]], columns: Sequence[tuple[str, str]]) -> str:
    headings = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key, _ in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"<div class=table-wrap><table><thead><tr>{headings}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def render_report(mapped_text: str, stats_text: str, pathway_text: str, study: Mapping) -> str:
    mapped_rows = parse_tsv(mapped_text)
    stats_rows = parse_tsv(stats_text)
    pathway_rows = parse_tsv(pathway_text)
    summary = mapping_summary(mapped_text)
    significant = [row for row in pathway_rows if float(row["fdr"]) <= float(study.get("fdr_threshold", 0.05))]
    changed = [row for row in stats_rows if row["changed"] == "yes"]
    ambiguous = [row for row in mapped_rows if row["mapping_status"] in UNUSABLE_MAPPING_STATUSES]
    status_cards = "".join(
        f"<article><strong>{html.escape(str(value))}</strong><span>{html.escape(label)}</span></article>"
        for value, label in [
            (summary["input_genes"], "input identifiers"),
            (summary["unique_entrez_ids"], "unique Entrez IDs"),
            (len(changed), "changed genes"),
            (len(significant), "pathways below FDR threshold"),
        ]
    )
    status_list = "".join(
        f"<li><code>{html.escape(key)}</code><b>{value}</b></li>" for key, value in summary["status_counts"].items()
    )
    pathway_table = render_table(
        pathway_rows,
        [
            ("pathway_id", "Pathway"),
            ("pathway_name", "Name"),
            ("measured_gene_count", "Measured"),
            ("changed_gene_count", "Changed"),
            ("changed_symbols", "Changed genes"),
            ("p_value", "p"),
            ("fdr", "FDR"),
        ],
    )
    gene_table = render_table(
        stats_rows,
        [
            ("symbol", "Gene"),
            ("entrez_id", "Entrez"),
            ("input_ids", "Input IDs"),
            ("control_mean", study.get("control_label", "Control")),
            ("case_mean", study.get("case_label", "Case")),
            ("log2_fold_change", "log2 FC"),
            ("changed", "Changed"),
        ],
    )
    receipt_table = render_table(
        mapped_rows,
        [
            ("input_id", "Input"),
            ("entrez_id", "Entrez"),
            ("symbol", "Symbol"),
            ("mapping_status", "Status"),
            ("candidate_entrez_ids", "Candidates"),
            ("mapping_source", "Source"),
            ("used_in_analysis", "Used"),
        ],
    )
    warnings = "" if ambiguous else "<p>No ambiguous or unmapped identifiers were detected.</p>"
    return f"""<!doctype html>
<html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>{html.escape(study['title'])}</title>
<style>
:root{{--ink:#1c2924;--muted:#607069;--paper:#f4f0e6;--panel:#fffdf6;--accent:#176b57;--line:#d8d2c5;--warn:#9d4c2d}}*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 ui-monospace,SFMono-Regular,Menlo,monospace}}main{{max-width:1120px;margin:auto;padding:64px 24px 96px}}header{{border-left:8px solid var(--accent);padding:8px 0 8px 24px;margin-bottom:42px}}h1{{font:700 clamp(2rem,5vw,4rem)/1.02 Georgia,serif;margin:0 0 14px;max-width:900px}}h2{{font:700 1.7rem Georgia,serif;margin-top:52px}}p{{max-width:780px}}code{{background:#e9e4d8;padding:2px 5px}}.eyebrow{{color:var(--accent);font-weight:700;text-transform:uppercase;letter-spacing:.12em}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}}.cards article{{background:var(--panel);border:1px solid var(--line);padding:20px}}.cards strong{{display:block;font:700 2rem Georgia,serif;color:var(--accent)}}.cards span{{color:var(--muted);font-size:.85rem}}ul.status{{list-style:none;padding:0;display:flex;flex-wrap:wrap;gap:8px}}ul.status li{{background:var(--panel);border:1px solid var(--line);padding:9px 12px}}ul.status b{{margin-left:12px}}.table-wrap{{overflow:auto;border:1px solid var(--line);background:var(--panel)}}table{{border-collapse:collapse;width:100%;font-size:.82rem}}th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid var(--line);white-space:nowrap}}th{{position:sticky;top:0;background:#e9e4d8}}tr:last-child td{{border-bottom:0}}.note{{border:1px solid #d8b6a7;background:#fff5ef;padding:16px;color:#69341f}}footer{{margin-top:56px;padding-top:18px;border-top:1px solid var(--line);color:var(--muted);font-size:.8rem}}
</style></head><body><main>
<header><div class=eyebrow>ppath · mapping provenance report</div><h1>{html.escape(study['title'])}</h1><p><b>{html.escape(study.get('organism_name', study['organism']))}</b> · study <code>{html.escape(study['study_id'])}</code> · mapping snapshot <code>{html.escape(study.get('mapping_snapshot', 'unspecified'))}</code></p></header>
<section class=cards>{status_cards}</section>
<h2>Mapping audit</h2><p>Every input identifier receives a receipt. Ambiguous and unmapped identifiers are retained in this report but excluded from numerical analysis.</p><ul class=status>{status_list}</ul>{warnings}{receipt_table}
<h2>Gene-level comparison</h2><p>Input rows resolving to the same Entrez ID are collapsed before calculating the mean expression and log2 fold change. The changed threshold is <b>|log2 FC| ≥ {html.escape(str(study['minimum_absolute_log2_fold_change']))}</b>.</p>{gene_table}
<h2>Pathway over-representation</h2><p>One-sided hypergeometric tests compare changed genes with the measured-gene universe. Benjamini–Hochberg correction is applied across the displayed pathways.</p>{pathway_table}
<p class=note><b>Scope:</b> {html.escape(study.get('disclaimer', 'Research use only.'))}</p>
<footer>Generated deterministically by ppath {VERSION}. No network resource was consulted during this run.</footer>
</main></body></html>
"""


def sha256_path(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_manifest(study: Mapping, labelled_paths: Sequence[str]) -> str:
    artifacts = []
    for item in labelled_paths:
        if "=" not in item:
            raise PPathError(f"manifest file must be LABEL=PATH: {item}")
        label, path = item.split("=", 1)
        artifact_path = Path(path)
        artifacts.append(
            {
                "label": label,
                "sha256": sha256_path(artifact_path),
                "bytes": artifact_path.stat().st_size,
            }
        )
    manifest = {
        "schema": "ppath-provenance-v1",
        "ppath_version": VERSION,
        "study_id": study["study_id"],
        "organism": study["organism"],
        "mapping_snapshot": study.get("mapping_snapshot", "unspecified"),
        "artifacts": sorted(artifacts, key=lambda artifact: artifact["label"]),
    }
    return json.dumps(manifest, indent=2, sort_keys=True) + "\n"


def emit(value: str) -> None:
    sys.stdout.write(value)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ppath")
    subparsers = parser.add_subparsers(dest="command", required=True)

    map_parser = subparsers.add_parser("map", help="map expression identifiers and emit receipts")
    map_parser.add_argument("--study", required=True)
    map_parser.add_argument("--expression", required=True)
    map_parser.add_argument("--mappings", required=True)

    summary_parser = subparsers.add_parser("mapping-summary", help="summarize mapped identifiers")
    summary_parser.add_argument("--study", required=True)
    summary_parser.add_argument("--mapped", required=True)

    analyze_parser = subparsers.add_parser("analyze", help="calculate collapsed gene statistics")
    analyze_parser.add_argument("--study", required=True)
    analyze_parser.add_argument("--mapped", required=True)

    pathways_parser = subparsers.add_parser("pathways", help="test pathway over-representation")
    pathways_parser.add_argument("--study", required=True)
    pathways_parser.add_argument("--stats", required=True)
    pathways_parser.add_argument("--pathways", required=True)

    report_parser = subparsers.add_parser("report", help="render a deterministic HTML audit")
    report_parser.add_argument("--study", required=True)
    report_parser.add_argument("--mapped", required=True)
    report_parser.add_argument("--stats", required=True)
    report_parser.add_argument("--pathway-results", required=True)

    manifest_parser = subparsers.add_parser("manifest", help="hash labelled inputs and outputs")
    manifest_parser.add_argument("--study", required=True)
    manifest_parser.add_argument("--file", action="append", default=[], required=True)

    args = parser.parse_args(argv)
    study = load_study(args.study)
    if args.command == "map":
        emit(map_expression(args.expression, args.mappings, study))
    elif args.command == "mapping-summary":
        mapped = Path(args.mapped).read_text(encoding="utf-8")
        emit(json.dumps(mapping_summary(mapped), indent=2, sort_keys=True) + "\n")
    elif args.command == "analyze":
        emit(analyze_expression(Path(args.mapped).read_text(encoding="utf-8"), study))
    elif args.command == "pathways":
        emit(analyze_pathways(Path(args.stats).read_text(encoding="utf-8"), args.pathways, study))
    elif args.command == "report":
        emit(
            render_report(
                Path(args.mapped).read_text(encoding="utf-8"),
                Path(args.stats).read_text(encoding="utf-8"),
                Path(args.pathway_results).read_text(encoding="utf-8"),
                study,
            )
        )
    elif args.command == "manifest":
        emit(build_manifest(study, args.file))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PPathError as error:
        print(f"ppath: {error}", file=sys.stderr)
        raise SystemExit(2) from error
