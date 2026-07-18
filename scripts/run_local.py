#!/usr/bin/env python3
"""Execute the reference pipeline without pp for development and testing."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ppath import (  # noqa: E402
    analyze_expression,
    analyze_pathways,
    build_manifest,
    load_study,
    map_expression,
    mapping_summary,
    render_report,
)


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=ROOT / "build")
    parser.add_argument("--study", type=Path, default=ROOT / "studies/gse50760-demo.json")
    args = parser.parse_args()
    expression = ROOT / "data/expression.tsv"
    mappings = ROOT / "data/gene-mappings.tsv"
    pathways = ROOT / "data/pathways.tsv"
    study = load_study(args.study)

    mapped = map_expression(expression, mappings, study)
    stats = analyze_expression(mapped, study)
    pathway_results = analyze_pathways(stats, pathways, study)
    report = render_report(mapped, stats, pathway_results, study)
    outputs = {
        "mapping-receipts.tsv": mapped,
        "mapping-summary.json": json.dumps(mapping_summary(mapped), indent=2, sort_keys=True) + "\n",
        "gene-statistics.tsv": stats,
        "pathway-results.tsv": pathway_results,
        "report.html": report,
    }
    for name, content in outputs.items():
        write(args.out / name, content)
    manifest_items = [
        f"study={args.study}",
        f"expression={expression}",
        f"mappings={mappings}",
        f"pathways={pathways}",
        f"python-tool={ROOT / 'src/ppath.py'}",
        *(f"{name}={args.out / name}" for name in outputs),
    ]
    write(args.out / "provenance.json", build_manifest(study, manifest_items))
    print(f"built {len(outputs) + 1} artifacts in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
