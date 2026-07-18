from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ppath import (  # noqa: E402
    analyze_expression,
    analyze_pathways,
    load_study,
    map_expression,
    mapping_summary,
    render_report,
)


class PPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.study = load_study(ROOT / "studies/gse50760-demo.json")
        cls.mapped = map_expression(
            ROOT / "data/expression.tsv",
            ROOT / "data/gene-mappings.tsv",
            cls.study,
        )
        cls.stats = analyze_expression(cls.mapped, cls.study)
        cls.pathways = analyze_pathways(cls.stats, ROOT / "data/pathways.tsv", cls.study)

    def test_mapping_receipts_preserve_failures(self) -> None:
        summary = mapping_summary(self.mapped)
        self.assertEqual(summary["input_genes"], 12)
        self.assertEqual(summary["status_counts"]["ambiguous"], 1)
        self.assertEqual(summary["status_counts"]["unmapped"], 1)
        self.assertEqual(summary["duplicate_rows_collapsed"], 1)

    def test_duplicate_entrez_rows_are_collapsed(self) -> None:
        rows = {row["entrez_id"]: row for row in csv.DictReader(self.stats.splitlines(), delimiter="\t")}
        self.assertEqual(rows["7157"]["input_ids"], "TP53;TP53_ALT")
        self.assertEqual(rows["7157"]["collapsed_input_count"], "2")

    def test_pathways_are_deterministic_and_corrected(self) -> None:
        first = analyze_pathways(self.stats, ROOT / "data/pathways.tsv", self.study)
        second = analyze_pathways(self.stats, ROOT / "data/pathways.tsv", self.study)
        self.assertEqual(first, second)
        rows = list(csv.DictReader(first.splitlines(), delimiter="\t"))
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(0 <= float(row["fdr"]) <= 1 for row in rows))

    def test_report_exposes_scope_and_receipts(self) -> None:
        report = render_report(self.mapped, self.stats, self.pathways, self.study)
        self.assertIn("Mapping audit", report)
        self.assertIn("AMBIG1", report)
        self.assertIn("not for clinical decision-making", report)

    def test_study_json_is_stable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "study.json"
            target.write_text(json.dumps(self.study, sort_keys=True), encoding="utf-8")
            self.assertEqual(load_study(target)["organism"], "hsa")


if __name__ == "__main__":
    unittest.main()
