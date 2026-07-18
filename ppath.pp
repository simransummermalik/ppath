# ppath.pp — the content-addressed orchestration layer for one pathway study.

def output-or-error(result) {
  if hash-map-get(result, "exit") = 0 {
    hash-map-get(result, "out")
  } else {
    error(hash-map-get(result, "err"))
  }
}

def map-identifiers(python, tool, study, expression, mappings) {
  node {
    output-or-error(perform run(
      python, tool, "map",
      "--study", study,
      "--expression", expression,
      "--mappings", mappings))
  }
}

def summarize-mappings(python, tool, study, mapped) {
  node {
    perform write-file("mapped.tsv", mapped)
    output-or-error(perform run(
      python, tool, "mapping-summary",
      "--study", study,
      "--mapped", "mapped.tsv"))
  }
}

def analyze-genes(python, tool, study, mapped) {
  node {
    perform write-file("mapped.tsv", mapped)
    output-or-error(perform run(
      python, tool, "analyze",
      "--study", study,
      "--mapped", "mapped.tsv"))
  }
}

def analyze-study-pathways(python, tool, study, pathways, statistics) {
  node {
    perform write-file("statistics.tsv", statistics)
    output-or-error(perform run(
      python, tool, "pathways",
      "--study", study,
      "--stats", "statistics.tsv",
      "--pathways", pathways))
  }
}

def render-study-report(python, tool, study, mapped, statistics, pathway-results) {
  node {
    perform write-file("mapped.tsv", mapped)
    perform write-file("statistics.tsv", statistics)
    perform write-file("pathway-results.tsv", pathway-results)
    output-or-error(perform run(
      python, tool, "report",
      "--study", study,
      "--mapped", "mapped.tsv",
      "--stats", "statistics.tsv",
      "--pathway-results", "pathway-results.tsv"))
  }
}

def create-manifest(python, tool, study, expression, mappings, pathways, mapped, summary, statistics, pathway-results, report) {
  node {
    perform write-file("mapping-receipts.tsv", mapped)
    perform write-file("mapping-summary.json", summary)
    perform write-file("gene-statistics.tsv", statistics)
    perform write-file("pathway-results.tsv", pathway-results)
    perform write-file("report.html", report)
    output-or-error(perform run(
      python, tool, "manifest",
      "--study", study,
      "--file", string-append("study=", study),
      "--file", string-append("expression=", expression),
      "--file", string-append("mappings=", mappings),
      "--file", string-append("pathways=", pathways),
      "--file", string-append("python-tool=", tool),
      "--file", "mapping-receipts.tsv=mapping-receipts.tsv",
      "--file", "mapping-summary.json=mapping-summary.json",
      "--file", "gene-statistics.tsv=gene-statistics.tsv",
      "--file", "pathway-results.tsv=pathway-results.tsv",
      "--file", "report.html=report.html"))
  }
}

let (root = car(argv()),
     python = car(cdr(argv())),
     tool = string-append(root, "/src/ppath.py"),
     study = string-append(root, "/studies/gse50760-demo.json"),
     expression = string-append(root, "/data/expression.tsv"),
     mappings = string-append(root, "/data/gene-mappings.tsv"),
     pathways = string-append(root, "/data/pathways.tsv"),
     mapped = map-identifiers(python, tool, study, expression, mappings),
     summary = summarize-mappings(python, tool, study, mapped),
     statistics = analyze-genes(python, tool, study, mapped),
     pathway-results = analyze-study-pathways(python, tool, study, pathways, statistics),
     report = render-study-report(python, tool, study, mapped, statistics, pathway-results),
     manifest = create-manifest(python, tool, study, expression, mappings, pathways, mapped, summary, statistics, pathway-results, report)) {
  {
    "mapping-receipts.tsv" -> mapped,
    "mapping-summary.json" -> summary,
    "gene-statistics.tsv" -> statistics,
    "pathway-results.tsv" -> pathway-results,
    "report.html" -> report,
    "provenance.json" -> manifest
  }
}
