# ppath

`ppath` is a small bioinformatics project about a step that is easy to overlook:
figuring out what our gene identifiers actually mean before we use them to make
biological claims.

I built it while exploring [`pp`](https://github.com/pranavra0/pp), a language
created by my friend. I wanted to use `pp` for something connected to the work I
already care about—bioinformatics, gene-expression analysis, and biological
pathways—instead of making a generic programming demo.

This is a for-fun research-software experiment. It mainly explores mapping
transparency, provenance, caching, and what a bioinformatics workflow can look
like in `pp`.

This repository is the result: an organism-aware gene-to-pathway workflow where
every mapping decision stays visible and every result can be traced back to the
files and code that produced it.

## Why I made this

A gene can be represented by a symbol, an Entrez ID, a KEGG identifier, a locus
tag, or another database-specific name. Translating between those systems is
not always as simple as looking up one value:

- an identifier may not exist for the selected organism;
- one identifier may point to several possible genes;
- two input names may resolve to the same gene;
- annotations can change between database releases;
- some genes can disappear silently before pathway analysis begins.

If those decisions are hidden inside a notebook, it can be difficult to explain
why a gene was included, why another was dropped, or why the result changed six
months later.

`ppath` treats mapping as part of the scientific evidence. Instead of quietly
discarding problems, it creates a receipt for every input identifier. The
receipt records what it mapped to, where the mapping came from, whether it was
ambiguous, and whether it was used in the analysis.

## What the project does

The demonstration follows a complete, intentionally small workflow:

```text
gene-expression table + organism-specific mapping snapshot
  -> mapping receipts
  -> duplicate-gene collapse
  -> gene-level fold changes
  -> pathway over-representation analysis
  -> multiple-testing correction
  -> readable HTML report
  -> SHA-256 provenance manifest
```

It starts with gene-expression measurements from two groups. Each gene is
resolved using a mapping table for the selected organism. Valid mappings move
forward, while ambiguous and unmapped identifiers remain visible in the audit.

When several input names resolve to the same Entrez gene, `ppath` combines them
before calculating statistics and records which original names were collapsed.
It then calculates group means and log2 fold changes, identifies genes that meet
the study's change threshold, and tests whether changed genes are
over-represented in each pathway. Pathway p-values are corrected with the
Benjamini–Hochberg procedure.

The final HTML report brings the mapping audit, gene-level results, and pathway
results together in one place.

## Why I used `pp`

I could have written the entire project as one Python script, but the point was
to find a meaningful use for a language my friend created. I wanted to see what
`pp` would look like in the kind of bioinformatics work I already enjoy rather
than using it for an example with no connection to my interests.

Bioinformatics workflows also happen to fit the ideas behind `pp`. An analysis
usually has several connected stages, and changing one input should not require
rerunning every unrelated piece. It is also useful to know exactly which code
and data produced a result. Because `pp` is content-addressed, it treats a
computation's code and inputs as part of that computation's identity.

For this project, that means `pp` can:

- remember work that has already been completed;
- reuse a result when its code and inputs have not changed;
- run independent parts of the workflow in parallel;
- require explicit permission before running programs or writing files;
- reconstruct generated files from stored results;
- explain with `pp why` why something was reused or calculated again.

Python still handles the tables and statistics because it is a natural tool for
that job. `pp` connects those Python operations into a reproducible workflow.
The mapping, gene analysis, pathway analysis, report, and provenance manifest
are separate cached computations, so repeating an unchanged study starts no new
analysis processes.

## What is included

The repository contains a small educational human dataset with normal and
primary-tumor-style sample groups. It is designed to exercise important edge
cases:

- exact mappings;
- an alias mapping;
- an identifier with multiple possible mappings;
- an unmapped identifier;
- two input identifiers that collapse to the same Entrez ID;
- pathways with complete, partial, and zero measured coverage.

This is curated demonstration data, not patient data or a clinical dataset. The
results are for education and software testing only and must not be used for
diagnosis or treatment decisions.

## Quick start

The Python reference workflow has no third-party Python dependencies. It needs
Python 3.11 or newer:

```sh
python3 scripts/run_local.py
```

Run the tests with:

```sh
python3 -m unittest discover -s tests -v
```

The generated report will be available at `build/report.html`.

## Run it through `pp`

First build the `pp` interpreter. On macOS, the required packages can be
installed with Homebrew:

```sh
brew install opam pkgconf zlib
scripts/bootstrap_pp.sh
```

Then run the content-addressed workflow:

```sh
scripts/run_pp.sh
```

The first run calculates and materializes the results. Repeating the command
should report zero created, updated, or deleted files because the desired output
already matches the build directory.

## Generated files

The workflow creates six artifacts:

```text
build/
├── mapping-receipts.tsv   # one mapping decision for every input identifier
├── mapping-summary.json   # counts of exact, alias, ambiguous, and failed mappings
├── gene-statistics.tsv    # collapsed genes, group means, and log2 fold changes
├── pathway-results.tsv    # pathway coverage, tests, and corrected p-values
├── report.html            # human-readable study report
└── provenance.json        # hashes of the study, code, inputs, and outputs
```

The provenance manifest does not merely list filenames. It records a SHA-256
hash and byte count for each artifact, the mapping-snapshot label, the organism,
the study ID, and the version of `ppath` that produced it.

## Project layout

```text
data/                       demonstration expression, mapping, and pathway data
studies/                    study configuration and analysis thresholds
src/ppath.py                mapping, statistics, pathway, report, and manifest code
ppath.pp                    content-addressed workflow written in pp
scripts/run_local.py        direct Python reference runner
scripts/run_pp.sh           pp workflow runner
scripts/bootstrap_pp.sh     local pp interpreter setup
tests/test_ppath.py         behavior and determinism tests
```

## Using different data

The demonstration can be replaced with another study by editing the files in
`data/` and the configuration in `studies/gse50760-demo.json`.

The expression table needs one gene-identifier column and one column for each
sample. The study configuration names the control and case columns.

The mapping table uses these columns:

```text
organism  input_id  entrez_id  symbol  mapping_status  source
```

The pathway table uses these columns:

```text
organism  pathway_id  pathway_name  entrez_ids
```

`entrez_ids` contains semicolon-separated identifiers.

Mapping and pathway database downloads are intentionally kept outside the
cached analysis. A real study should download and review a database snapshot,
save it as a versioned input, and record its name in the study configuration.
This prevents an invisible live-database update from changing an old result.

## Current scope

This repository is a working proof of concept, not a replacement for mature
workflow systems or established differential-expression packages. The current
statistics are appropriate for demonstrating provenance and workflow behavior,
but a real expression study should connect the workflow to tools such as
DESeq2, edgeR, or Pathview and use a properly designed biological dataset.

The main contribution here is the transparent chain from identifier mapping to
pathway result: nothing is silently guessed, discarded, or separated from its
provenance.

## License

MIT
