# Context — pqverse

Ubiquitous language for this project: the shared vocabulary used by the
developer and AI agents. Each entry maps a domain term or piece of
jargon to a plain definition, so conversations and code stay consistent
and concise.

Idea from <https://github.com/mattpocock/skills>.

## Glossary

**pqverse** — Monorepo workspace of R packages for metabarcoding
analysis, all centred on the `phyloseq` object. Named by analogy with
the tidyverse.

**phyloseq object** — The central S4 data structure (from the `phyloseq`
R package) grouping an OTU/ASV table, a taxonomy table, sample metadata,
and optionally a reference sequence set and a phylogenetic tree. All
pqverse packages manipulate or extend it.

**OTU (Operational Taxonomic Unit)** — A cluster of sequences grouped at
a similarity threshold (usually 97–99%) used as a proxy for species in
metabarcoding studies. Contrast with ASV.

**ASV (Amplicon Sequence Variant)** — An exact, denoised representative
sequence produced by DADA2 or similar tools; finer resolution than an
OTU and fully reproducible across studies.

**metabarcoding** — High-throughput DNA sequencing of a taxonomic marker
gene (e.g. ITS for fungi, 16S for bacteria) from an environmental sample
to characterise community composition.

**eDNA (environmental DNA)** — DNA extracted directly from an
environmental sample (water, soil, air) without isolating individual
organisms. The umbrella term covering metabarcoding and related amplicon
approaches.

**ITS (Internal Transcribed Spacer)** — The standard DNA barcode region
for fungi, located between ribosomal RNA genes. ITS1 and ITS2
sub-regions are amplified with primers such as ITS1-F / ITS2.

**MiscMetabar** — The flagship package of the pqverse, providing
miscellaneous helper functions for description, transformation,
exploration, and reproducibility of metabarcoding analyses. Published on
CRAN and JOSS.

**tidypq** — pqverse package offering tidyverse-style verbs
(`filter_samples_pq`, `mutate_taxa_pq`, …) for manipulating phyloseq
objects at four scales: samples, taxa, occurrences, and tree.

**contaminant** — Broad umbrella for any taxon flagged for removal
because it does not represent genuine biological signal in the sample.
Sub-types: *technical* (chimeras, primer read-through), *reagent/lab*
(`lab_contaminant`), *cross-sample* (`sample_contaminant`), and
`artifact`. The narrow sense (lab/sample/artifact, as classified by
`neg_control_classify_pq()`) is one branch of this umbrella.

**detector** — A contamination-detection method in tidypq that flags
suspect taxa without removing them. Pure: it returns a `contam_tbl` and
never mutates the phyloseq object; removal is a separate step (the
filter seam).

**contam_tbl** — Canonical object returned by every tidypq detector: a
class-tagged tibble with one row per flagged taxon, keyed by a `taxon`
column (matching `taxa_names()`), plus a `method` column and
method-specific evidence columns. Row-bindable across detectors, and
consumed both by the removal seam (`filter_contam_pq()`) and by
contamination diagnostics.

**dbpq** — pqverse package for downloading, formatting, and filtering
FASTA reference databases used in taxonomic assignment (UNITE,
Eukaryome, SILVA, …).

**comparpq** — pqverse package for comparing alternative representations
of a metabarcoding community, across two axes: alternative taxonomic
assignations within one phyloseq (`tc_` family — `tc_metrics_mock()`,
`tc_heatmap()`, `tc_circle()`, `tc_linked_trees()`), or alternative
phyloseq objects (`list_phyloseq` S7 class + `_lpq` family). The
unifying question is “how do ≥2 views of the same community differ?”;
functions that operate on a single view without comparing it (single-pq
differential abundance, bubble plots, primer detection, mock-community
creation, taxonomy-tree plotting) belong elsewhere. *Avoid*: treating
“single vs multiple objects” as the boundary — the boundary is “compares
≥2 views.”

**view** — One representation of a metabarcoding community being
compared in comparpq. Either an *assignation view* (a `tax_table` rank
column, e.g. `Genus_SILVA` vs `Genus_KSGP`) or an *object view* (a
phyloseq object held in a `list_phyloseq`). comparpq functions must
compare ≥2 views; single-view functions are out of scope.

**taxinfo** — pqverse package for augmenting phyloseq objects with
external taxonomy-based information from GBIF, Wikipedia, and GloBI.

**`augment_tax_table()`** — taxinfo-internal deep module (the merge-back
seam) shared by the `tax_*_pq` family. Takes a phyloseq object and a
per-taxon `info_tbl` keyed by a column (`info_key`, default
`taxa_name`), and returns the phyloseq with the tibble’s columns
left-joined into its `tax_table`. Owns key construction, `col_prefix`
collision handling, the `left_join`, the `as.matrix() |> tax_table()`
round-trip, and the `rownames <- taxa_names()` restore. Each `tax_*_pq`
function only builds its `info_tbl`; the external fetch stays in the
caller.

**`taxnames_from_rank()`** — taxinfo-internal key-builder shared by
`taxonomic_rank_to_taxnames()` (query side) and `augment_tax_table()`
(merge side). Pastes the `taxonomic_rank` column(s) per taxon into a
single name and applies the `"NA NA"` / `" NA"` cleanup, guaranteeing
the API query key and the merge-back join key are identical (closing a
latent silent-drop bug for multi-column ranks with genus-only taxa).

**`resolve_taxa_input()`** — taxinfo-internal front-matter helper shared
by the `tax_*_pq` family (the counterpart of the `augment_tax_table()`
merge-back). Validates the mutually exclusive `physeq` / `taxnames`
input, resolves the `add_to_phyloseq` default, and extracts taxon names
from `physeq` via `taxonomic_rank_to_taxnames()`; returns
`list(taxnames, add_to_phyloseq)`. Each `tax_*_pq` function calls it
instead of repeating the ~22-line validate-and-extract block.

**greenAlgoR** — pqverse package for estimating the carbon footprint of
R computations, including targets pipelines, based on the Green
Algorithms framework.

**ggplotpq** — pqverse package providing ggplot2-based visualisation
functions for phyloseq objects. Separates data grouping (`fact`) from
colour aesthetics (`color_by`) and faceting (`facet_by`), unlike
MiscMetabar where these overlap.

**`fact` / `bifactor`** — Canonical parameter names across pqverse
functions: `fact` = primary grouping column in `sample_data()`,
`bifactor` = secondary grouping column. Variants such as `modality`,
`variable`, `var_to_test`, and `by.fact` are being unified toward these
names.

**`merge_sample_by`** — Canonical parameter name for the sample-data
column used to aggregate (merge) samples within a phyloseq object.
Already consistent across MiscMetabar, comparpq, and ggplotpq.

**`pq_to_tidy()`** — Canonical phyloseq-to-tidy data-preparation module,
exported from **tidypq**. Converts a phyloseq object into an ungrouped
tibble with columns `sample_id`, `taxon_id`, `abundance`, `fact`,
`bifactor`, selected taxonomic rank columns (wide, NA→“Unknown”), and
all `sample_data` columns. Owns the full pipeline: `verify_pq` →
`clean_pq` → `pivot_longer` → join sample_data + tax_table → optional
`merge_sample_by` aggregation (tibble space) → optional `transform` (per
sample_id, keeps `abundance_raw`) → resolve fact/bifactor (enforces
2-level on bifactor) → optional `filter_zero`. Does NOT wrap
[`phyloseq::psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html) —
it builds from scratch. Does NOT replicate
[`MiscMetabar::psmelt_samples_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/psmelt_samples_pq.md)
(sample-level Hill summary). See ADR 0002. **Prefer `pq_to_tidy()` over
reinventing the pipeline in each plot/analysis function.**

**`_pq` suffix** — Naming convention for all exported pqverse functions
that take a phyloseq object as their primary argument
(e.g. `alpha_div_pq()`, `filter_samples_pq()`).

**pqverse meta-package** — A thin wrapper package (following the
tidyverse pattern) that installs and loads all pqverse packages with a
single [`library(pqverse)`](https://adrientaudiere.github.io/pqverse/)
call.

**`build_all_packages.R`** — The pqverse build engine: runs
`air format → build_readme → document → test → build_news → build_site → build → codecov`
for each package and writes CSV metrics and an issues log under
`pqverse_metrics/`.

**pqverse_metrics/** — Directory under the workspace root collecting
per-run and latest-run CSVs from `build_all_packages.R`, build-issue
logs, and the `report/` sub-directory with PNG figures and `report.html`
produced by `pqverse_build_report.R`.

**merge-release.sh** — Workspace-level script that merges the `dev`
branch into `main`, bumps the version in `DESCRIPTION` and `NEWS.md`,
and pushes. Accepts `--bump major|minor|patch` (default: `minor`) and
`--dry-run`.

**ROADMAP_kanban.html** — Auto-generated static kanban board from
`ROADMAP.md`, served locally by `scripts/run_kanban.py`. Items are
tagged with priority (Critical/High/Medium/Low) and facility
(easy/moderate/hard).

**`\donttest{}` / `\dontrun{}`** — roxygen2 example-wrapping strategies.
`\donttest{}` runs in CRAN’s `--run-donttest` pass (use for primary
examples); `\dontrun{}` never runs automatically (use to preserve extra
examples without contributing to check time). Never delete examples —
wrap them instead.

## Key decisions

**phyloseq as the central data structure** — Every pqverse package
builds on the `phyloseq` S4 class rather than introducing a new class.
This keeps the ecosystem interoperable with the large existing phyloseq
user base and avoids the cost of an S4/S7 migration.

**MiscMetabar scope guard** — MiscMetabar stays lean: new features that
require a heavy or brand-new dependency are dispatched to a specialist
sub-package (comparpq, taxinfo, dbpq, greenAlgoR, tidypq, ggplotpq). The
pqverse meta-package absorbs growth instead of MiscMetabar.

**`Depends` vs `Imports` in MiscMetabar** — Heavy packages (notably
`dada2`) live in `Imports`, not `Depends`, to prevent them from being
attached (and paying their load cost) at every CRAN check stage.
`phyloseq`, `ggplot2`, and `dplyr` remain in `Depends` for now; moving
them to `Imports` is deferred to 0.17.0 because it is a breaking change
for users.

**CRAN check-time budget** — MiscMetabar was archived on 2026-05-19 for
exceeding CRAN’s 10-minute check limit. All new examples must be wrapped
in `\donttest{}` or `\dontrun{}`; `data_fungi_mini` examples must be
restricted to 5 samples; `skip_on_cran()` must appear as the very first
line of test files to skip them on CRAN.

**Consolidated build/release workflow** — The pqverse uses two skills
(`/pqverse-build`, `/pqverse-fix`) plus `build_all_packages.R` and
`merge-release.sh` rather than the older per-package `r-pkg-release` /
`r-pkg-bump-version` scripts. The latter are retained only for
`winbuilder`/`rhub` remote pre-CRAN checks.

**Latest-CSV merge behaviour** — `build_all_packages.R` run on a subset
of packages carries forward previous rows for the omitted packages in
the `_latest.csv` files, so `pqverse_build_report.R --latest` always
renders the full picture. The per-run timestamped CSVs reflect only the
packages actually built in that run.
