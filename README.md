# IBD immunoglobulin isotype balance

R workflow accompanying **Cross-cohort evaluation of rectal immunoglobulin isotype balance in ulcerative colitis and Crohn disease**.

Repository: <https://github.com/arhajialiyari/ibd-isotype-balance>

The workflow reconstructs donor-level single-cell analyses, bulk participant scores, external/context analyses, sensitivities, supplementary tables and computational figures from the declared native expression inputs. All scientific calculations, workbook writing and computational figure rendering are in R. Figure 1 is the separately authored AI-generated graphical abstract and is outside this workflow.

## Validation status

The complete source-to-result workflow passed on the author's Windows 11 computer with R 4.4.2 and Bioconductor 3.20 on 2026-10-25. Seven analysis stages and a final comparison ran in separate fresh R processes using an independently restored, exact-version 19-package library. Analysis directories began empty; previous fitted caches and publication products were not analysis inputs. The isolated run took 15 minutes 52 seconds on that computer. This is an observed runtime, not a hardware-independent estimate.

Every stage and comparison reported zero warnings/failures. The returned workbook and figures were independently inspected. See [VALIDATION.md](VALIDATION.md) for check counts, comparison scope, visual review and limitations. Other operating systems have not been tested.

This package preserves the seven R files from the successful stage10.0 run. Documentation, Git file rules and packaging records were updated afterward; no statistical implementation or R-generated publication product changed. Included products retain the reviewed stage07.3 byte identities in `config/reviewed_outputs.tsv`. The isolated rerun reproduced their contents; timestamps and a few PNG rendering pixels can differ without a numerical change.

## Files and scope

`run.R` is the entry point. Six modules in `R/` implement discovery, replication, development, external/context analyses, sensitivity analyses and outputs. Public configuration is in `config/`. Four exact source metadata snapshots are in `inputs/metadata/`. Six archives in `reference/` contain historical comparison targets, opened after computation rather than used to select or calculate results.

`renv.lock` pins 19 packages. `software_versions.tsv` records versions/sources; `source_identity.tsv` records native input identities. `release_file_manifest.tsv` records package member sizes, MD5 and SHA-256, excluding itself. `.gitattributes` disables Git line-ending conversion to preserve recorded source/metadata bytes.

Large expression matrices, private `config/local.R`, installed package libraries, diagnostic ZIPs and fitted RDS caches are excluded. Bundled source metadata and historical reference archives do contain public sample identifiers and derived comparison rows; this is not an assertion that the repository contains no sample-level information.

## Obtain and configure inputs

Use the exact files, byte counts and MD5 values in `config/inputs.tsv`. Four metadata files are supplied and automatically resolved through `config/bundled_inputs.tsv`. Explicit user file overrides take precedence. Analysis commands do not download or modify native inputs.

| Dataset | External native expression files | Recorded source |
| --- | --- | --- |
| GSE282122 / TAURUS | `myeloid_final.h5ad`, `fibperi_final.h5ad` | [Zenodo 14007626](https://zenodo.org/records/14007626) |
| GSE214695 | `GSE214695_RAW.tar`, `GSE214695_cell_annotation.csv.gz` | [GEO GSE214695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695) |
| GSE193677 | `GSE193677_MSCCR_Biopsy_counts.txt.gz` | [GEO GSE193677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193677) |
| GSE117993 | `GSE117993_RAW.tar` | [GEO GSE117993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117993) |
| GSE111889 / HMP2 | `GSE111889_host_tx_counts.tsv.gz` | [GEO GSE111889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111889) |
| GSE57945 | `GSE57945_all_samples_RPKM.txt.gz` | [GEO GSE57945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57945) |
| GSE83687 | `GSE83687_RAW.tar` | [GEO GSE83687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83687) |

These are recorded acquisition sources. The Zenodo record currently lists both required H5AD files and optional `paired_sample_list.csv` with MD5 values matching the manifest. All large files have not been downloaded again for packaging; availability and byte identity are not guaranteed solely by a landing-page link. Match the manifest before analysis. The optional pairing CSV is inventoried for provenance, while donor identities come from native H5AD metadata. Its absence produces a preflight warning.

The bundled HMP2 snapshot originated from this [commit-pinned mirror](https://raw.githubusercontent.com/lotas-ai/hutlab-demo/9189026f2cb66883690c57ad73ad9d96b02a861c/data/hmp2_metadata_2018-08-20.csv), with archived MD5 `465421f86e6eb104c9f98675ab209a35`. Bundling preserves original bytes rather than substituting a newly formatted export. Original dataset attribution and source terms remain applicable.

Choose an absolute input root and arrange external files under the relative paths in `config/inputs.tsv`. Alternatively, copy `config/local.example.R` to the ignored `config/local.R` and specify absolute `file_overrides` keyed by input ID:

```r
legacy_root <- "D:/IBD_inputs"
verify_md5 <- TRUE
file_overrides <- c(
  replication_annotations = "D:/IBD_inputs/GSE214695_cell_annotation.csv.gz"
)
```

Root precedence is `--source-root=PATH`, then `IBD_DATA_ROOT`, then `legacy_root` in `config/local.R`. Bundled metadata paths are automatic unless explicitly overridden. Keep external expression files outside the cloned repository.

## Restore a separate environment

Use R 4.4.2. Source builds can require Rtools44 on Windows; the completed restoration compiled packages including stringi and Matrix. Internet access, temporary storage and the relevant build tools are necessary. Linux/macOS build prerequisites and execution have not been validated here.

There is no automatic `.Rprofile` activation. The following explicit setup requests a new private restoration through the existing R administrative action. It does not require the author's old work directory or caches. The bootstrap version is available from the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/renv/); the supplied lock controls the separate analysis target.

In a fresh R 4.4.2 session, set the working directory to the clone, configure the input paths, then run:

```r
stopifnot(getRversion() == package_version("4.4.2"))
project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
bootstrap <- file.path(project, "renv", "bootstrap")
dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
utils::install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.2.4.tar.gz",
  repos = NULL, type = "source", lib = bootstrap
)
utils::install.packages("zip", repos = "https://cloud.r-project.org", lib = bootstrap)
.libPaths(c(bootstrap, .Library), include.site = FALSE)
stopifnot(utils::packageVersion("renv") == package_version("1.2.4"))

Sys.setenv(IBD_DATA_ROOT = "D:/IBD_inputs")
options(ibd.restore.fresh = TRUE, ibd.stage = "restore")
restore_status <- source("run.R", echo = FALSE)$value
stopifnot(identical(restore_status, 0L))
state <- readRDS("work/fresh_restore_state.rds")
stopifnot(isTRUE(state$validation_passed))
cat("Analysis directory:", state$request$project, "\n")
cat("Package library:", state$request$library, "\n")
```

Replace `D:/IBD_inputs` with the configured root. Bootstrap installation is separate from the locked target. That target must physically contain all 19 exact-version packages. The restore action uses a private target, completes missing locked CRAN packages without upgrades, and probes namespaces/source records and compiled runtime fixtures in a second clean process. This includes Matrix and lattice even when they also occur in R's system library. An unresolved source/version mismatch stops restoration. Network installation can take substantially longer than the subsequent analysis.

Require the current `work/restore_report.txt` to show `Fresh dependency restoration: PASS` and zero warnings/failures. A failed attempt cannot certify an older saved state. Keep the restoration directory. Its printed analysis directory contains copied code, private pointers to the declared inputs, and empty work/results directories.

## Run the seven scientific stages

Use R 4.4.2's `Rscript` from PowerShell. Replace the two paths below with the analysis directory and package library printed by the successful restore. The new empty analysis copy is used, preserving the reviewed publication snapshots at the repository root.

```powershell
Set-Location "D:/path/to/work/fresh_restore_ID/project"
$env:R_LIBS = ""
$env:R_LIBS_SITE = ""
$env:R_LIBS_USER = "D:/path/to/work/fresh_restore_ID/package_library"

Rscript --vanilla run.R --stage=preflight
Rscript --vanilla run.R --stage=discovery
Rscript --vanilla run.R --stage=replication
Rscript --vanilla run.R --stage=development
Rscript --vanilla run.R --stage=external
Rscript --vanilla run.R --stage=sensitivity
Rscript --vanilla run.R --stage=outputs
```

Execute one stage at a time. After each command require `$LASTEXITCODE` to be zero and inspect its current `work/*_report.txt` before continuing. An absent optional pairing file can produce a documented preflight warning; other warnings require investigation. Stages cover input readiness, single-cell discovery, single-cell candidate replication, bulk development, external/context reconstruction, sensitivities and output generation, respectively.

Downstream stages accept only validated calculated upstream caches with a consistent dependency chain. A fatal run does not refresh a preceding successful cache: an old RDS alone is not evidence that the current attempt passed. Checksum skipping is refused for analysis. Historical fitted objects and an undeclared original output directory are not required.

For interactive use, select the same restored library explicitly, set `options(ibd.stage = "STAGE")`, and source the copied `run.R`. The stage option is consumed each time; omitting it on a later invocation defaults to preflight.

The `release` and `verify` actions retain the author's guarded packaging and returned-evidence audit protocols. Their historical snapshot reports are not the current validation summary; `verify` expects specific reviewed author evidence. These actions are not prerequisites for the seven-stage route above. No command commits or publishes to GitHub.

## Publication products

`results/Supplementary_Tables.xlsx` has six ordered sheets, `Table_S1` through `Table_S6`, containing seven blocks. S3A and S3B share `Table_S3` with independent schemas. Saved values, types, IDs, flags and missingness are checked against computed values. Display formatting does not replace stored precision.

`results/Manuscript_Tables.html` contains two main tables for review/transfer rather than journal typesetting. `results/figures/` contains Figures 2, 3, 4, 5 and S1, each as a 300-dpi PNG and vector PDF. Figure 1 is not generated, copied or modified. Neither manuscript DOCX nor supplementary DOCX is edited by these scripts.

The outputs stage also creates `results/Publication_Outputs.zip` for transport. This duplicate archive, working reports, installed libraries and RDS caches are ignored by Git and excluded from this package. All 12 individual publication products were generated by the author's R execution and copied without alteration.

## Scientific interpretation

The seven-gene isotype score is `log2((G + 0.5)/(A + M + 0.5))`: four IgG genes in G, two IgA genes in A, and IGHM in M. Ratios use native cohort units; abundance transforms are cohort-specific. Development GSE193677 pools 679 biopsies into 665 participants. External GSE117993 includes 43 UC, 92 CD and 55 controls; controls are excluded from UC-versus-CD AUCs. HMP2/GSE111889 pools 100 rectal samples into 90 participants and remains context only because of prior project exposure.

Score exploration, later comparator additions and external subtype analyses retain exploratory/post hoc status. IGHG3 passes the implemented three-gene macrophage BH family but not pooled adjustment across all 13 candidates. Historical evidence does not independently establish prospective specification of those families. GSE83687 and GSE57945 are nonestimable for the required components; no substitute-score AUC is assigned.

The 18,940-gene AUC background uses observed labels. It is descriptive, not a permutation null or correction for score exploration. DeLong intervals are pointwise; paired comparisons use the declared retrospective Holm families. Nonsignificance does not establish equivalence. Component decomposition is descriptive algebra, not causal or AUC decomposition. Cell co-detection and donor influence do not establish macrophage-intrinsic transcription.

Historical exports are incomplete for the full replication tested-gene universe, every background-gene AUC, donor-deletion rows and pseudocount analyses. The new 24-deletion analysis is compared with its 2,101-gene baseline and manuscript-reported ranges; full old row-by-row agreement is not claimed where no old export exists. Numerical reproduction does not remove these evidence limits or establish clinical utility.

## License and citation

The original workflow code and project-authored documentation are licensed under the [MIT License](LICENSE), selected by the author on 2026-10-25. The license notice identifies the IBD Isotype Balance contributors. Bundled third-party source metadata, historical reference materials and package dependencies retain their original attribution and terms; the software license does not relicense those materials or assign a manuscript publication license to the research tables and figures.

The manuscript title above identifies the associated study. Final citation metadata remain pending author confirmation; no article DOI, release DOI or software version tag is asserted. This package was prepared before the first repository upload.
