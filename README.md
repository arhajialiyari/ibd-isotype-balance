# IBD immunoglobulin isotype balance

R workflow for **Cross-cohort evaluation of rectal immunoglobulin isotype balance in ulcerative colitis and Crohn disease**.

The workflow analyzes donor-level single-cell pseudobulk expression, participant-level bulk scores, cross-cohort discrimination and sensitivity analyses. It generates supplementary tables and computational figures in R.

## Data

| Dataset | Role | Source |
| --- | --- | --- |
| GSE282122 / TAURUS | Single-cell discovery | [Zenodo](https://zenodo.org/records/14007626) |
| GSE214695 | Single-cell candidate replication | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695) |
| GSE193677 | Exploratory bulk development | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193677) |
| GSE117993 | External score evaluation; post hoc analysis | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117993) |
| GSE111889 / HMP2 | Context cohort | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111889) |
| GSE57945 and GSE83687 | Score-component availability audits | [GSE57945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57945), [GSE83687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83687) |

Download expression files separately. Required filenames, relative paths, byte counts and MD5 values are listed in [config/inputs.tsv](config/inputs.tsv). Four source metadata files are bundled under `inputs/metadata/`; the HMP2 metadata snapshot comes from a [commit-pinned mirror](https://raw.githubusercontent.com/lotas-ai/hutlab-demo/9189026f2cb66883690c57ad73ad9d96b02a861c/data/hmp2_metadata_2018-08-20.csv).

Keep expression files outside the repository. Copy `config/local.example.R` to `config/local.R`, set an absolute input root and supply file overrides when needed:

```r
legacy_root <- "D:/IBD_inputs"
verify_md5 <- TRUE
file_overrides <- c(
  replication_annotations = "D:/IBD_inputs/GSE214695_cell_annotation.csv.gz"
)
```

Alternatively, set `IBD_DATA_ROOT` or pass `--source-root=PATH`. Bundled metadata are selected automatically unless explicitly overridden. The optional pairing CSV is not required for analysis.

## Environment setup

Use **R 4.4.2** and **Bioconductor 3.20**. Package versions are pinned in `renv.lock`. The commands below use Windows PowerShell; source compilation may require Rtools44.

From a fresh R session in the repository directory, configure the inputs above and run:

```r
stopifnot(getRversion() == package_version("4.4.2"))
bootstrap <- file.path(getwd(), "renv", "bootstrap")
dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
utils::install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.2.4.tar.gz",
  repos = NULL, type = "source", lib = bootstrap
)
utils::install.packages("zip", repos = "https://cloud.r-project.org", lib = bootstrap)
.libPaths(c(bootstrap, .Library), include.site = FALSE)

Sys.setenv(IBD_DATA_ROOT = "D:/IBD_inputs")
options(ibd.restore.fresh = TRUE, ibd.stage = "restore")
restore_status <- source("run.R", echo = FALSE)$value
stopifnot(identical(restore_status, 0L))
state <- readRDS("work/fresh_restore_state.rds")
stopifnot(isTRUE(state$validation_passed))
cat("Analysis directory:", state$request$project, "\n")
cat("Package library:", state$request$library, "\n")
```

Replace `D:/IBD_inputs` with your input root. This installs the locked dependencies in a separate library and prepares an analysis directory with configured input paths.

## Run

In PowerShell, substitute the analysis directory and package-library paths printed above. Use R 4.4.2's `Rscript` and execute one stage at a time:

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

Check `$LASTEXITCODE` and the current `work/*_report.txt` after each stage. Continue only when the stage succeeds; an absent optional pairing CSV can produce a preflight warning. Downstream stages require validated upstream results.

For interactive use, select the restored library, set `options(ibd.stage = "STAGE")`, and source `run.R`. Set the stage option for each invocation.

## Outputs

`results/Supplementary_Tables.xlsx` contains six sheets, `Table_S1` through `Table_S6`; S3A and S3B share `Table_S3` with separate schemas. `results/Manuscript_Tables.html` contains the two main tables. `results/figures/` contains Figures 2-5 and S1 as PNG and vector PDF files. Figure 1 is outside computational figure generation.

Working files and reports remain in `work/`. The outputs stage also creates `results/Publication_Outputs.zip`. Local configuration, installed libraries, working files and transport archives are excluded from Git.

## Analysis scope

The isotype score is `log2((G + 0.5)/(A + M + 0.5))`, using four IgG genes, two IgA genes and IGHM in native cohort units. Development and later comparator/subtype analyses retain their exploratory or post hoc roles. HMP2 is context only because of prior project exposure. GSE57945 and GSE83687 are nonestimable for the required score components.

Reference archives provide comparison targets. They do not supply fitted models or participant scores for the calculations. Discrimination, transcriptome context and component decomposition do not establish causal mechanisms or clinical utility.

## License

Original workflow code and project-authored documentation are distributed under the [MIT License](LICENSE). Third-party datasets, reference materials and dependencies retain their original attribution and terms.
