# Validation record

This record summarizes the author's completed R executions and subsequent independent review of returned files. Packaging was performed afterward and is not an additional R analysis. The tested platform was R 4.4.2 on Windows 11 x64, with Bioconductor 3.20 and the 19 versions in `renv.lock`. No execution on another operating system is claimed.

## Environment restoration

The locked packages were downloaded/installed into a separate private library. The final stage09.3 restoration/probe checked and retained all 19 exact-version installations from the completed remote attempt; it did not reinstall them or reuse the author's working library. The parent report had 11 PASS, 0 WARN and 0 FAIL; restore and independent probe workers had 2 and 23 PASS checks. All locked namespace paths, source/dependency records, HDF5 integer roundtrip, sparse Matrix operation and fixed-direction AUC fixture passed. Source, lock and package-file identities were checked.

## Complete isolated analysis

The stage10.0 run began at 2026-10-25 08:19:54 UTC and finished at 08:35:46 UTC: 952 seconds on the author's computer. Seven scientific stages and a comparison ran in independent fresh R processes. The analysis copy began with empty work/results directories. Native expression files were shared read-only as declared inputs; old fitted caches and publication products were not used for scientific computation. All five prerequisite caches used by output generation came from this run.

| Stage | Revision | PASS | WARN | FAIL |
| --- | --- | ---: | ---: | ---: |
| Preflight | stage01.3 | 34 | 0 | 0 |
| Discovery | stage02.0 | 20 | 0 | 0 |
| Replication | stage03.0 | 22 | 0 | 0 |
| Development | stage04.2 | 40 | 0 | 0 |
| External/context | stage05.0 | 42 | 0 | 0 |
| Sensitivity | stage06.0 | 30 | 0 | 0 |
| Publication outputs | stage07.3 | 38 | 0 | 0 |
| Independent R comparison | stage10.0 | 7 | 0 | 0 |

The parent verification had 14 PASS, 0 WARN and 0 FAIL. Each worker checked all 19 locked versions, source/dependency metadata and namespace paths. All 14 declared native input identities agreed after execution. Original author caches, products, code and package files were unchanged. The returned diagnostic archive contained 55 members; all 54 declared captured files were byte-verified. Raw diagnostics and machine-specific private state are not included here.

## Numerical and artifact review

The R comparison checked the newly calculated and previously reviewed workbooks against new full-precision values, including types, IDs, flags and missingness. Maximum absolute workbook roundtrip difference was `4.2632564e-14`. Independent read-only inspection found six sheets/seven named blocks, with zero old/new cell-value, cell-type or style differences. Workbook XML differed only in its creation timestamp. The two-main-table HTML was byte-identical.

All five PDFs had identical decompressed page drawing streams and extracted text to the accepted stage07.3 versions; each had one vector page and no out-of-page text. Raw PDF bytes can differ because of timestamps. Figures 2 and 3 had byte-identical PNGs. Figures 4, 5 and S1 differed in 76, 128 and 92 pixels, respectively, across images of millions of pixels. These three images were visually inspected and accepted with no material content or layout difference.

The original R report's automatic visual flag remains `REVIEW_REQUIRED` because PNG bytes differed. The subsequent visual assessment completes that review; the original report was not edited to manufacture an automatic PASS. This package retains original reviewed stage07.3 product bytes, identified in `config/reviewed_outputs.tsv`, rather than overwriting those identities with equivalent rerun files.

Ordinary reference comparisons use absolute tolerance `1e-10` plus relative tolerance `1e-6`; probabilities use `1e-14` absolute plus `1e-6` relative. IDs, counts and categories are exact. Published rounding is applied only to relevant comparison views. Workbook roundtrip tolerances are `1e-12` absolute plus `1e-12` relative, with `1e-14` absolute for probabilities. Stored calculations retain full precision.

## Boundaries

This establishes source-to-result reproduction on the tested Windows host with declared input bytes, available comparisons and a locked environment. It does not establish cross-platform execution, permanent remote availability, prospective specification, causal interpretation, clinical utility or equivalence from nonsignificance.

Historical exports are incomplete for several secondary comparisons. In particular, the old full 24-donor-deletion table is absent: checks cover the complete 2,101-gene baseline and manuscript-reported ranges, not unavailable deletion rows. Further interpretation limits remain in the README and scientific reports.

Figure 1 and both manuscript DOCX files are outside computational output generation. No manuscript editing, Git commit, upload or publication was performed by R validation or packaging. The package's seven R files match the successful stage10.0 source manifest; only documentation, Git file rules and administrative packaging records were added or updated afterward. The author selected MIT on 2026-10-25, after scientific validation; LICENSE and the distribution inventory were then updated without changing code or research products.
