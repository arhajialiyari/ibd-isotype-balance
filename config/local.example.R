# Copy this file to local.R and set the existing project input directory.
# This local configuration and large expression inputs stay outside version control.
# Four exact public source metadata snapshots are bundled and selected by run.R.
legacy_root <- Sys.getenv("IBD_DATA_ROOT", unset = "")
verify_md5 <- TRUE
file_overrides <- character()
# Example override for an input moved outside the original directory:
# file_overrides <- c(replication_annotations = "D:/IBD_inputs/annotations.csv.gz")
