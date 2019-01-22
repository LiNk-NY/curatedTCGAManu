## ExperimentList data types from curatedTCGAData metadata file
.splitCut <- function(strvec, splitter, position) {
    unlist(
        lapply(strsplit(strvec, splitter), `[`, position)
    )
}
meta <- system.file("extdata/metadata.csv", package = "curatedTCGAData",
    mustWork = TRUE)
metadata <- read.csv(meta, stringsAsFactors = FALSE)
datanames <- metadata[["Title"]]
datanames <- .splitCut(datanames, "-", 1L)
assaynames <- gsub("^[A-Z]*_(.*)", "\\1", datanames)
assaynames <- unique(assaynames)

assaynames[
    !assaynames %in% c("metadata", "sampleMap", "colData", "GISTIC", "CNACGH")]
