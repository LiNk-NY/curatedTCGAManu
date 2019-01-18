## Code to obtain hit matrix of available subtypes by cancer
cname <- names(subtypeMaps)

allTypes <- unique(unlist(lapply(subtypeMaps, `[`, 1L)))
allTypes <- allTypes[allTypes != "Patient_ID"]

rnametypes <- gsub("_subtype[s]*", "", allTypes)
rnametypes <- gsub("copynumber", "Copy Number", rnametypes)
rnametypes <- gsub("_", " ", rnametypes)
toCap <- c("rna", "scna", "msi", "cimp", "rppa")
for (i in seq_along(toCap))
    rnametypes<- gsub(toCap[i], toupper(toCap[i]), rnametypes)

firstLetter <- function(charvect) {
    vapply(charvect, function(x) {
        if (x %in% c("mRNA", "microRNA")) return(x)
        s <- strsplit(x, " ")[[1]]
        paste0(toupper(substring(s, 1,1)), substring(s, 2), collapse=" ")
    }, character(1L))
}

rnametypes <- firstLetter(rnametypes)


hitMat <- vapply(subtypeMaps,
    function(x) allTypes %in% x[[1L]], logical(length(allTypes)))
rownames(hitMat) <- rnametypes

subtypesAvail <- ifelse(hitMat, "X", "")

subtypesAvail

write.csv(subtypesAvail, "data/Available_Subtypes_BY_CA.csv")


