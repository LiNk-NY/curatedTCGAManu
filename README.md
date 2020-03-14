
# curatedTCGAManu

## Multi-omic integration of public oncology databases in Bioconductor

This repository contains scripts and datasets for the curatedTCGAData +
cBioPortalData manuscript.

GitHub Repository: <https://doi.org/XXX>

## Overview

We provide several examples to demonstrate the powerful and flexible
analysis environment provided. These analyses, previously only
achievable through a significant investment of time and bioinformatic
training, become straight- forward analysis exercises provided as
vignettes.

## Key Points

  - Key objective: To provide flexible, integrated multi-omic
    representations of public oncology databases in R/Bioconductor with
    grealy reduced data management overhead.

  - Knowledge generated: Our Bioconductor software packages provide a
    novel approach to lower barriers to analysis and tool development
    for the TCGA and cBioPortal databases.

  - Relevance: Our tools provide flexible, programmatic analysis of
    hundreds of fully integrated multi’omic oncology datasets within an
    ecosystem of multi-omic analysis tools.

## Key Packages

  - [MultiAssayExperiment](http://waldronlab.io/MultiAssayExperiment/index.html)
  - [curatedTCGAData](https://bioconductor.org/packages/curatedTCGAData)
  - [TCGAutils](https://bioconductor.org/packages/TCGAutils)
  - [cBioPortalData](https://github.com/waldronlab/cBioPortalData)
  - [EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
  - [GSEABenchmarkeR](https://bioconductor.org/packages/GSEABenchmarkeR)
  - [RTCGAToolbox](https://bioconductor.org/packages/RTCGAToolbox)

## Installation

Until the release of Bioconductor 3.11 (scheduled for April 28, 2020),
running the code in this repository requires the *devel* version of
Bioconductor. This can be installed the [traditional
way](https://bioconductor.org/developers/how-to/useDevel/) or by using
the [Docker container](https://bioconductor.org/help/docker/).

Additionally until the release of Bioconductor 3.11, cBioPortalData must
be installed from GitHub as shown in the following code chunk which
installs all necessary packages either directly or as dependencies. Note
that this code chunk is not evaluated, because installation only needs
to be performed once.

``` r
BiocManager::install(c("waldronlab/cBioPortalData", "LiNk-NY/curatedTCGAManu"))
```

## Loading packages

``` r
library(curatedTCGAData)
library(cBioPortalData)
library(TCGAutils)
library(GenomicDataCommons)
library(rtracklayer)
```

*Note*. For clarity, we include `library` commands within the
supplemental code chunks.

# Supplement Reference

## Figure S1

<a href="https://github.com/waldronlab/schematics/"><img src="https://raw.githubusercontent.com/waldronlab/schematics/master/pngs/TCGAMAEPipelineSchematic.png" width="989" height="1091"/></a>

## Figure S2A

### Example code for installing and downloading TCGA data using curatedTCGAData.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# BiocManager::install("curatedTCGAData")

library(curatedTCGAData)

## Glioblastoma Multiforme (GBM)

curatedTCGAData(diseaseCode = "GBM", assays = "*", dry.run = FALSE)
#> snapshotDate(): 2020-03-12
#> # [truncated...]
#> loading from cache
#> A MultiAssayExperiment object of 18 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 18:
#>  [1] GBM_CNACGH_CGH_hg_244a-20160128: RaggedExperiment with 81512 rows and 438 columns
#>  [2] GBM_CNACGH_CGH_hg_415k_g4124a-20160128: RaggedExperiment with 57975 rows and 338 columns
#>  [3] GBM_CNASNP-20160128: RaggedExperiment with 602338 rows and 1104 columns
#>  [4] GBM_CNVSNP-20160128: RaggedExperiment with 146852 rows and 1104 columns
#>  [5] GBM_GISTIC_AllByGene-20160128: SummarizedExperiment with 24776 rows and 577 columns
#>  [6] GBM_GISTIC_ThresholdedByGene-20160128: SummarizedExperiment with 24776 rows and 577 columns
#>  [7] GBM_miRNAArray-20160128: SummarizedExperiment with 534 rows and 565 columns
#>  [8] GBM_miRNASeqGene-20160128: SummarizedExperiment with 1046 rows and 0 columns
#>  [9] GBM_mRNAArray_huex-20160128: SummarizedExperiment with 18632 rows and 431 columns
#>  [10] GBM_mRNAArray_TX_g4502a_1-20160128: SummarizedExperiment with 17814 rows and 401 columns
#>  [11] GBM_mRNAArray_TX_g4502a-20160128: SummarizedExperiment with 17814 rows and 101 columns
#>  [12] GBM_mRNAArray_TX_ht_hg_u133a-20160128: SummarizedExperiment with 12042 rows and 528 columns
#>  [13] GBM_Mutation-20160128: RaggedExperiment with 22073 rows and 290 columns
#>  [14] GBM_RNASeq2GeneNorm-20160128: SummarizedExperiment with 20501 rows and 166 columns
#>  [15] GBM_RPPAArray-20160128: SummarizedExperiment with 208 rows and 244 columns
#>  [16] GBM_GISTIC_Peaks-20160128: RangedSummarizedExperiment with 68 rows and 577 columns
#>  [17] GBM_Methylation_methyl27-20160128: SummarizedExperiment with 27578 rows and 285 columns
#>  [18] GBM_Methylation_methyl450-20160128: SummarizedExperiment with 485577 rows and 154 columns
#> Features:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample availability DFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

## Figure S2B

### Example cBioPortalData code for downloading TCGA data from cBioPortal.org and via the cBioPortal API

``` r
## installation
# BiocManager::install("WaldronLab/cBioPortalData")

library(cBioPortalData)

## https://cbioportal.org/datasets (Bulk data method)
gbm <- cBioDataPack("gbm_tcga")
#> Study file in cache: gbm_tcga
#> # [truncated...]
#> harmonizing input:
#>   removing 7 colData rownames not in sampleMap 'primary'

## https://cBioPortal.org/api (API method)
cBio <- cBioPortal()
cBioPortalData(cBio, studyId = "gbm_tcga", genePanelId = "IMPACT341")
#> harmonizing input:
#>   removing 7 colData rownames not in sampleMap 'primary'
#> A MultiAssayExperiment object of 13 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 13:
#>  [1] gbm_tcga_rppa: SummarizedExperiment with 67 rows and 244 columns
#>  [2] gbm_tcga_rppa_Zscores: SummarizedExperiment with 67 rows and 244 columns
#>  [3] gbm_tcga_gistic: SummarizedExperiment with 339 rows and 577 columns
#>  [4] gbm_tcga_mrna_U133: SummarizedExperiment with 311 rows and 528 columns
#>  [5] gbm_tcga_mrna_U133_Zscores: SummarizedExperiment with 308 rows and 528 columns
#>  [6] gbm_tcga_mrna: SummarizedExperiment with 334 rows and 401 columns
#>  [7] gbm_tcga_mrna_median_Zscores: SummarizedExperiment with 329 rows and 401 columns
#>  [8] gbm_tcga_rna_seq_v2_mrna: SummarizedExperiment with 341 rows and 166 columns
#>  [9] gbm_tcga_rna_seq_v2_mrna_median_Zscores: SummarizedExperiment with 333 rows and 166 columns
#>  [10] gbm_tcga_linear_CNA: SummarizedExperiment with 339 rows and 577 columns
#>  [11] gbm_tcga_methylation_hm27: SummarizedExperiment with 282 rows and 285 columns
#>  [12] gbm_tcga_methylation_hm450: SummarizedExperiment with 282 rows and 153 columns
#>  [13] gbm_tcga_mutations: RangedSummarizedExperiment with 810 rows and 271 columns
#> Features:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DFrame
#>  sampleMap() - the sample availability DFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

## Figure S2C

### Example hg19 to hg38 liftover procedure using Bioconductor tools

``` r
liftchain <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
cloc38 <- file.path(tempdir(), gsub("\\.gz", "", basename(liftchain)))
dfile <- tempfile(fileext = ".gz")
download.file(liftchain, dfile)
R.utils::gunzip(dfile, destname = cloc38, remove = FALSE)

library(rtracklayer)
chain38 <- suppressMessages(
    import.chain(cloc38)
)
## Note. Run bulk data download example (in S2B) first
mutations <- gbm[["mutations_extended"]]
seqlevelsStyle(mutations) <- "UCSC"

ranges38 <- liftOver(rowRanges(mutations), chain38)
```

## Figure S3

### Example code for downloading data via GenomicDataCommons and loading with TCGAutils

``` r
library(TCGAutils)
library(GenomicDataCommons)

query <- files(legacy = TRUE) %>%
    filter( ~ cases.project.project_id == "TCGA-COAD" &
        data_category == "Gene expression" &
        data_type == "Exon quantification" )

fileids <- manifest(query)$id[1:4]
exonfiles <- gdcdata(fileids)

makeGRangesListFromExonFiles(exonfiles, nrows = 4)
#> Parsed with column specification:
#> cols(
#>   exon = col_character(),
#>   raw_counts = col_double(),
#>   median_length_normalized = col_double(),
#>   RPKM = col_double()
#> )
#> Parsed with column specification:
#> cols(
#>   exon = col_character(),
#>   raw_counts = col_double(),
#>   median_length_normalized = col_double(),
#>   RPKM = col_double()
#> )
#> Parsed with column specification:
#> cols(
#>   exon = col_character(),
#>   raw_counts = col_double(),
#>   median_length_normalized = col_double(),
#>   RPKM = col_double()
#> )
#> Parsed with column specification:
#> cols(
#>   exon = col_character(),
#>   raw_counts = col_double(),
#>   median_length_normalized = col_double(),
#>   RPKM = col_double()
#> )
#> GRangesList object of length 4:
#> $`TCGA-5M-AAT4-01A-11R-A41B-07`
#> GRanges object with 4 ranges and 3 metadata columns:
#>       seqnames      ranges strand | raw_counts median_length_normalized
#>          <Rle>   <IRanges>  <Rle> |  <numeric>                <numeric>
#>   [1]     chr1 11874-12227      + |          1                 0.135977
#>   [2]     chr1 12595-12721      + |          2                 0.547619
#>   [3]     chr1 12613-12721      + |          2                 0.472222
#>   [4]     chr1 12646-12697      + |          1                 0.529412
#>            RPKM
#>       <numeric>
#>   [1] 0.0228442
#>   [2] 0.1273517
#>   [3] 0.1483822
#>   [4] 0.1555160
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>
#> $`TCGA-A6-6782-01A-11R-1839-07`
#> GRanges object with 4 ranges and 3 metadata columns:
#>       seqnames      ranges strand | raw_counts median_length_normalized
#>          <Rle>   <IRanges>  <Rle> |  <numeric>                <numeric>
#>   [1]     chr1 11874-12227      + |         35                 0.784702
#>   [2]     chr1 12595-12721      + |          9                 0.873016
#>   [3]     chr1 12613-12721      + |          9                 0.851852
#>   [4]     chr1 12646-12697      + |          8                 0.843137
#>            RPKM
#>       <numeric>
#>   [1]  0.691243
#>   [2]  0.495456
#>   [3]  0.577274
#>   [4]  1.075605
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>
#> $`TCGA-AA-3678-01A-01R-0905-07`
#> GRanges object with 4 ranges and 3 metadata columns:
#>       seqnames      ranges strand | raw_counts median_length_normalized
#>          <Rle>   <IRanges>  <Rle> |  <numeric>                <numeric>
#>   [1]     chr1 11874-12227      + |          4                 0.492918
#>   [2]     chr1 12595-12721      + |          2                 0.341270
#>   [3]     chr1 12613-12721      + |          2                 0.398148
#>   [4]     chr1 12646-12697      + |          2                 0.372549
#>            RPKM
#>       <numeric>
#>   [1]  0.322477
#>   [2]  0.449436
#>   [3]  0.523655
#>   [4]  1.097661
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>
#> $`TCGA-AA-3955-01A-02R-1022-07`
#> GRanges object with 4 ranges and 3 metadata columns:
#>       seqnames      ranges strand | raw_counts median_length_normalized
#>          <Rle>   <IRanges>  <Rle> |  <numeric>                <numeric>
#>   [1]     chr1 11874-12227      + |          0                        0
#>   [2]     chr1 12595-12721      + |          0                        0
#>   [3]     chr1 12613-12721      + |          0                        0
#>   [4]     chr1 12646-12697      + |          0                        0
#>            RPKM
#>       <numeric>
#>   [1]         0
#>   [2]         0
#>   [3]         0
#>   [4]         0
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
