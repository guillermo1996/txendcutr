
<!-- README.md is generated from README.Rmd. Please edit that file -->

# txendcutr

<!-- badges: start -->

[![R build
status](https://github.com/guillermo1996/txendcutr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/guillermo1996/txendcutr/actions)
<!-- badges: end -->

## Overview

Various mRNA sequencing library preparation methods generate sequencing
reads from the transcript ends. Quantification of isoform usage can be
improved by using truncated versions of transcriptome annotations when
assigning such reads to isoforms. The `txendcutr` package implements
some convenience methods for readily generating such truncated
annotations from either their 5’ or 3’ transcript ends and their
corresponding sequences.

`txendcutr` is a standalone fork of
[`txcutr`](https://github.com/mfansler/txcutr), maintained independently
of the upstream Bioconductor package. See [Differences from
`txcutr`](#differences-from-txcutr) below for what has changed.

## Installation instructions

`txendcutr` is not published on Bioconductor or CRAN. Install it
directly from GitHub. Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/), then install `txendcutr` with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## BiocManager::install() also installs from GitHub, and correctly
## resolves this package's Bioconductor dependencies
BiocManager::install("guillermo1996/txendcutr")
```

or, equivalently, with `remotes`:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("guillermo1996/txendcutr")
```

## Example

A typical workflow for `txendcutr` involves

- loading an existing annotation as a `TxDb` object
- truncating the annotation from the 3’ or 5’ end
- exporting the truncated annotation (GTF)
- exporting supporting files (FASTA, merge TSV)

``` r
library(rtracklayer)
library(txendcutr)
library(BSgenome.Hsapiens.UCSC.hg38)

## load human genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## load human GENCODE annotation
txdb <- makeTxDbFromGFF("gencode.v38.annotaton.gtf.gz", organism="Homo sapiens")

## truncate to maximum of 500 nts from the 3' end
txdb_w500 <- truncate3primeTxome(txdb, maxTxLength=500)

## ...or from the 5' end
txdb_5p_w500 <- truncate5primeTxome(txdb, maxTxLength=500)

## export annotation
exportGTF(txdb_w500, file="gencode.v38.txendcutr_w500.gtf.gz")

## export FASTA
exportFASTA(txdb_w500, genome=hg38, file="gencode.v38.txendcutr_w500.fa.gz")

## export merge-table
exportMergeTable(txdb_w500, minDistance=200,
                 file="gencode.v38.txendcutr_w500.merge.tsv.gz")
```

## Differences from `txcutr`

`txendcutr` began as a fork of
[`mfansler/txcutr`](https://github.com/mfansler/txcutr) and has since
diverged as a standalone package. Key changes:

- **Native 5’ truncation.** Upstream `txcutr` only truncates from the 3’
  end. `txendcutr` generalizes `truncateTxome()` with a `txEnd` argument
  (`"3prime"` or `"5prime"`), and adds `truncate3primeTxome()` /
  `truncate5primeTxome()` convenience wrappers (each with a `quiet`
  option to suppress progress messages).
- **Overlap export.** `truncateTxome()` (and both wrappers) accept an
  `overlapFile` argument to export a TSV of transcript pairs that were
  collapsed as duplicates during truncation, for auditing and debugging.
- **Faster truncation pipeline.** The internal clipping step no longer
  schedules one `BiocParallel` task per transcript; transcripts are now
  batched across workers and processed with a vectorized `dplyr`
  pipeline, which is substantially faster on large transcriptomes.

## Citation

`txendcutr` is a derivative work, built directly on the methodology and
original implementation of
[`txcutr`](https://github.com/mfansler/txcutr) by Mervin Fansler. **If
you use `txendcutr`, please cite both packages.**

Below is the citation output from using `citation('txendcutr')` in R.
Please run this yourself to check for any updates on how to cite
**txendcutr**.

``` r
print(citation('txendcutr'), bibtex = TRUE)
#> To cite package 'txendcutr' in publications use:
#> 
#>   Rocamora Pérez G, Fansler M (2026). _txendcutr: Transcriptome CUTteR
#>   with 5' and 3' End Support_. R package version 1.0.0,
#>   <https://github.com/guillermo1996/txendcutr>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {txendcutr: Transcriptome CUTteR with 5' and 3' End Support},
#>     author = {Guillermo {Rocamora Pérez} and Mervin Fansler},
#>     year = {2026},
#>     note = {R package version 1.0.0},
#>     url = {https://github.com/guillermo1996/txendcutr},
#>   }
```

Please also cite the original `txcutr` package, without which this
project would not exist:

    #> Fansler M (2025). _txcutr: Transcriptome CUTteR_. R package version
    #> 1.15.2, <https://github.com/mfansler/txcutr>.
    #> 
    #> A BibTeX entry for LaTeX users is
    #> @Manual{,
    #>   title = {txcutr: Transcriptome CUTteR},
    #>   author = {Mervin Fansler},
    #>   year = {2025},
    #>   note = {R package version 1.15.2},
    #>   url = {https://github.com/mfansler/txcutr},
    #> }

Note that `txendcutr` was only made possible thanks to Mervin Fansler’s
original work on `txcutr`, as well as many other R and bioinformatics
software authors, which are cited either in the vignettes and/or the
paper(s) describing this package.

## Code of Conduct

Please note that the `txendcutr` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  Actions](https://github.com/guillermo1996/txendcutr/actions), running
  on [Bioconductor’s Docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.23/BiocCheck)*,
  adapted from the workflow originally generated by
  *[biocthis](https://bioconductor.org/packages/3.23/biocthis)* for
  upstream `txcutr`.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.
- Unit testing is powered by
  *[testthat](https://CRAN.R-project.org/package=testthat)* (3rd
  edition).

This package originates from Mervin Fansler’s `txcutr`; see the
[upstream repository](https://github.com/mfansler/txcutr) for its own,
separately maintained development infrastructure (including a `pkgdown`
documentation site).
