#' @rdname truncateTxome
#' @param txdb an object representing a transcriptome
#' @param maxTxLength the maximum length of resulting transcripts
#' @param txEnd transcript truncation end
#' @param overlapFile (optional) path to export a CSV file containing transcript
#'   overlaps (query_transcript, subject_transcript) post-truncation. If NULL,
#'   no file is exported.
#' @param ... additional arguments
#'
#' @return a \code{TxDb} object
#' @export
#'
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setGeneric("truncateTxome",
  signature = c("txdb", "maxTxLength", "txEnd", "overlapFile"),
  function(txdb, maxTxLength = 500, txEnd = "3prime", overlapFile = NULL, BPPARAM = bpparam(), ...) standardGeneric("truncateTxome")
)

#' Truncate Transcriptome
#'
#' @rdname truncateTxome
#'
#' @param txdb a \code{TxDb} object
#' @param maxTxLength the maximum length of transcripts
#' @param txEnd transcript truncation end
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   and how the method should be parallelized.
#' @param overlapFile (optional) path to export a CSV file containing transcript
#'   overlaps (query_transcript, subject_transcript) post-truncation. If NULL,
#'   no file is exported.
#' @return a \code{TxDb} object
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#'
#' ## load annotation
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#'
#' ## restrict to 'chrI' transcripts
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncateTxome(txdb)
#' txdb_w500
#'
#' ## last 100 nts per tx
#' txdb_w100 <- truncateTxome(txdb, maxTxLength = 100)
#' txdb_w100
#'
#' ## first 500 nts per tx (5' truncation)
#' txdb_5p_w500 <- truncateTxome(txdb, txEnd = "5prime")
#' txdb_5p_w500
#'
#' @importFrom GenomicRanges GRangesList GRanges mcols
#' @importFrom GenomicFeatures exonsBy
#' @importFrom txdbmaker makeTxDbFromGRanges
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom AnnotationDbi select taxonomyId
#' @importFrom S4Vectors queryHits subjectHits split
#' @importFrom methods setMethod
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @export
setMethod("truncateTxome", "TxDb", function(txdb,
                                            maxTxLength = 500,
                                            txEnd = "3prime",
                                            overlapFile = NULL,
                                            BPPARAM = bpparam()) {
  ############################################################################
  # Ensure correct values of `txEnd`
  valid_3prime <- c("3", "3'", "3p", "3prime", "3_prime")
  valid_5prime <- c("5", "5'", "5p", "5prime", "5_prime")

  if (txEnd %in% valid_3prime) txEnd <- "3prime"
  if (txEnd %in% valid_5prime) txEnd <- "5prime"

  if (!txEnd %in% c("3prime", "5prime")) stop("txEnd parameter not valid - only '3prime' or '5prime' parameters are accepted.")

  ############################################################################
  # Split exons by transcripts and create a mapping dictionary from
  # transcript_id to gene_id
  grlExons <- exonsBy(txdb, use.names = TRUE)
  dfTxGene <- suppressMessages(select(txdb, keys = names(grlExons), keytype = "TXNAME", columns = "GENEID"))
  mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)

  ############################################################################
  # Transcript truncation
  message("Truncating transcripts...")
  clipped = .clipTranscript(grlExons, maxTxLength = maxTxLength, txEnd = txEnd, BPPARAM = BPPARAM)
  message("Done.")
  
  ############################################################################
  # Removal of duplicated transcripts
  
  ## Split  exons by transcript to facilitate finding overlaps
  grlC_clipped <- S4Vectors::split(clipped, mcols(clipped)["transcript_id"])
  grlC_clipped_names <- names(grlC_clipped)
  
  ## Generate the overlap: look for exact matches between transcripts
  message("Checking for duplicate transcripts...")
  overlaps <- findOverlaps(grlC_clipped, minoverlap=maxTxLength, 
                           ignore.strand=FALSE,
                           drop.self=TRUE, drop.redundant=TRUE)
  
  ## The previous overlap approach to find duplicates will fail if the
  ## transcript is less than `maxTxLength`. The following code will fix that
  ## issue.
  
  # overlaps <- findOverlaps(grlC_clipped, type="equal", 
  #                          ignore.strand=FALSE,
  #                          drop.self=TRUE, drop.redundant=TRUE)
  
  ##  Ensure that overlaps are from the same gene
  matched_overlaps <- tibble::as_tibble(overlaps) %>%
    dplyr::mutate(queryTx = grlC_clipped_names[queryHits],
                  subjectTx = grlC_clipped_names[subjectHits]) %>% 
    dplyr::mutate(queryGene = mapTxToGene[queryTx],
                  subjectGene = mapTxToGene[subjectTx]) %>% 
    dplyr::filter(queryGene == subjectGene)
  
  ## export overlap data.frame
  if (!is.null(overlapFile) && overlapFile != "") {
    ### create parent directory
    output_dir <- dirname(overlapFile)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ### write in disk
    write.table(matched_overlaps, overlapFile, sep="\t", row.names=FALSE, quote=FALSE)
    message("Post-truncation transcript overlaps exported to: ", overlapFile)
  }
  
  ## Remove overlaps
  duplicates <- unique(matched_overlaps$queryHits)
  if (length(duplicates) > 0) {
    grlC_clipped <- grlC_clipped[-duplicates]
  }
  message(sprintf("Removed %d duplicates.", length(duplicates)))

  ############################################################################
  # Create the final exon ranges
  message("Creating exon ranges...")

  ## flatten with tx_id in metadata
  grExons <- slot(grlC_clipped, "unlistData")
  mcols(grExons)["type"] <- "exon"

  ## add gene id
  mcols(grExons)["gene_id"] <- mapTxToGene[mcols(grExons)$transcript_id]

  ## reindex exon info
  grExons <- sort(grExons)
  mcols(grExons)["exon_id"] <- seq_along(grExons)
  mcols(grExons)["exon_name"] <- NULL
  ## TODO: include `exon_rank`
  message("Done.")

  ############################################################################
  # Create the final transcript ranges
  message("Creating transcript ranges...")
  grTxs <- unlist(range(grlC_clipped))
  mcols(grTxs)["transcript_id"] <- factor(names(grlC_clipped), levels = levels(grExons$transcript_id))
  mcols(grTxs)["type"] <- "transcript"
  
  ## add gene id
  mcols(grTxs)["gene_id"] <- mapTxToGene[grTxs$transcript_id]
  
  grTxs <- grTxs[order(grTxs$transcript_id)]
  message("Done.")

  ############################################################################
  # Create the final gene ranges
  message("Creating gene ranges...")
  grGenes <- unlist(range(S4Vectors::split(grTxs, mcols(grTxs)$gene_id)))
  mcols(grGenes)["gene_id"] <- names(grGenes)
  mcols(grGenes)["type"] <- "gene"
  message("Done.")
  
  ############################################################################
  # Generate the final TxDb object
  dfMetadata <- data.frame(
    name=c("Truncated by", "Maximum Transcript Length"),
    value=c("txendcutr", maxTxLength)
  )
  
  .suppressTxDbGenomeWarning(
    makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                        taxonomyId = taxonomyId(txdb),
                        metadata = dfMetadata)
  )
})

#' Clip Transcript to Given Length
#'
#' Internal function for operating on \code{CompressedGRangesList}, where each
#' element represents a transcript and its exons.
#'
#' @param grl a \code{CompressedGRangesList} object
#' @param maxTxLength a positive integer
#' @param txEnd transcript truncation end
#'
#' @return the truncated \code{GRanges} object
#'
#' @importFrom GenomicRanges invertStrand makeGRangesFromDataFrame seqinfo
#' @importFrom methods slot
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% distinct mutate ntile row_number left_join group_by
#'   group_split bind_rows arrange filter if_else
#' 
.clipTranscript <- function(grl, maxTxLength, txEnd, BPPARAM) {
  # Merge all exons into a single GRanges and add the transcript ID
  exons_gr <- slot(.mutateEach(grl, transcript_id = names(grl)), "unlistData")
  exons_gr$transcript_id <- factor(exons_gr$transcript_id, levels = names(grl))
  
  # 5' truncation is the same a 3' truncation if we invert the strand
  if(txEnd == "5prime") exons_gr <- invertStrand(exons_gr)
  
  # Convert GRanges into tibble and remove transcript with inconsistent strands
  exons_df <- as_tibble(exons_gr[, "transcript_id"])
  exons_df <- .pruneInconsistentStrand(exons_df)
  
  # Split tibble into batches based on the number of threads registered in
  # BiocParallel and the strand
  tx_to_workers <- exons_df %>% 
    dplyr::distinct(transcript_id) %>% 
    dplyr::mutate(worker_id = dplyr::ntile(dplyr::row_number(), BPPARAM$workers))
  
  exons_by_strand_workers <- exons_df %>% 
    dplyr::left_join(tx_to_workers, by = "transcript_id", relationship = "many-to-one") %>% 
    dplyr::group_by(worker_id, strand) %>% 
    dplyr::group_split()
  
  # Apply the truncation pipeline to each group individually
  exons_truncated_df <- bplapply(exons_by_strand_workers, function(exons_split){
    split_strand = unique(exons_split$strand)
    
    if(split_strand == "+"){
      # Truncation pipeline consist in extracting exons with a cumulative width
      # less than the desired `maxTxLength` and then prune the following exon so
      # that the total width per transcript is exactly `maxTxLength`.
      exons_truncated_split <- exons_split %>% 
        dplyr::group_by(transcript_id) %>% 
        dplyr::arrange(-end, .by_group = TRUE) %>% 
        dplyr::mutate(cumLength = cumsum(width)) %>% 
        # Filter the first N + 1 exons, where N is the number of exons with cumulative width less than `maxTxLength`
        dplyr::filter(dplyr::row_number() <= sum(cumLength < maxTxLength) + 1) %>%  
        # Modify only the exon with a cumulative width higher than `maxTxLength`
        dplyr::mutate(start = dplyr::if_else(cumLength > maxTxLength, start + (cumLength - maxTxLength), start)) %>%  
        dplyr::arrange(start, .by_group = TRUE) %>% 
        dplyr::select(-cumLength, -worker_id) %>% 
        dplyr::ungroup() 
    }else if(split_strand == "-"){
      # If exons are reversed stranded, we simply need to apply the truncation
      # to the other transcript end (same logic)
      exons_truncated_split <- exons_split %>% 
        dplyr::group_by(transcript_id) %>% 
        dplyr::arrange(start, .by_group = TRUE) %>% 
        dplyr::mutate(cumLength = cumsum(width)) %>% 
        # Filter the first N + 1 exons, where N is the number of exons with cumulative width less than `maxTxLength`
        dplyr::filter(dplyr::row_number() <= sum(cumLength < maxTxLength) + 1) %>%  
        # Modify only the exon with a cumulative width higher than `maxTxLength`
        dplyr::mutate(end = dplyr::if_else(cumLength > maxTxLength, end - (cumLength - maxTxLength), end)) %>%  
        dplyr::arrange(start, .by_group = TRUE) %>% 
        dplyr::select(-cumLength, -worker_id) %>% 
        dplyr::ungroup() 
    }
    
    return(exons_truncated_split)
  }, BPPARAM = BPPARAM) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(transcript_id)
  
  exons_truncated_gr <- makeGRangesFromDataFrame(exons_truncated_df, keep.extra.columns = T, seqinfo = seqinfo(exons_gr))
  exons_truncated_gr <- sort(exons_truncated_gr)
  
  if(txEnd == "5prime") exons_truncated_gr <- invertStrand(exons_truncated_gr)
  
  return(exons_truncated_gr)
}

#' Prune GRanges of transcripts with inconsistent strand
#'
#' Internal function operating on a \code{tibble}, with rows representing exons.
#'
#' @param exons_df a \code{tibble} of exons with \code{strand} and
#'   \code{transcript_id} fields.
#'
#' @return \code{tibble} with inconsistent transcripts removed
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% group_by summarise n_distinct filter pull
#' 
.pruneInconsistentStrand <- function(exons_df){
  # GR: Group by transcripts and filter based on the number of distinct strands.
  multistrand_tx <- exons_df %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::summarise(n = dplyr::n_distinct(strand)) %>% 
    dplyr::filter(n > 1) %>% 
    dplyr::pull(transcript_id)
  
  if(length(multistrand_tx) > 0){
    warning("Some transcripts have inconsistend strand annotation! These will be ignored")
    remove_tx <- unique(multistrand_tx)
    exons_df <- exons_df %>% dplyr::filter(!transcript_id %in% remove_tx)
  }
  
  return(exons_df)
}

#' #' Invert start and end positions
#' #'
#' #' Internal function operating on a \code{tibble}, with rows representing exons.
#' #' The function inverts the starting and ending positions of exons. The pipeline
#' #' to truncate the exons can be the same to both strands if we make the
#' #' following modifications to the reverse strand:
#' #'
#' #' 1. Swap start and ending positions. 
#' #' 2. Set the genomic coordinates as negatives. 
#' #' 3. Apply the truncation process. 
#' #' 4. Revert the sign of the genomic coordinates and the start/end positions.
#' #'
#' #' @param exons_df a \code{tibble} of exons with \code{strand} and
#' #'   \code{transcript_id} fields.
#' #'
#' #' @return \code{tibble} with inverted starting and ending positions.
#' #' 
#' .strandTransform <- function(exons_df){
#'   inv_exons_df <- exons_df %>% 
#'     dplyr::mutate(start = -start, end = -end) %>% 
#'     dplyr::rename(start = end, end = start)
#'   
#'   return(inv_exons_df)
#' }