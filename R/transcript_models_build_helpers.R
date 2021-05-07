#' Create bins for multiple groups of overlapped transcripts
#'
#' Function for generating bins for mutiple groups of overlapped transcripts.
#'
#' @param transcript_groups A \code{\link[GenomicRanges]{GRangesList-class}}
#' object,each element in the list contains transcripts from one gene or several
#' genes which are overlapped.
#' @param bin_size An integer, used to tile the gene region. Default is 250bp.
#' @return A \code{\link[GenomicRanges]{GRangesList-class}} object, containing
#' the binned region.
#' @export
create_bins <- function(transcript_groups, bin_size = 250) {
  # check input class
  if (!methods::is(transcript_groups, "GRangesList")) {
    stop("transcript_groups is not a GRangesList object")
  }
  if (!methods::is(bin_size, "numeric") || bin_size < 0) {
    stop("bin_size is not a positive number")
  }

  # First reduce transcript groups to single range using large gapwidth to
  # ensure single range
  red_tx_grps <- unlist(GenomicRanges::reduce(transcript_groups, min.gapwidth = 1e9))
  # Pre-compute number of bins in each range
  bin_count <- ceiling(GenomicRanges::width(red_tx_grps) / bin_size)
  # Pre-build rle for each element
  chrom <-
    S4Vectors::Rle(as.vector(GenomicRanges::seqnames(red_tx_grps)), bin_count)
  strnd <-
    S4Vectors::Rle(as.vector(GenomicRanges::strand(red_tx_grps)), bin_count)
  # Create Granges as single GRanges
  tiles <- GenomicRanges::GRanges(
    chrom,
    IRanges::IRanges(start = unlist(
      seqv(from = GenomicRanges::start(red_tx_grps),
           to = GenomicRanges::end(red_tx_grps),
           by = bin_size)),
      width = bin_size),
    strand = strnd)
  # Return split GRangesList
  return(GenomicRanges::split(x = tiles, rep(names(transcript_groups), bin_count))
  )
}

#' Vectorized seq.default
#'
#' Version of the \code{seq} function that takes vectorized arguments
#'
#' @inheritParams base::seq
seqv <- Vectorize(seq.default, vectorize.args = c("from", "to"))

#' Create model masks
#'
#' Function for generating masks
#'
#' @param transcript_models A list of matricies containing transcript models.
#' @param mask_start_bins A numeric vertor with length two which giving the
#' number of bins should be masked near the start of a transcript. The first
#' integer is the number of bins will be masked inside the transcript, while
#' the second interger is the number outside the transcript. Default c(0, 0).
#' @param mask_end_bins A numeric vertor with length two which giving the
#' number of bins should be masked near the end of a transcript. The first
#' integer is the number of bins will be masked inside the transcript, while
#' the second interger is the number outside the transcript. Default c(0, 0).
#' @param strand a vector containing one entry per group of transcript models
#' indicating the strand.
#'
#' @rdname create_model_masks
#'
#' @include type_checkers.R
#' @return A list of vectors, where each vector indicates the bins that should
#' be masked
#' @export
#'
#'
create_model_masks <- function(transcript_models, strand,
                               mask_start_bins = NULL,
                               mask_end_bins = NULL) {
  if (is.null(mask_start_bins)) mask_start_bins <- c(0, 0)
  if (is.null(mask_end_bins)) mask_end_bins <- c(0, 0)
  # ***Checks***
  if (!is_strand_vector(strand, allow_star = FALSE)) {
    stop("strand must be a vector containing only '+','-'")
  }
  if (!is_matrix_list(transcript_models)) {
    stop("transcript models must be a list of matricies")
  }
  if (!(is_length_two_vector(mask_start_bins) &&
        is_length_two_vector(mask_end_bins))) {
    stop("the length of input vector should be two,
         and the numbers in the vertor should be non-negative integers")
  }
  if (length(strand) != length(transcript_models)) {
    stop("strand and transcript_models must be the same length")
  }
  # ***End checks***
  mask_start_in <- mask_start_bins[1]
  mask_start_out <- mask_start_bins[2]
  mask_end_in <- mask_end_bins[1]
  mask_end_out <- mask_end_bins[2]

  # Get indicies of masked bins per transcript group
  all_masks <- mapply(function(tx_models, strand) {
    group_masks <- apply(tx_models, 2, function(tx) {
      first_non_zero <- min(which(tx != 0))
      last_non_zero <- max(which(tx != 0))
      # Orient the first/last axis along the 5' to 3' axis
      if (strand == "+") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_start_in,
                             by = 1),
                         seq(from = first_non_zero - 1,
                             length.out = mask_start_out,
                             by = -1),
                         seq(from = last_non_zero,
                             length.out = mask_end_in,
                             by = -1),
                         seq(from = last_non_zero + 1,
                             length.out = mask_end_out,
                             by = 1)))
      } else if (strand == "-") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_end_in,
                             by = 1),
                         seq(from = first_non_zero - 1,
                             length.out = mask_end_out,
                             by = -1),
                         seq(from = last_non_zero,
                             length.out = mask_start_in,
                             by = -1),
                         seq(from = last_non_zero + 1,
                             length.out = mask_start_out,
                             by = 1)))
      }
      # Ensure all masked bins are within the tx range being considered
      mask <- mask[mask >= 1 & mask <= length(tx)]
      return(mask)
    })
    # Ensure masks are returned as one 1D vector per transcript group
    return(sort(unique(as.vector(unlist(group_masks)))))
  }, tx_models = transcript_models, strand = as.list(strand), SIMPLIFY = F)

  return(all_masks)
  }

#' Function for generating masks for additional genomic regions
#'
#' @param add_mask A \code{\link[GenomicRanges]{GRanges-class}}
#' object records genomic regions being masked. If seqnames in add_mask are not
#' available in transcripts, entries associated with these seqnames will be removed.
#' @rdname create_additional_masks
#' @inheritParams create_transcript_models
#' @return A list of vectors, where each vector indicates the bins that should
#' be masked
#' @export
#'
create_additional_masks <- function(bins, add_mask) {
    if (!is.null(add_mask)) {
        # Get non-redundant regions
        add_mask <- GenomicRanges::reduce(add_mask)
        # Filter seqnames which are in mask regions but not in bins
        seq_filter <- as.vector(GenomeInfoDb::seqnames(add_mask)) %in%
            GenomeInfoDb::seqlevels(bins)
        add_mask <- add_mask[seq_filter, ]
        # Find overlaps between bins and masks
        if (length(add_mask) > 0) {
            add_masks <-
                lapply(bins, function(x) {
                    S4Vectors::queryHits(
                    GenomicRanges::findOverlaps(x, add_mask, select = "all"))})
            } else {
                add_masks <- list(NULL)
        }
    } else {
        add_masks <- list(NULL)
    }
    return(add_masks)
}

#' Function for combining masks
#'
#' @param model_masks,add_masks A list of vectors, where each vector indicates
#' the bins that should be masked.
#' @rdname combine_masks
#' @return A list of vectors, where each vector has combined the masks, and indicates
#' the bins that should be masked.
#' @export
#'
combine_masks <- function(model_masks, add_masks) {
    # if both masks exist, combine them
    if (any(lengths(model_masks) > 0) && any(lengths(add_masks) > 0)) {
        all_masks <- mapply(function(x, y) {
            unique(sort(c(x, y)))
        }, model_masks, add_masks[names(model_masks)], SIMPLIFY = FALSE)
        return(all_masks)
    }
    # if only additional masks, returns it
    if (all(lengths(model_masks) == 0) && any(lengths(add_masks) > 0)) return(add_masks)
    # else return model masks
    return(model_masks)
}

#' Function for scaling bins for masking additional genomic regions
#'
#' @param add_mask A \code{\link[GenomicRanges]{GRanges-class}}
#' object records genomic regions being masked. If seqnames in add_mask are not
#' available in transcripts, entries associated with these seqnames will be removed.
#' @rdname scale_additional_masks
#' @inheritParams create_transcript_models
#' @return A list of vectors, where each vector indicates the fractions to scale
#' the bins
#' @export
#'
scale_additional_masks <- function(bins, add_mask, bin_size) {
    # Get non-redundant masking regions
    add_mask <- GenomicRanges::reduce(add_mask)
    # Keep same seqlevels in both bins and masks for IRanges::Views
    if (!identical(sort(GenomeInfoDb::seqlevels(add_mask)),
                   sort(GenomeInfoDb::seqlevels(bins)))) {
        seq_filter <-
            intersect(GenomeInfoDb::seqlevels(add_mask),
                      GenomeInfoDb::seqlevels(bins))
        add_mask <-
            add_mask[as.vector(GenomicRanges::seqnames(add_mask)) %in% seq_filter]
        GenomeInfoDb::seqlevels(add_mask) <- GenomeInfoDb::seqlevelsInUse(add_mask)
        unl_bins <- unlist(bins)
        unl_bins <-
            unl_bins[as.vector(GenomicRanges::seqnames(unl_bins)) %in% seq_filter]
        GenomeInfoDb::seqlevels(unl_bins) <-
            GenomeInfoDb::seqlevelsInUse(unl_bins)
        gp <- names(unl_bins)
        names(unl_bins) <- NULL
        bins <- split(unl_bins, gp)
    }
    # return NA if no masking is left
    if (length(add_mask) == 0) return(NA)
    # get coverage for masked regions
    cov_mask <- GenomicRanges::coverage(add_mask)
    # 1 indicates maksed regions, while 0 means not masked, here all regions
    # without masks
    cov_log <- cov_mask == 0L
    add_scale <- lapply(bins, function(x) {
        # view not masked regions
        views <- IRanges::Views(cov_log, x)
        # pick up chromosome being hit
        views <- views[[which(lengths(views) > 0)]]
        return(IRanges::viewSums(views) / bin_size)
    })
    return(add_scale)
}

#' Mask transcripts
#'
#' Function for generating masks
#'
#' @param transcript_models A list of matricies containing transcript models
#' @param masks a list of
#'
#' @rdname mask_transcripts
#'
#' @include type_checkers.R
#' @return A list of matrices where each row in the matrix corresponds to a
#' transcript and each column is a bin
#' @export
mask_transcripts <- function(transcript_models, masks) {
  # ***Checks***
  if (!is_matrix_list(transcript_models)) {
    stop("transcript models must be a list of matricies")
  }
  if (length(masks) != length(transcript_models)) {
    stop("masks and transcript_models must be the same length")
  }

  # ***End checks***

  # Mask transcript models
  masked_transcripts <- mapply(function(tx_models, masks) {
    tx_models[masks, ] <- 0
    return(tx_models)
  }, tx_models = transcript_models, masks = masks, SIMPLIFY = F)
  return(masked_transcripts)
}

#' Create transcript models
#'
#' Function for generating binned transcript models given a set of bins
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param bins A \code{\link[GenomicRanges]{GRangesList-class}} object,
#' each element in the list contains the bins for each corresponding group of
#' transcripts
#' @param bin_size width of bins
#' @param transcript_name_column A string that indicates which column in the
#' GRanges object contain the transcript names
#' @return A list of matrices where each row in the matrix corresponds to a
#' bin and each column is a transcript
#' @export
create_transcript_models <- function(transcripts, bins, bin_size,
                                     transcript_name_column) {
  # check input class
  if (!methods::is(transcripts, "GRanges")) {
    stop("transcripts is not a GRanges object")
  }
  # Compute percent overlap of each transcript per bin and cast to matrix
  # First find all pairwise overlaps
  ovr <- IRanges::findOverlapPairs(transcripts, bins)
  # Compute bin percent overlap for all intersections
  ovr_intersect <- IRanges::pintersect(ovr, drop.nohit.ranges = FALSE)
  ovr_val <- GenomicRanges::width(ovr_intersect) / bin_size
  # Compute the dimensions of model matrix for each loci
  d <- S4Vectors::elementNROWS(ovr_val)
  # Pre-extract transcript names for ease of access
  transcript_names <-
    GenomicRanges::values(transcripts)[, transcript_name_column]
  # Reshape intersection output into matricies of the correct dimension
  tx_matrix_models <- lapply(S4Vectors::split(d, names(d))[names(bins)], function(x) {
      # Retrive rows in transcripts that each intersection pair came from
      lookup_ind <- which(names(ovr_val) == names(x[1]))
      # Place into matrix of correct dimension
      m <- matrix(data = BiocGenerics::unlist(ovr_val[lookup_ind]),
                  nrow = x[1],
                  ncol = length(x))
      # Set column names to correct transcript names
      colnames(m) <- transcript_names[lookup_ind]
      return(m)
    }
  )
  return(tx_matrix_models)
}

#' @title Group transcripts
#'
#' @description
#' Creates groups of transcripts on the same strand based on
#' proximity. The groups are constructed of connected transcripts
#' such that a pair of transcripts are considered connected if
#' they are within a given distance \emph{d} of each other. The
#' group then consists of all transcripts that are connected to
#' at least one other member.
#'
#' @param transcript_granges \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param distance the distance within which two transcripts are
#' considered connected
#'
#' @return A \code{\link[GenomicRanges]{GenomicRangesList-class}}
#' object
#'
#' @rdname group_transcripts
#' @export
group_transcripts <- function(transcript_granges, distance = 0) {
  # Check that txdb is an actual transcript database object
  if (!methods::is(transcript_granges, "GRanges")) {
    stop("transcript is not a TxDb object")
  }
  # Ensure everything is sorted, otherwise later assignment assumptions fail
  transcript_granges <- GenomeInfoDb::sortSeqlevels(transcript_granges)
  transcript_granges <- GenomicRanges::sort(transcript_granges)
  # Assign a unique id to each row
  transcript_granges$unique_id <- seq_len(length(transcript_granges))
  # Create expanded granges
  tx_granges_expand <- GenomicRanges::resize(
    transcript_granges,
    width = GenomicRanges::width(transcript_granges) + floor(distance / 2),
    fix = "end")
  tx_granges_expand <- GenomicRanges::resize(
    tx_granges_expand,
    width = GenomicRanges::width(tx_granges_expand) + ceiling(distance / 2),
    fix = "start")
  tx_granges_expand <- GenomicRanges::trim(tx_granges_expand)
  # Reduce granges to get unions of overlapping ranges which will become the
  # transcript groups. The -/+ 1 bit is to assure that perfectly adjoining
  # groups are not merged
  GenomicRanges::end(tx_granges_expand) <- GenomicRanges::end(tx_granges_expand) - 1
  tx_reduce <- GenomicRanges::reduce(tx_granges_expand, ignore.strand = FALSE)
  GenomicRanges::end(tx_reduce) <- GenomicRanges::end(tx_reduce) + 1
  # Get the group assignments
  group_assignment <- GenomicRanges::findOverlaps(query = transcript_granges,
                                                  subject = tx_reduce)
  # Create group vector string
  uid <- paste0(S4Vectors::subjectHits(group_assignment), "_",
                GenomicRanges::strand(transcript_granges))

  # Split the transcripting into their groups
  gr_groups <- GenomicRanges::split(transcript_granges, f = uid)
  return(gr_groups)
}

#' Reduce transcript models
#'
#' Combine identical transcript models and get non-redundent and reduced
#' transcript models
#'
#' @param transcript_models_ls A list of matrices where each row in the matrix
#' corresponds to a bin and each column is a transcript.
#' @param transcripts a \link[GenomicRanges]{GRanges-class} object that must
#' contain a metadata column with a transcript id and may contain an additional
#' column with a gene id.
#' @param transcript_name_column A string that indicates which column in the
#' GRanges object contain the transcript names
#' @param bin_operation Three different modes to deal with decimals in the
#' transript model (due to partial overlap of the first or last exon and bins).
#' Either "ceiling", "floor", or "round" (default: "round").
#' @return  A list of matrices holding the reduced transcript models and a
#' dataframe holding each transcript belongs to which group and reduced model.
#' @export
reduce_transcript_models <-
  function(transcript_models_ls,
           transcripts,
           transcript_name_column,
           bin_operation = c("round", "floor", "ceiling")
           ) {
    # check transcript class
    tx_model_class_matrix <- unique(sapply(transcript_models_ls, is.matrix))
    if (!all(tx_model_class_matrix)) {
      stop("invalid input for transcript models")
    }
    # check bin_operation
    if (length(bin_operation) > 1) bin_operation <- bin_operation[[1]]
    if (!bin_operation %in% c("ceiling", "floor", "round")) {
      stop("invalid arguement for bin_operation, it should be 'ceiling', 'floor', ",
           "or 'round'.")
    }
    # deal with decimals in transcript models (only modify the starts and ends)
    integer_models <-
        lapply(transcript_models_ls, function(txm) {
            # handle matrix with row number equal to 1 otherwise apply will return
            # vector rather than matrix and cause subseting issues later
            if (nrow(txm) == 1) {
                txm <- do.call(bin_operation, list(txm))
                return(txm)
                }
            txm <- apply(txm, 2, function(x) {
                idx <- which(x > 0)
                if (length(idx) > 0) {
                    x[idx[1]] <- do.call(bin_operation, list(x[idx[1]]))
                    x[idx[length(idx)]] <-
                        do.call(bin_operation, list(x[idx[length(idx)]]))
                }
                return(x)
            })
            return(txm)
        })

    # compute reduced transcript groups
    tx_groups <- lapply(integer_models, group_trancript_models)

    # get reduced tx models
    get_reduced_tx_models <- function(tx_models, tx_groups) {
      tx_models[, unlist(lapply(tx_groups, function(x) x[[1]])), drop = FALSE]
    }

    reduced_tx_models <- mapply(get_reduced_tx_models,
                                integer_models,
                                tx_groups, SIMPLIFY = FALSE)

    # create group and model identifier for each transcript
    group_len <- lapply(tx_groups, lengths)
    tx_num <- sapply(group_len, sum)
    tx_group_model <- data.frame(tx_name = unlist(tx_groups),
                                 group = rep(seq_along(tx_num), tx_num),
                                 model = unlist(mapply(rep, lapply(group_len, seq_along),
                                                       group_len, SIMPLIFY = FALSE)))
    rownames(tx_group_model) <- NULL
    # create groups based on TSS and TTS positions
    get_start_end_set <- function(tx_models) {
        # get first and last non-zero values in model
        tx_models <- tx_models != 0
        tx_start_end <- apply(tx_models, 2, function(x) {
            vec <- which(x)
            # deal with vector with all 0s
            if (is.integer(vec) && length(vec) == 0L) vec <- c(0, 0)
            return(c(vec[1], vec[length(vec)]))
        })
        # get unique start and end positions
        start_bin <- unique(tx_start_end[1, ])
        end_bin <- unique(tx_start_end[2, ])
        # a dataframe record the grouping for start and end
        start_end_df <-
            data.frame(tx_name = colnames(tx_start_end),
                       start_set = match(tx_start_end[1, ], start_bin),
                       end_set = match(tx_start_end[2, ], end_bin))
        return(start_end_df)
    }
    start_end_dfs <- lapply(integer_models, get_start_end_set)
    start_end_df <- data.table::rbindlist(start_end_dfs)
    tx_group_model <-
        data.table::merge.data.table(tx_group_model, start_end_df,
                                     by = "tx_name", all.x = TRUE, sort = FALSE)
    # get strand info
    tx_group_model <-
        data.table::merge.data.table(
            tx_group_model,
            data.table::data.table(
                tx_name = GenomicRanges::values(transcripts)[, transcript_name_column],
                strand = as.vector(GenomicRanges::strand(transcripts))
            ), by = "tx_name", all.x = TRUE, sort = FALSE
        )
    # get tss and tts for transcripts on different strand
    tx_group_model$tss_set <-
        ifelse(tx_group_model$strand == "+",
               tx_group_model$start_set,
               tx_group_model$end_set)
    tx_group_model$tts_set <-
        ifelse(tx_group_model$strand == "+",
               tx_group_model$end_set,
               tx_group_model$start_set)
    tx_group_model <-
        tx_group_model[c("tx_name", "group", "model", "tss_set", "tts_set")]
    return(list(reduced_tx_models, tx_group_model))
}

#' Group transcripts
#'
#' Group transripts if they have identical transcript model.
#'
#' @param tx_models A processed matrix without decimals in transcript model
#' @seealso \code{\link{reduce_transcript_models}},
#' @return  A list of grouped transcript names.

group_trancript_models <- function(tx_models) {
  pool <- colnames(tx_models)
  groups <- list()
  while (length(pool) > 1) {
    identical_models <- which(
      colSums(abs(tx_models[, pool[-1], drop = FALSE] - tx_models[, pool[1]])) == 0
    )
    groups <- append(groups, list(c(pool[1], names(identical_models))))
    pool <- pool[-c(1, identical_models + 1)]
  }
  if (length(pool) == 1) groups <- append(groups, pool)
  return(groups)
}
