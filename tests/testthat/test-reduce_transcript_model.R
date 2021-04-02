################################################################################
library(GenomicRanges)
library(DENR)
gr_ds <-
    GRanges(
        seqnames = Rle("1", 8),
        ranges = IRanges(start = c(100, 100, 350, 350,
                                   1300, 1300, 1300, 1600),
                         end = c(1300, 1300, 1300, 1100,
                                 2500, 2500, 2250, 2250)),
        strand = Rle(strand(c("+", "-")), c(4, 4)),
        tx_name = c("t1.1", "t1.2", "t1.3", "t1.4",
                    "t2.1", "t2.2", "t2.3", "t2.4"),
        gene_id = as.character(Rle(c("g1", "g2"), c(4, 4))))

bsize <- 250
tq <- transcript_quantifier(gr_ds, bin_size = bsize,
                            transcript_name_column = "tx_name",
                            gene_name_column = "gene_id",
                            mask_start_bins = c(0, 0),
                            mask_end_bins = c(0, 0))

# model 1 on positive strand
tx_model_1 <- matrix(1, nrow = 5, ncol = 4,
                     dimnames = list(as.character(1:5), paste0("t1.", 1:4)))
tx_model_1[1, 3:4] <- 0
tx_model_1[5, 4] <- 0

# model 2 on negative strand
tx_model_2 <- matrix(1, nrow = 5, ncol = 4,
                     dimnames = list(as.character(1:5), paste0("t2.", 1:4)))
tx_model_2[5, 3:4] <- 0
tx_model_2[1, 4] <- 0

# model 3 with decimals and all 0s, granges only provides strand info and
# ranges are not made to be matched with model_3
gr_model_3 <-
    GRanges(
        seqnames = Rle("1", 4),
        ranges = IRanges(start = c(100, 100, 350, 350),
                         end = c(1300, 1300, 1300, 1100)),
        strand = Rle(strand(c("+")), c(4)),
        tx_name = c("t3.1", "t3.2", "t3.3", "t3.4"),
        gene_id = as.character(Rle(c("g3"), c(4))))

tx_model_3 <- matrix(1, nrow = 5, ncol = 4,
                     dimnames = list(as.character(1:5), paste0("t3.", 1:4)))

tx_model_3[c(1, 5), 2] <- 0.3
tx_model_3[c(1, 5), 3] <- 0.7
tx_model_3[, 4] <- 0

# test inputs
test_that("Inputs are correct for reduce_transcript_models", {
    # not a list
    expect_error(reduce_transcript_models(tx_model_1, gr_ds,
                                          "tx_name"), "invalid input")
    # element in list is not matrix
    expect_error(reduce_transcript_models(list(tx_model_1, NA), gr_ds,
                                          "tx_name"), "invalid input")
    expect_error(reduce_transcript_models(list(tx_model_1, 3), gr_ds,
                                          "tx_name"), "invalid input")
    # invalid argument
    expect_error(reduce_transcript_models(list(tx_model_1), gr_ds,
                                          "tx_name", "Round"),
                 "invalid arguement")
})

test_that("Transcript models are reduced correctly", {
    true_reduced_modle_1 <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t1.1 = rep(1, 5),
            t1.3 = c(0, rep(1, 4)),
            t1.4 = c(0, rep(1, 3), 0)
        ))

    true_reduced_modle_2 <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t2.1 = rep(1, 5),
            t2.3 = c(rep(1, 4), 0),
            t2.4 = c(0, rep(1, 3), 0)
        ))

    true_reduced_modle_3rf <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t3.1 = rep(1, 5),
            t3.2 = c(0, rep(1, 3), 0),
            t3.4 = rep(0, 5)
        ))

    true_reduced_modle_3c <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t3.1 = rep(1, 5),
            t3.4 = rep(0, 5)
        ))

    expect_equal(
        reduce_transcript_models(list(tx_model_1, tx_model_2),
                                 gr_ds, "tx_name")[[1]],
        list(true_reduced_modle_1, true_reduced_modle_2)
    )

    expect_equal(
        reduce_transcript_models(list(tx_model_3), gr_model_3,
                                 "tx_name", "round")[[1]],
        list(true_reduced_modle_3rf)
    )

    expect_equal(
        reduce_transcript_models(list(tx_model_3), gr_model_3,
                                 "tx_name", "floor")[[1]],
        list(true_reduced_modle_3rf)
    )

    expect_equal(
        reduce_transcript_models(list(tx_model_3), gr_model_3,
                                 "tx_name", "ceiling")[[1]],
        list(true_reduced_modle_3c)
    )
})

test_that("Transcript groupings are correct", {
    true_group_1 <-
        data.frame(
            tx_name = c("t1.1", "t1.2", "t1.3", "t1.4",
                        "t2.1", "t2.2", "t2.3", "t2.4"),
            group = c(rep(1, 4), rep(2, 4)),
            model = c(1, 1, 2, 3, 1, 1, 2, 3),
            tss_set = c(1, 1, 2, 2, 1, 1, 2, 2),
            tts_set = c(1, 1, 1, 2, 1, 1, 1, 2)
        )

    expect_equal(reduce_transcript_models(list(tx_model_1, tx_model_2),
                                          gr_ds, "tx_name")[[2]],
                 true_group_1)

})
