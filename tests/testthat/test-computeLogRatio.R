

atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = t(mat_counts_atac_tumor),
    allele_counts = allele_counts_atac_tumor,
    features = peaks
)
rna <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = t(mat_counts_rna_tumor),
    allele_counts = allele_counts_rna_tumor,
    features = genes
)

atac_ref <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = t(mat_counts_atac_ref),
    allele_counts = allele_counts_atac_ref,
    features = peaks
)
rna_ref <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = t(mat_counts_rna_ref),
    allele_counts = allele_counts_rna_ref,
    features = genes
)

muscadet <- CreateMuscadetObject(
    omics = list(atac, rna),
    bulk.lrr = bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38"
)
muscadet_ref <- CreateMuscadetObject(
    omics = list(atac_ref, rna_ref),
    genome = "hg38"
)



test_that("computeLogRatioATAC gives a correct object as output", {

    obj_atac <- computeLogRatioATAC(
        matTumor = matCounts(muscadet)$ATAC,
        matRef = matCounts(muscadet_ref)$ATAC,
        peaksCoord = coordFeatures(muscadet)$ATAC,
        genome = slot(muscadet, "genome"),
        minReads = 1,
        minPeaks = 1,
        quiet = TRUE
    )

    expect_identical(names(obj_atac), c("matTumor", "matRef", "params", "coord"))
})

test_that("ComputeLogRatioRNA gives a correct object as output", {

    obj_rna <- computeLogRatioRNA(
        matTumor = matCounts(muscadet)$RNA,
        matRef = matCounts(muscadet_ref)$RNA,
        genesCoord = coordFeatures(muscadet)$RNA,
        genome = slot(muscadet, "genome"),
        refReads = 2,
        quiet = TRUE
    )

    expect_identical(names(obj_rna), c("matTumor","matRef","params","coord"))
})

test_that("ComputeLogRatioRNA gives a muscadet object as output", {

    obj_atac_all <- computeLogRatioATAC(
        matTumor = matCounts(muscadet)$ATAC,
        matRef = matCounts(muscadet_ref)$ATAC,
        peaksCoord = coordFeatures(muscadet)$ATAC,
        genome = slot(muscadet, "genome"),
        minReads = 1,
        minPeaks = 10,
        all_steps = TRUE,
        quiet = TRUE
    )
    expect_identical(names(obj_atac_all), c("step01","step02","step03","step04","step05","step07","step08","params","coord"))
    expect_identical(names(obj_atac_all$step08), c("matTumor","matRef","name"))
})

test_that("ComputeLogRatioRNA gives a muscadet object as output", {

    obj_rna_all <- computeLogRatioRNA(
        matTumor = matCounts(muscadet)$RNA,
        matRef = matCounts(muscadet_ref)$RNA,
        genesCoord = coordFeatures(muscadet)$RNA,
        genome = slot(muscadet, "genome"),
        refReads = 2,
        all_steps = TRUE,
        quiet = TRUE
    )
    expect_identical(names(obj_rna_all), c("step01","step02","step03","step04","step05","step06","step07","step08","params","coord"))
    expect_identical(names(obj_rna_all$step08), c("matTumor","matRef","name"))
})

test_that("ComputeLogRatio gives a correct muscadet object as output", {
    muscadet <- computeLogRatio(
        x = muscadet,
        reference = muscadet_ref,
        omic = "ATAC",
        method = "ATAC",
        minReads = 1,
        minPeaks = 1,
        quiet = TRUE
    )
    muscadet <- computeLogRatio(
        x = muscadet,
        reference = muscadet_ref,
        omic = "RNA",
        method = "RNA",
        refReads = 2,
        quiet = TRUE
    )
    expect_identical(as.character(class(muscadet)), "muscadet")
    expect_true(is.matrix(muscadet@omics[["ATAC"]]@coverage[["logratio"]][["mat"]]))
    expect_true(is.matrix(muscadet@omics[["RNA"]]@coverage[["logratio"]][["mat"]]))
})


test_that("getLogRatioBulk gives a correct output", {

    features_bulk_lrr <- getLogRatioBulk(
      x = muscadet_obj$ATAC,
      bulk.lrr = muscadet_obj$bulk.data$logratio
    )
    expect_equal(ncol(features_bulk_lrr), 4)
    expect_true(is.data.frame(features_bulk_lrr))
})






