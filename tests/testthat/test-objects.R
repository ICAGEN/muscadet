# CreateMuscomicObject ---------------------------------------------------------

test_that("CreateMuscomicObject returns a correct muscomic object", {
  atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat_counts_atac_tumor,
    allele_counts = allele_counts_atac_tumor,
    features = peaks
  )
  rna <- CreateMuscomicObject(
      type = "RNA",
      mat_counts = mat_counts_rna_tumor,
      allele_counts = allele_counts_rna_tumor,
      features = genes
  )
  # class
  expect_identical(as.character(class(atac)), "muscomic")
  expect_identical(as.character(class(rna)), "muscomic")
  # slots
  expect_identical(slotNames(atac), c("type","label.omic","coverage","allelic"))
  expect_identical(slotNames(rna), c("type","label.omic","coverage","allelic"))
  # type
  expect_identical(slot(atac, "type"), "ATAC")
  expect_identical(slot(rna, "type"), "RNA")
  # label
  expect_identical(slot(atac, "label.omic"), "scATAC-seq")
  expect_identical(slot(rna, "label.omic"), "scRNA-seq")
  # cov
  expect_identical(class(slot(atac, "coverage")), "list")
  expect_identical(class(slot(rna, "coverage")), "list")
  expect_identical(names(slot(atac, "coverage")), c("mat.counts","table.counts","coord.features","label.features"))
  expect_identical(names(slot(rna, "coverage")), c("mat.counts","table.counts","coord.features","label.features"))
  # dgCMatrix
  expect_identical(as.character(class(slot(atac, "coverage")[["mat.counts"]])), "dgCMatrix")
  expect_identical(as.character(class(slot(rna, "coverage")[["mat.counts"]])), "dgCMatrix")
})

test_that("CreateMuscomicObject without a matrix in input: error", {
  expect_error(
    CreateMuscomicObject(
      type = "ATAC",
      mat_counts = "mat",
      allele_counts = allele_counts_atac_tumor,
      features = peaks
    )
  )
})

test_that(
    "CreateMuscomicObject with incorrect `type` argument: error",
    {
        expect_error(
            CreateMuscomicObject(
                type = "",
                mat_counts = mat_counts_atac_tumor,
                allele_counts = allele_counts_atac_tumor,
                features = peaks
            )
        )
    }
)

# CreateMuscomicObject ---------------------------------------------------------

test_that("CreateMuscadetObject returns a correct muscadet object", {
    atac <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_tumor,
        allele_counts = allele_counts_atac_tumor,
        features = peaks
    )
    rna <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_tumor,
        allele_counts = allele_counts_rna_tumor,
        features = genes
    )
    muscadet <- CreateMuscadetObject(
        omics = list(atac, rna),
        bulk.lrr = bulk_lrr,
        bulk.label = "WGS",
        genome = "hg38"
    )
    # class
    expect_identical(as.character(class(muscadet)), "muscadet")
    # slots
    expect_identical(slotNames(muscadet), c("omics", "bulk.data", "clustering", "cnacalling", "genome"))
    # omics
    expect_identical(class(slot(muscadet, "omics")), "list")
    expect_identical(as.character(class(slot(muscadet, "omics")[[1]])), "muscomic")
})

test_that("CreateMuscadetObject with identical muscomic: error of identical labels", {
    atac <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_tumor,
        allele_counts = allele_counts_atac_tumor,
        features = peaks
    )
    rna <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_tumor,
        allele_counts = allele_counts_rna_tumor,
        features = genes
    )
    expect_error(
        muscadet <- CreateMuscadetObject(
            omics = list(atac, atac),
            bulk.lrr = bulk_lrr,
            bulk.label = "WGS",
            genome = "hg38"
        )
    )
})

test_that("CreateMuscadetObject with incorrect genome : error of incorrect genome", {
    atac <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_tumor,
        allele_counts = allele_counts_atac_tumor,
        features = peaks
    )
    rna <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_tumor,
        allele_counts = allele_counts_rna_tumor,
        features = genes
    )
    expect_error(
        muscadet <- CreateMuscadetObject(
            omics = list(atac, rna),
            bulk.lrr = bulk_lrr,
            bulk.label = "WGS",
            genome = "test"
        )
    )
})


# Methods ----------------------------------------------------------------------


test_that("Check access methods", {

    data("muscadet_obj")

    obj <- muscadet_obj
    obj["clustering"] <- muscadet_obj["clustering"]

    expect_length(obj["clustering"], 9)
    expect_error(muscadet_obj["abc"])

    expect_length(obj["ATAC"], 1)
    expect_length(obj["ATAC"]["type"], 1)
    expect_error(obj["ATAC"]["abc"])

    expect_length(obj$ATAC, 1)
    expect_length(obj$ATAC$type, 1)
    expect_error(obj$ATAC$abc)


})


test_that("Check show methods", {

    data("muscadet_obj")
    expect_output(show(muscadet_obj), regexp = "muscadet")

    atac <- muscadet_obj$ATAC
    expect_output(show(atac), regexp = "muscomic")
})


test_that("Check Cells/Features methods outputs", {
    atac <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_tumor,
        allele_counts = allele_counts_atac_tumor,
        features = peaks
    )
    rna <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_tumor,
        allele_counts = allele_counts_rna_tumor,
        features = genes
    )
    muscadet <- CreateMuscadetObject(
        omics = list(atac, rna),
        bulk.lrr = bulk_lrr,
        bulk.label = "WGS",
        genome = "hg38"
    )

    expect_identical(Cells(atac), colnames(slot(atac, "coverage")[["mat.counts"]]))
    expect_identical(Cells(rna), colnames(slot(rna, "coverage")[["mat.counts"]]))
    expect_identical(Cells(muscadet), lapply(slot(muscadet, "omics"), function(omic) {
        return(colnames(slot(omic, "coverage")[["mat.counts"]]))
    }))

    expect_identical(Features(atac), rownames(slot(atac, "coverage")[["mat.counts"]]))
    expect_identical(Features(rna), rownames(slot(rna, "coverage")[["mat.counts"]]))
    expect_identical(Features(muscadet), lapply(slot(muscadet, "omics"), function(omic) {
        return(rownames(slot(omic, "coverage")[["mat.counts"]]))
    }))

    muscadet_monoomic <- CreateMuscadetObject(
        omics = list(atac),
        bulk.lrr = bulk_lrr,
        bulk.label = "WGS",
        genome = "hg38"
    )
    expect_identical(Cells(muscadet_monoomic), lapply(slot(muscadet_monoomic, "omics"), function(omic) {
        return(colnames(slot(omic, "coverage")[["mat.counts"]]))
    }))
})

