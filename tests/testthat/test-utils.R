

test_that("addAlleleCounts() returns a muscadet object with allele table counts,
          identically as if allele counts were added through CreateMuscomicObject()" , {

    # Create muscomic objects
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
    atac_ref <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_ref,
        allele_counts = allele_counts_atac_ref,
        features = peaks
    )
    rna_ref <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_ref,
        allele_counts = allele_counts_rna_ref,
        features = genes
    )

    # Create muscadet objects
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

    # Create muscomic objects without allele data
    atac2 <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_tumor,
        features = peaks
    )
    rna2 <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_tumor,
        features = genes
    )
    atac2_ref <- CreateMuscomicObject(
        type = "ATAC",
        mat_counts = mat_counts_atac_ref,
        features = peaks
    )
    rna2_ref <- CreateMuscomicObject(
        type = "RNA",
        mat_counts = mat_counts_rna_ref,
        features = genes
    )

    # Create muscadet objects without allele data
    muscadet2 <- CreateMuscadetObject(
        omics = list(atac2, rna2),
        bulk.lrr = bulk_lrr,
        bulk.label = "WGS",
        genome = "hg38"
    )
    muscadet2_ref <- CreateMuscadetObject(
        omics = list(atac2_ref, rna2_ref),
        genome = "hg38"
    )

    # Add allele data
    obj <- addAlleleCounts(
        muscadet2,
        allele_counts = list(allele_counts_atac_tumor, allele_counts_rna_tumor))
    obj_ref <- addAlleleCounts(
        muscadet2_ref,
        allele_counts = list(allele_counts_atac_ref, allele_counts_rna_ref))

    # Check slot added
    expect_true(!is.null(obj@omics$ATAC@allelic$table.counts))
    expect_true(!is.null(obj_ref@omics$ATAC@allelic$table.counts))
    expect_true(!is.null(obj@omics$RNA@allelic$table.counts))
    expect_true(!is.null(obj_ref@omics$RNA@allelic$table.counts))

    # Compare with objects with allele data
    expect_identical(obj, muscadet)
    expect_identical(obj_ref, muscadet_ref)

})


test_that("process_allele() returns a correct data frame", {

    sc_data <- data.frame(
        ReadGroup = c("cell1", "cell2"),
        CHROM = c(1, 1),
        POS = c(1001, 1002),
        REF = c("A", "C"),
        ALT = c("G", "T"),
        SNVCount = c(3, 5),
        RefCount = c(7, 5),
        GoodReads = c(10, 10)
    )

    vcf_data <- data.frame(
        CHROM = c(1, 1),
        POS = c(1001, 1002),
        ID = c(".", "."),
        REF = c("A", "C"),
        ALT = c("G", "T"),
        QUAL = c(".", "."),
        FILTER = c(".", "."),
        INFO = c(".", "."),
        FORMAT = c("GT", "GT"),
        sample1 = c("0|1", "1|0"),
        stringsAsFactors = FALSE
    )

    expect_s3_class(process_allele(sc_data, vcf = vcf_data),
                    "data.frame")
    expect_length(process_allele(sc_data, vcf = vcf_data),
                  10)
    # expects no GT column when no vcf is provided
    expect_s3_class(process_allele(sc_data), "data.frame")
    expect_length(process_allele(sc_data), 9)
})


test_that("process_allele() with scReadcounts data input missing columns", {

    sc_data <- data.frame(
        ReadGroup = c("cell1", "cell2"),
        CHROM = c(1, 1),
        POS = c(1001, 1002),
        REF = c("A", "C"),
        ALT = c("G", "T"),
        SNVCount = c(3, 5)
    )

    vcf_data <- data.frame(
        CHROM = c(1, 1),
        POS = c(1001, 1002),
        ID = c(".", "."),
        REF = c("A", "C"),
        ALT = c("G", "T"),
        QUAL = c(".", "."),
        FILTER = c(".", "."),
        INFO = c(".", "."),
        FORMAT = c("GT", "GT"),
        sample1 = c("0|1", "1|0"),
        stringsAsFactors = FALSE
    )

    expect_error(process_allele(sc_data))
    expect_error(process_allele(sc_data, vcf = vcf_data))
})


test_that("save_as_vcf() creates a vcf file", {
    skip_on_cran()

    vcf_data <- data.frame(
        CHROM = c("chr1", "chr2"),
        POS = c(12345, 67890),
        ID = c(".", "."),
        REF = c("A", "T"),
        ALT = c("G", "C"),
        QUAL = c(".", "."),
        FILTER = c("PASS", "PASS"),
        INFO = c("AF=0.5", "AF=0.3"),
        FORMAT = c("GT", "GT"),
        SAMPLE = c("0/1", "1/1")
    )

    output_filename <- tempfile(fileext = ".vcf")
    # temp file in session, so automatically cleaned up after

    save_as_vcf(vcf_data, output_filename)
    expect_true(file.exists(output_filename))

    # Clean up
    unlink(output_filename)
})

