

test_that("assignClusters() returns an updated muscadet object", {

    data("muscadet_obj", package = "muscadet")
    # remove cna calling result
    muscadet_obj@cnacalling <- list()

    # Select clustering result for partition 0.6
    obj1 <- assignClusters(muscadet_obj, partition = 0.6)

    expect_true(length(obj1@cnacalling) != 0)
    expect_identical(length(table(obj1@cnacalling$clusters)), as.integer(2)) # 2 clusters

    # Assign custom clusters
    cell_names <- Reduce(union, SeuratObject::Cells(muscadet_obj))
    n1 <- sample(1:length(cell_names), 1)
    n2 <- length(cell_names) - n1
    custom_clusters <- c(rep.int(1, n1), rep.int(2, n2))
    names(custom_clusters) <- cell_names

    obj2 <- assignClusters(muscadet_obj, clusters = custom_clusters, redo_imputation = FALSE)

    expect_true(length(obj2@cnacalling) != 0)
    expect_length(table(obj2@cnacalling$clusters), 2) # 2 clusters

    # Assign clusters with remapping
    # example to remap from 5 clusters to 4 by merging clusters 1 and 2
    clusters <- muscadet_obj$clustering$clusters[["1"]] # res = 1
    mapping <- c("1" = 1, "2" = 1, "3" = 2, "4" = 3, "5" = 4)

    obj3 <- assignClusters(muscadet_obj, clusters = clusters, mapping = mapping)

    expect_true(length(obj3@cnacalling) != 0)
    expect_length(table(obj3@cnacalling$clusters), 4) # 4 clusters

    obj4 <- assignClusters(muscadet_obj, partition = 1, mapping = mapping)

    expect_true(length(obj4@cnacalling) != 0)
    expect_length(table(obj4@cnacalling$clusters), 4) # 4 clusters
})


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


test_that("mergeCounts() returns an updated muscadet object with data frames in the cnacalling slot", {

    # Load example muscadet objects
    data(muscadet_obj)
    data(muscadet_obj_ref)

    # remove cna calling result
    muscadet_obj@cnacalling <- list()
    muscadet_obj <- assignClusters(muscadet_obj, partition = 0.6)

    # Merge counts from all omics from both sample and reference
    obj <- mergeCounts(muscadet_obj, muscadet_obj_ref)

    expect_true(length(obj@cnacalling) != 1)
    expect_identical(class(obj@cnacalling$allelic.counts), "data.frame")
    expect_identical(class(obj@cnacalling$coverage.counts), "data.frame")
    expect_identical(class(obj@cnacalling$combined.counts), "data.frame")

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

