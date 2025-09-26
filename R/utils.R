#' Add allele counts to a `muscadet` object
#'
#' This function adds allele counts data to the `omics` of a
#' [`muscadet`][muscadet-class] object. The data frames in the `allele_counts` list
#' are assigned to the `allelic` slots of `omics`.
#'
#' @param x A [`muscadet`][muscadet-class] object.
#' @param allele_counts A list of data frames where each data frame contains
#'   allele counts for a specific omic (`list`). The list must have the same
#'   length and order as the number of omics in the `x` object. Each data frames
#'   must contain the following columns : `cell`, `id`, `CHROM`, `POS`, `REF`,
#'   `ALT`, `RD`, `AD`, `DP`, `GT`. See [allele_counts] for details.
#'
#' @return
#' A modified [`muscadet`][muscadet-class] object with updated allele counts in the
#' `allelic` slot of each `muscomic` object in the `omics` slot.
#'
#' @note
#' As the allele counts are not used for computing log R
#' ratios and cell clustering, they are not mandatory at the creation of
#' `muscomic` and `muscadet` objects. The allele counts data can be added to
#' objects later with this `addAlleleCounts` function, before using the
#' [mergeCounts()] function.
#'
#' This function is also useful to add allele counts for individual-specific
#' variant positions to a common `muscadet` object, for example for the
#' reference `muscadet` object: a common `muscadet` object with reference
#' coverage data can be stored as a unique object to use as reference for
#' computing log R ratios for all samples, and then it can updated with allele
#' counts at individual-specific variant positions (e.g. found by bulk
#' sequencing) before Copy Number Alterations (CNAs) calling.
#'
#' @seealso [CreateMuscomicObject()], [mergeCounts()]
#'
#' @importFrom stringr str_remove
#' @importFrom gtools mixedsort
#'
#' @export
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#' # data("muscadet_obj_ref")
#'
#' # Add allele counts data frames to muscadet objects
#' muscadet_obj <- addAlleleCounts(
#'     muscadet_obj,
#'     allele_counts = list(allele_counts_atac_tumor, allele_counts_rna_tumor))
#' muscadet_obj_ref <- addAlleleCounts(
#'     muscadet_obj_ref,
#'     allele_counts = list(allele_counts_atac_ref, allele_counts_rna_ref))
#'
addAlleleCounts <- function(x, allele_counts) {

    # Check that x is a muscadet object
    stopifnot("The argument 'x' must be a muscadet object." = inherits(x, "muscadet"))

    # Check that allele_counts is a list of dataframes
    stopifnot("The argument 'allele_counts' must be a list." = inherits(allele_counts, "list"))
    stopifnot("The argument 'allele_counts' must be a list of data frames." = all(unlist(
        lapply(allele_counts, function(x)
            inherits(x, "data.frame"))
    )))

    # Ensure allele_counts list matches the number of omics in x
    if (length(allele_counts) != length(x@omics)) {
        stop("The 'allele_counts' list must have the same length as the number of omics in 'x'.")
    }

    # Assign names to allele_counts based on omics names in x
    names(allele_counts) <- names(x@omics)
    # Extract types from omics in x
    type.omic <- unlist(lapply(x@omics, function(x) x$type))

    # Format allele count data frames
    allele_counts <- lapply(allele_counts, function(allele_df) {
        # Remove "chr" if necessary
        allele_df$CHROM <- stringr::str_remove(allele_df$CHROM, "chr")
        # Ordered chromosomes
        allele_df$CHROM <- ordered(allele_df$CHROM, levels = gtools::mixedsort(unique(allele_df$CHROM)))
        # Reorder data
        allele_df <- allele_df[order(allele_df$CHROM, allele_df$POS), ]
        return(allele_df)

    })

    # Add allele counts to each omic in the muscadet object
    for (i in names(x@omics)) {
        slot(x@omics[[i]], "allelic") <- list(table.counts = data.frame(omic = type.omic[[i]], allele_counts[[i]]))
    }

    return(x)
}


#' Save a data frame as a VCF File
#'
#' This function saves a given data frame in [variant calling format
#' (VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format) with necessary
#' headers.
#'
#' @param data A data frame containing variant data (`data.frame`). The first
#'   column should represent the chromosome.
#' @param file A character string specifying the output file path with the
#'   `.vcf` or `.vcf.gz` extension (`character string`). If the file has the
#'   `.vcf.gz` extension it will be compressed.
#' @param header An optional custom header for the VCF format (`character`
#'   vector). By default, `NULL` sets a minimal header.
#'
#' @details
#' The function ensures that the first column is named `#CHROM` and writes standard VCF headers.
#' It then appends the variant data in a tab-separated format.
#'
#' @return The function does not return a value but writes a VCF file to the specified path.
#'
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' vcf_data <- data.frame(
#'   CHROM = c("chr1", "chr2"),
#'   POS = c(12345, 67890),
#'   ID = c(".", "."),
#'   REF = c("A", "T"),
#'   ALT = c("G", "C"),
#'   QUAL = c(".", "."),
#'   FILTER = c("PASS", "PASS"),
#'   INFO = c("AF=0.5", "AF=0.3"),
#'   FORMAT = c("GT", "GT"),
#'   SAMPLE = c("0/1", "1/1")
#' )
#' save_as_vcf(vcf_data, "output.vcf")
#' }
#'
#' @export
save_as_vcf <- function(data, file, header = NULL) {
    # Check if input is a data frame
    stopifnot("The `data` must be a data frame." = is.data.frame(data))

    # Check if file extension
    stopifnot("The output file must have a .vcf or .vcf.gz extension." = grepl("\\.vcf(\\.gz)?$", file, ignore.case = TRUE))

    # Ensure the first column is named "#CHROM" (VCF standard)
    colnames(data)[1] <- "#CHROM"

    # VCF headers lines
    if (is.null(header)) {
        header <- c(
            "##fileformat=VCFv4.1",
            "##INFO=<ID=AF,Number=.,Type=Integer>",
            "##FORMAT=<ID=GT,Number=1,Type=String>"
        )
    }

    # Determine if compression is needed
    compress <- grepl("\\.vcf.gz$", file, ignore.case = TRUE)

    # If compressing, write to a temporary .vcf file first
    if (compress) {
        temp_vcf <- tempfile(fileext = ".vcf")
    } else {
        temp_vcf <- file
    }

    # Write VCF headers
    writeLines(header, temp_vcf)
    # Append variant data in tab-separated format without quotes
    suppressWarnings(
        utils::write.table(
            data,
            temp_vcf,
            sep = "\t",
            quote = FALSE,
            append = TRUE,
            row.names = FALSE,
            col.names = TRUE
        )
    )

    # Compress if needed
    if (compress) {
        if (!requireNamespace("Rsamtools", quietly = TRUE)) {
            stop("The 'Rsamtools' package is required for bgzip compression and tabix indexing. Please install it.")
        }

        # Compress with bgzip
        Rsamtools::bgzip(temp_vcf, dest = file, overwrite = TRUE)
        unlink(temp_vcf)  # Remove temporary file

        # Index the VCF with tabix
        Rsamtools::indexTabix(file, format = "vcf")
    }
}



#' Process allele counts results from SCReadCounts format
#'
#' This function reformats SCReadCounts output into a standardized VCF-like data
#' frame and optionally merges unphased/phased genotype information from a
#' VCF-format data frame.
#'
#' @param data A data frame of SCReadCounts output table (`data.frame`).
#' @param vcf An optional VCF-format data frame with unphased/phased genotype
#'   information to add it in output (`data.frame`).
#' @param samplename An optional sample name to prefix cell barcodes
#'   (`character`).
#'
#' @return A cleaned data frame with VCF-like standardized columns:
#' \describe{
#'   \item{`cell`}{Barcodes of cells (`character`).}
#'   \item{`id`}{Variant unique identifier defined as
#'   CHROM_POS_REF_ALT, e.g. "1_920949_C_G" (`character`).}
#'   \item{`CHROM`}{Chromosome in integer format, e.g. 15 (X and Y chromosomes
#'   are not included) (`integer`).}
#'   \item{`POS`}{Position of the variant (1-base positions) (`integer`).}
#'   \item{`REF`}{Reference allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`ALT`}{Alternative allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`RD`}{Reference allele depth/count (`integer`).}
#'   \item{`AD`}{Alternative allele depth/count (`integer`).}
#'   \item{`DP`}{Total depth/count (`integer`).}
#'   \item{`GT`}{Genotype (if `vcf` provided): "0/1" or "1/0" if unphased; "0|1"
#'   or "1|0" if phased. (`character`).}
#' }
#'
#' @note
#' - [SCReadCounts website](https://horvathlab.github.io/NGS/SCReadCounts/)
#' - [Variant Call Format (VCF) format](https://en.wikipedia.org/wiki/Variant_Call_Format)
#'
#' @references
#' \describe{
#'   \item{[SCReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/)}{
#'   Liu, H., Bousounis, P., Movassagh, M., Edwards, N., and Horvath, A.
#'   SCReadCounts: estimation of cell-level SNVs expression from scRNA-seq data.
#'   BMC Genomics 22, 689 (2021).
#'   doi: [10.1186/s12864-021-07974-8](https://doi.org/10.1186/s12864-021-07974-8)
#'   }
#' }
#'
#' @examples
#'
#' # Example data
#' sc_data <- data.frame(
#'   ReadGroup = c("cell1", "cell2"),
#'   CHROM = c(1, 1),
#'   POS = c(1001, 1002),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T"),
#'   SNVCount = c(3, 5),
#'   RefCount = c(7, 5),
#'   GoodReads = c(10, 10)
#' )
#'
#' vcf_data <- data.frame(
#'   CHROM = c(1, 1),
#'   POS = c(1001, 1002),
#'   ID = c(".", "."),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T"),
#'   QUAL = c(".", "."),
#'   FILTER = c(".", "."),
#'   INFO = c(".", "."),
#'   FORMAT = c("GT", "GT"),
#'   sample1 = c("0|1", "1|0"),
#'   stringsAsFactors = FALSE
#' )
#'
#' data <- process_allele(sc_data, vcf = vcf_data, samplename = "sampleA")
#'
#' @export
process_allele <- function(data, vcf = NULL, samplename = NULL) {

    # Check for inputs
    stopifnot("The `data` must be a data frame." = is.data.frame(data))

    required_cols <- c("CHROM", "POS", "REF", "ALT", "ReadGroup", "SNVCount", "RefCount", "GoodReads")
    missing <- setdiff(required_cols, colnames(data))
    if (length(missing) > 0) stop("Missing required columns in `data`: ", paste(missing, collapse = ", "))

    # Create ID column
    data$id <- paste(data$CHROM, data$POS, data$REF, data$ALT, sep = "_")

    if (!is.null(vcf)) {
        stopifnot("The `vcf` must be a data frame." = is.data.frame(vcf))
        # Rename the column 10 to standard (instead of the sample name)
        colnames(vcf)[10] <- "GT"
        # Create ID column
        vcf$id <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = "_")
        # Add genotype information (possibly phased genotype)
        df <- left_join(data, vcf[, c("id", "GT")], by = "id")
    } else {
        df <- data
        df$GT <- "" # temporary column
    }

    # Chromosomes as factor
    df$CHROM <- factor(df$CHROM, sort(unique(df$CHROM)))

    # Rename columns
    colnames(df)[colnames(df) == "ReadGroup"] <- "cell"
    colnames(df)[colnames(df) == "SNVCount"] <- "AD"
    colnames(df)[colnames(df) == "RefCount"] <- "RD"
    colnames(df)[colnames(df) == "GoodReads"] <- "DP"

    # Select final columns
    df <- df[, c("cell", "id", "CHROM", "POS", "REF", "ALT", "RD", "AD", "DP", "GT")]
    # Remove temporary GT column if vcf not provided
    if (is.null(vcf)) {
        df$GT <- NULL
    }

    # Order by position
    df <- df[order(df$CHROM, df$POS), ]

    # Add sample name to the cell barcode if samplename provided
    if (!is.null(samplename)) {
        df$cell <- paste(samplename, df$cell, sep = "_")
    }

    # Make sure unique data only remain
    df <- unique(df)

    return(df)
}
