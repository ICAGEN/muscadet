# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions ------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The muscomic class
#'
#' A `muscomic` object encapsulates coverage and allelic data from one
#' single-cell omic dataset used as primary input for muscadet analysis.
#'
#' @slot type Type of single-cell omic. Either ATAC or RNA currently supported
#'   (`character`).
#' @slot label.omic Label to display for the single-cell omic (`character`).
#' @slot coverage Coverage data on features (`list`).
#' @slot allelic Allelic data at variant positions (common SNPs or
#'   individual-specific heterozygous positions) (`list`).
#'
#' @aliases muscomic
#'
#' @seealso
#' Functions related to `muscomic` objects:
#' * [CreateMuscomicObject()] (create `muscomic` objects)
#' * [muscadet-access] (simplified access and assignment methods)
#' * [.DollarNames.muscadet] (autocompletion for `$`)
#' * [`show,muscomic-method`] (`show` method)
#'
#' @importFrom methods setClass
#' @exportClass muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' muscadet_obj$ATAC
#' muscadet_obj$RNA
#'
#' str(muscadet_obj$ATAC, max.level = 2)
#'
methods::setClass(
    "muscomic",
    slots = c(
        type = "character",
        label.omic = "character",
        coverage = "list",
        allelic = "list"
    )
)


#' The muscadet class
#'
#' A `muscadet` object encapsulates data from different single-cell omics as
#' [`muscomic`][muscomic-class] objects as well as downstream analysis results
#' after clustering and CNA calling.
#'
#' @slot omics List of [`muscomic`][muscomic-class] objects, one per single-cell omic
#'   (`list`).
#' @slot bulk.data List of objects containing data from paired bulk sequencing
#'   (`list`).
#' @slot clustering List of objects containing data associated with the
#'   clustering of cells based on coverage log R ratio values (`list`).
#' @slot cnacalling List of objects containing data associated with the calling
#'   of copy number alterations (CNAs) (`list`).
#' @slot genome Reference genome name among: "hg38", "hg19" and "mm10" (`character`).
#'
#' @aliases muscadet
#'
#' @seealso
#' Functions related to `muscadet` objects:
#' * [CreateMuscadetObject()] (create `muscadet` objects)
#' * [muscadet-access] (simplified access and assignment methods)
#' * [.DollarNames.muscadet] (autocompletion for `$`)
#' * [`show,muscadet-method`] (`show` method)
#'
#' @importFrom methods setClass
#' @exportClass muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' muscadet_obj
#'
#' str(muscadet_obj, max.level = 2)
#'
#'
methods::setClass(
    "muscadet",
    slots = c(
        omics = "list",
        bulk.data = "list",
        clustering = "list",
        cnacalling = "list",
        genome = "character"
    )
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create objects functions -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a `muscomic` object
#'
#' Create a [`muscomic`][muscomic-class] object containing coverage and
#' (optionally) allelic count information for a single-cell omic dataset.
#'
#' @param type Type of single cell omic, either "ATAC" or "RNA"
#'   (`character` string).
#' @param mat_counts Matrix of raw counts *cells x features* (`matrix` or
#'   `dgCMatrix`). Rows are cells and columns are features. (they must
#'   correspond to the id column of `features`)
#' @param features Data frame of features (peaks, genes...) coordinates on
#'   genome (`data.frame`). It should contain 4 columns:
#'   \describe{
#'   \item{`CHROM`}{Chromosome names in character format, e.g. "15", "X" (`character`).}
#'   \item{`start`}{Start positions (`integer`).}
#'   \item{`end`}{End positions (`integer`).}
#'   \item{`id`}{Unique identifiers, e.g. gene name "CDH1" or peak identifier
#'   CHROM_start_end "1_1600338_1600838" (`character`). It should match the
#'   feature identifiers as row names of `mat_counts`.}
#' }
#' @param allele_counts Data frame of allele counts at variant positions per
#'   cell (`data.frame`). Variant positions can be either common single
#'   nucleotide polymorphisms (SNPs) positions or individual-specific
#'   heterozygous positions retrieved from bulk sequencing. The data frame format
#'   is based on the Variant Calling Format (VCF), thereby it must contain the
#'   following columns : `cell`, `CHROM`, `POS`, `REF`, `ALT`, `RD`, `AD`. See
#'   [allele_counts()] for details.
#' @param label.omic Label for the single cell omic (`character` string). By
#'   default "scATAC-seq" is used for "ATAC" type and "scRNA-seq" for "RNA"
#'   type.
#' @param label.features Label for features (`character` string). By default
#'   "peaks" is used for "ATAC" type and "genes" for "RNA" type.
#'
#' @return
#' A [`muscomic`][muscomic-class] object.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom stringr str_remove
#' @importFrom methods new
#' @importFrom rlang .data
#'
#' @export
#'
#' @seealso
#' * [CreateMuscadetObject()]
#' * [muscomic-class], [muscadet-class]
#'
#' @examples
#'
#' # Minimal example
#' mat <- Matrix::sparseMatrix(
#'     i = c(1, 1, 2),
#'     j = c(1, 2, 1),
#'     x = c(5, 3, 2),
#'     dims = c(2, 2),
#'     dimnames = list(c("cell1", "cell2"), c("f1", "f2"))
#' )
#'
#' features <- data.frame(
#'     CHROM = c("1", "1"),
#'     start = c(100, 200),
#'     end   = c(150, 250),
#'     id    = c("f1", "f2")
#' )
#'
#' allele <- data.frame(
#'     cell = c("cell1", "cell2"),
#'     CHROM = c("1", "1"),
#'     POS = c(100, 200),
#'     REF = c("A", "G"),
#'     ALT = c("C", "T"),
#'     RD = c(10, 5),
#'     AD = c(3, 1)
#' )
#'
#' muscomic <- CreateMuscomicObject(
#'     type = "ATAC",
#'     mat_counts = mat,
#'     features = features,
#'     allele_counts = allele
#' )
#'
#' # On the example dataset
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = t(mat_counts_atac_tumor),
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' atac
#'
#' # without allele counts data (not required for clustering step)
#' atac2 <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = t(mat_counts_atac_tumor),
#'   features = peaks
#' )
#' atac2
#'
CreateMuscomicObject <- function(type = c("ATAC", "RNA"),
                                 mat_counts,
                                 features,
                                 allele_counts = NULL,
                                 label.omic = NULL,
                                 label.features = NULL) {
    # Check for type of omic
    type.omic <- match.arg(arg = type)

    # Default labels
    if (type.omic == "ATAC") {
        if (is.null(label.omic)) {
            label.omic <- "scATAC-seq"
        }
        if (is.null(label.features)) {
            label.features <- "peaks"
        }
    }
    if (type.omic == "RNA") {
        if (is.null(label.omic)) {
            label.omic <- "scRNA-seq"
        }
        if (is.null(label.features)) {
            label.features <- "genes"
        }
    }

    # COVERAGE -----------------------------------------------------------------

    # Check matrix count class
    stopifnot("`mat_counts` must be either a matrix or a dgCMatrix." =
                  any(class(mat_counts) %in% c("matrix", "dgCMatrix")))
    # Convert to dgCMatrix class objects if matrix
    if (any(class(mat_counts) == "matrix")) {
        mat_counts <- Matrix::Matrix(mat_counts, sparse = TRUE)
    }
    # Check mat counts has row names (cell names)
    stopifnot(
        "The count matrix `mat_counts` must have row names (cell names)." =
            !is.null(rownames(mat_counts))
    )

    # Check coordinates of features
    stopifnot("`features` must be a data frame." = class(features) == "data.frame")
    stopifnot(
        "`features` must be a data frame with 4 columns (CHROM, start, end, id)." = ncol(features) == 4
    )

    # Remove potential NAs
    features <- features[complete.cases(features), ]

    # Force column names of features table
    colnames(features) <- c("CHROM", "start", "end", "id")
    # Modify class if needed
    features$CHROM <- as.character(features$CHROM)
    features$start <- as.integer(features$start)
    features$end <- as.integer(features$end)
    features$id <- as.character(features$id)
    # Remove "chr" in CHROM if necessary
    features$CHROM <- stringr::str_remove(features$CHROM, "chr")
    # Order chromosomes
    chromorder <- c(as.character(1:22), "X", "Y")
    features$CHROM <- ordered(features$CHROM, levels = chromorder[chromorder %in% unique(features$CHROM)])
    # Reorder data
    features <- features[order(features$CHROM, features$start), ]

    # Check col names of matrix matching features id
    stopifnot(
        "Column names of count matrix `mat_counts` must not be NULL and must match `features` id column." =
            !is.null(colnames(mat_counts))
    )
    stopifnot(
        "Column names of count matrix `mat_counts` must match `features` id column." =
            any(features$id %in% colnames(mat_counts))
    )

    # Sort and filter matrix based on provided features
    mat_counts <- mat_counts[, features$id[features$id %in% colnames(mat_counts)], drop = FALSE]

    # Filter features to match mat_counts features
    features <- features[features$id %in% colnames(mat_counts), ]

    # Remove cells with zero counts
    nonzero <- rowSums(mat_counts) > 0
    if (any(!nonzero)) {
        warning(
            "Removing ",
            sum(!nonzero),
            " cells with zero counts: ",
            paste(rownames(mat_counts)[!nonzero], collapse = ", ")
        )
    }
    mat_counts <- mat_counts[nonzero, , drop = FALSE]

    # Coverage data
    coverage <- list(
        counts = list(
            mat = mat_counts,
            coord.features = features,
            label.features = label.features
        )
    )

    # ALLELIC ------------------------------------------------------------------

    if (!is.null(allele_counts)) {
        required_cols <- c("cell", "CHROM", "POS", "REF", "ALT", "RD", "AD")
        stopifnot(
            "The data frame `allele_counts` must contain the required columns: cell, CHROM, POS, REF, ALT, RD, AD." =
                all(required_cols %in% colnames(allele_counts))
        )

        # Remove potential NAs
        allele_counts <- allele_counts[complete.cases(allele_counts), ]
        # Transform to sparse matrices
        allelic <- makeAllelicSparse(allele_counts)

    } else {
        allelic <- list()
    }

    # Create object
    obj <- new(
        Class = "muscomic",
        type = type.omic,
        label.omic = label.omic,
        coverage = coverage,
        allelic = allelic
    )
    return(obj)
}


#' Create a `muscadet` object
#'
#' Create a [`muscadet`][muscadet-class] object.
#'
#' @param omics A list of [`muscomic`][muscomic-class] objects (`list`). The
#'   names of the list will set the names of omics in the final object, if the
#'   list is unamed, the type is taken instead.
#' @param bulk.lrr A data frame containing log R ratio per genomic segments from
#'   bulk sequencing data (`data.frame`). One row per segment and 4 columns
#'   ordered as followed: chromosome (`integer`), start position (`integer`),
#'   end position (`integer`), and Log R ratio value (`numeric`).
#' @param bulk.label Label for bulk data (`character` string).
#' @param genome Reference genome name among: "hg38", "hg19" and "mm10"
#'   (`character` string). "hg38" by default.
#'
#' @return
#' A [`muscadet`][muscadet-class] object.
#'
#' @importFrom stats setNames ave
#' @importFrom methods new
#'
#' @export
#'
#' @note A [`muscadet`][muscadet-class] object can contain several
#'   [`muscomic`][muscomic-class] objects of the same type (`ATAC` or `RNA` for
#'   slot `type`) but they can't have identical labels (slot `label.omic`).
#'
#' @seealso
#' * [CreateMuscomicObject()]
#' * [muscomic-class], [muscadet-class]
#'
#' @examples
#' # Create muscomic objects
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' rna <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_tumor,
#'   allele_counts = allele_counts_rna_tumor,
#'   features = genes
#' )
#' atac_ref <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_ref,
#'   allele_counts = allele_counts_atac_ref,
#'   features = peaks
#' )
#' rna_ref <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_ref,
#'   allele_counts = allele_counts_rna_ref,
#'   features = genes
#' )
#'
#' # Create muscadet objects
#' muscadet <- CreateMuscadetObject(
#'   omics = list(atac, rna),
#'   bulk.lrr = bulk_lrr,
#'   bulk.label = "WGS",
#'   genome = "hg38"
#' )
#' muscadet_ref <- CreateMuscadetObject(
#'   omics = list(atac_ref, rna_ref),
#'   genome = "hg38"
#' )
#'
CreateMuscadetObject <- function(omics,
                                 bulk.lrr = NULL,
                                 bulk.label = NULL,
                                 genome = "hg38") {
  # Check for omic objects
  for (i in seq_along(omics)) {
    stopifnot(
      "`omics` list elements must be of class `muscomic` (use CreateMuscomicObject())." =
        as.character(class(omics[[i]])) == "muscomic"
    )
  }


  # Set names of omics based on type if unamed
  if (is.null(names(omics))) {
    omics <- setNames(omics, unlist(lapply(omics, function(o) {
      o@type
    })))

    omics <- setNames(omics, ave(
      names(omics),
      names(omics),
      FUN = function(i) {
        if (length(i) > 1) {
          paste0(i, "_", seq_along(i))
        } else {
          i
        }
      }
    ))
  }

  # Labels of omics can't be identical
  stopifnot(
    "Identical omic labels found in the `omics` (muscomic object list) provided.
    You can check the labels with `muscomic_obj$label.omic`." =
      !any(duplicated(unlist(
        lapply(omics, function(o) {
          o@label.omic
        })
      )))
  )

  # Check at least one common cell name between all omics
  common_cells <- Reduce(intersect, lapply(omics, Cells))
  stopifnot(
      "No common cells found between omics. Check matrices columns (cell names) for the different omics." =
          length(common_cells) != 0
  )

  # Check bulk data
  if (!is.null(bulk.lrr)) {
      # default names for bulk df columns
      colnames(bulk.lrr) <- c("CHROM", "start", "end", "lrr")
      stopifnot(
          "Label for bulk data (`bulk.label`) is required when `bulk.lrr` is provided." = !is.null(bulk.label)
      )
  }

  # Check for genome
  stopifnot(
    "Genome must be either 'hg38', 'hg19' or 'mm10'." =
      genome %in% c("hg38", "hg19", "mm10")
  )

  # Add bulk log R ratio data
  bulk.data <- list(log.ratio = bulk.lrr, label = bulk.label)

  # Create object
  obj <- new(
    Class = "muscadet",
    omics = omics,
    bulk.data = bulk.data,
    genome = genome
  )
  return(obj)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods to access objects ----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Autocompletion for `$` access on `muscadet` or `muscomic` objects
#'
#' @description Enable autocompletion for `$` access for
#'   [`muscadet`][muscadet-class] or [`muscomic`][muscomic-class] objects. For
#'   `muscadet` objects, it also lists omic datasets contained inside the
#'   `omics` slot.
#'
#' @aliases muscadet-auto
#' @rdname musc-auto
#'
#' @inheritParams utils::.DollarNames
#' @param x A [`muscadet`][muscadet-class] or [`muscomic`][muscomic-class]
#'   object.
#'
#' @return Character vector of matching element names.
#'
#' @importFrom utils .DollarNames
#' @importFrom methods slotNames
#'
#' @order 1
#'
#' @export
#' @method .DollarNames muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#' muscadet_obj$ATAC
#'
".DollarNames.muscadet" <- function(x, pattern = "") {
    # Combine available omics names and slot names
    available <- c(names(x@omics), slotNames(x))
    return(available[grep(pattern, available)])
}

#' @rdname musc-auto
#'
#' @importFrom utils .DollarNames
#' @importFrom methods slotNames
#'
#' @order 2
#'
#' @export
#' @method .DollarNames muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#' muscadet_obj$ATAC$type
#'
".DollarNames.muscomic" <- function(x, pattern = "") {
    available <- c(slotNames(x))
    return(available[grep(pattern, available)])
}


#' @title Access and assignment methods for `muscadet` objects
#'
#' @aliases muscadet-access
#' @rdname musc-access
#'
#' @description Simplified access to omic datasets and slots inside
#'   \code{\link{muscadet}} objects.
#'
#' @inheritParams .DollarNames.muscadet
#' @param i The name of the slot (or omic).
#' @param ... Other arguments (ignored).
#'
#' @return The selected slot or the omic dataset (\code{[muscadet::muscomic()]}
#'   object) for `muscadet` objects. The selected slot for `muscomic` objects.
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [ muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Access to muscadet omics or slots
#' muscadet_obj["ATAC"]
#' muscadet_obj["genome"]
#'
"[.muscadet" <- function(x, i, ...) {
    if (i %in% names(x@omics)) {
        return(x@omics[[i]])
    } else if (i %in% slotNames(x)) {
        return(slot(x, i))
    } else {
        stop(paste("No element named", i, "in the muscadet object."))
    }
}

#' @rdname musc-access
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [ muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Access to muscomic slots
#' muscadet_obj["ATAC"]["label.omic"]
#'
"[.muscomic" <- function(x, i, ...) {
    if (i %in% slotNames(x)) {
        return(slot(x, i))
    } else {
        stop(paste("No element named", i, "in the muscomic object."))
    }
}

#' @rdname musc-access
#'
#' @inheritParams .DollarNames.muscadet
#' @param name The name of the slot (or omic).
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method $ muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Access to muscadet omics or slots
#' muscadet_obj$ATAC
#' muscadet_obj$genome
#'
"$.muscadet" <- function(x, name) {
    if (name %in% names(x@omics)) {
        return(x@omics[[name]])
    } else if (name %in% slotNames(x)) {
        return(slot(x, name))
    } else {
        stop(paste("No element named", name, "in the muscadet object."))
    }
}

#' @rdname musc-access
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method $ muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Access to muscomic slots
#' muscadet_obj$ATAC$label.omic
#'
"$.muscomic" <- function(x, name) {
    if (name %in% slotNames(x)) {
        return(slot(x, name))
    } else {
        stop(paste("No element named", name, "in the muscomic object."))
    }
}

#' @rdname musc-access
#'
#' @description Assign new data in a \code{[muscadet::muscadet()]} or
#'   \code{[muscadet::muscomic()]} object. For `muscadet` objects, the omic datasets in
#'   the `omics` slot can be directly reassigned.
#'
#' @param value The new value to assign.
#'
#' @return The updated \code{[muscadet::muscadet()]} or
#'   \code{[muscadet::muscomic()]} object.
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [<- muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Assign new data in muscadet object
#' muscadet_obj["genome"] <- "hg38"
#'
"[<-.muscadet" <- function(x, i, value) {
    if (i %in% names(x@omics)) {
        x@omics[[i]] <- value
    } else if (i %in% slotNames(x)) {
        slot(x, i) <- value
    } else {
        stop(paste("Cannot assign to", i, "- it does not exist in the muscadet object."))
    }
    return(x)
}


#' @rdname musc-access
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [<- muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Assign new data in muscomic object
#' muscadet_obj["ATAC"]["label.omic"] <- "scATAC-seq"
#'
"[<-.muscomic" <- function(x, i, value) {
    if (i %in% slotNames(x)) {
        slot(x, i) <- value
    } else {
        stop(paste("Cannot assign to", i, "- it does not exist in the muscomic object."))
    }
    return(x)
}

#' @rdname musc-access
#'
#' @export
#' @method $<- muscadet
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Assign new data in muscadet object
#' muscadet_obj$genome <- "hg38"
#'
"$<-.muscadet" <- function(x, i, value) {
    # Use the [<- method for assignment
    x <- `[<-.muscadet`(x, i, value)
    return(x)
}

#' @rdname musc-access
#'
#' @export
#' @method $<- muscomic
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Assign new data in muscomic object
#' muscadet_obj$ATAC$label.omic <- "scATAC-seq"
#'
"$<-.muscomic" <- function(x, i, value) {
    # Use the [<- method for assignment
    x <- `[<-.muscomic`(x, i, value)
    return(x)
}





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods -------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# show -------------------------------------------------------------------------

#' muscomic object overview
#'
#' Overview of a [`muscomic`][muscomic-class] object.
#'
#' @param object A [`muscomic`][muscomic-class] object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @keywords internal
#'
#' @seealso [muscomic-class], [CreateMuscomicObject()]
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = t(mat_counts_atac_tumor),
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' atac
setMethod(
    f = "show",
    signature = signature(object = "muscomic"),
    definition = function(object) {

        get_muscomic_summary <- function(i) {
            list(
                type = i@type,
                label = i@label.omic,
                cells = nrow(i@coverage$counts$mat),
                features_counts = nrow(i@coverage$counts$coord.features),
                feature_lab_counts = i@coverage$counts$label.features,
                features_logratio = nrow(i@coverage$logratio$mat),
                feature_lab_logratio = i@coverage$logratio$label.features,
                vars = nrow(i@allelic$coord.vars)
            )
        }

        summary <- get_muscomic_summary(object)

        line1 <- paste0("A muscomic object of type ", summary$type, ", labelled ", summary$label, ", containing:")
        line_counts <- paste0("counts data: ", summary$features_counts, " features (", summary$feature_lab_counts, ")")
        line_logratio <- paste0("logratio data: ", summary$features_logratio, " features (", summary$feature_lab_logratio, ")")


        if (!is.null(summary$features_logratio)) {
            cat(line1, "\n",
                summary$cells, "cells", "\n",
                line_counts, "\n",
                line_logratio, "\n",
                summary$vars, "variant positions", "\n")

        } else {
            cat(line1, "\n",
                summary$cells, "cells", "\n",
                line_counts, "\n",
                summary$vars, "variant positions", "\n")
        }
    }
)

#' muscadet object overview
#'
#' Overview of a [`muscadet`][muscadet-class] object.
#'
#' @param object A [`muscadet`][muscadet-class] object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @importFrom methods slot
#'
#' @keywords internal
#'
#' @seealso [muscadet-class], [CreateMuscadetObject()]
#'
#' @examples
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' # Overview of the muscadet object
#' show(muscadet_obj)
#'
#' # Overview of the muscomic objects within
#' show(slot(muscadet_obj, "omics"))
#'
#' # Overview of the first muscomic objects within
#' show(slot(muscadet_obj, "omics")[[1]])
#'
setMethod(
    f = "show",
    signature = signature(object = "muscadet"),
    definition = function(object) {

        get_muscomic_summary <- function(i) {
            # Determine which matrix to use (mat.counts or log.ratio)
            matrix_type <- if ("log.ratio" %in% names(i@coverage)) {
                "log.ratio"
            } else {
                "mat.counts"
            }
            list(
                type = i@type,
                label = i@label.omic,
                cells = ncol(i@coverage[[matrix_type]]),
                features = nrow(i@coverage[[matrix_type]]),
                vars = length(unique(i@allelic$table.counts$id)),
                feature_labels = i@coverage$label.features,
                matrix_used = matrix_type
            )
        }

        # Prepare summary data for each muscomic object inside muscadet
        omic_summary <- lapply(object@omics, get_muscomic_summary)
        omic_types <- sapply(omic_summary, function(x) x$type)
        omic_labels <- sapply(omic_summary, function(x) x$label)
        omic_cells <- sapply(omic_summary, function(x) x$cells)
        omic_features <- sapply(omic_summary, function(x) x$features)
        omic_vars <- sapply(omic_summary, function(x) x$vars)
        omic_feature_labels <- sapply(omic_summary, function(x) x$feature_labels)
        omic_matrix_used <- sapply(omic_summary, function(x) x$matrix_used)


        if (!is.null(object@cnacalling$consensus.segs)) {
            cnacall_txt <- paste(
                length(unique(object@cnacalling$clusters)),
                "clusters ;",
                nrow(object@cnacalling$consensus.segs),
                "consensus segments including",
                sum(object@cnacalling$consensus.segs$cna, na.rm = TRUE),
                "CNA segments"
            )
        }

        cat(
            "A muscadet object", "\n",
            length(object@omics), "omics:", paste(names(object@omics), collapse = ", "), "\n",
            "types:", paste(omic_types, collapse = ", "), "\n",
            "labels:", paste(omic_labels, collapse = ", "), "\n",
            "coverage data matrix:", paste(omic_matrix_used, collapse = ", "), "\n",
            "cells:", paste(omic_cells, collapse = ", "), paste0("(common: ", length(Reduce(intersect, Cells(object))), ", total: ", length(Reduce(union, Cells(object))), ")"), "\n", "features:", paste(omic_features, collapse = ", "), "\n",
            "feature labels:", paste(omic_feature_labels, collapse = ", "), "\n",
            "variant positions:", paste(omic_vars, collapse = ", "), "\n",
            "data from paired bulk sequencing:", ifelse(
                is.null(object@bulk.data$label),
                "None", object@bulk.data$label), "\n",
            "clustering:", ifelse(
                length(object@clustering) == 0,
                "None",
                paste("partitions =", paste(names(object@clustering$clusters), collapse = ", "),
                      "; optimal partition =", object@clustering$partition.opt)
            ), "\n",
            "CNA calling:", ifelse(is.null(object@cnacalling$consensus.segs), "None", cnacall_txt), "\n",
            "genome:", object@genome, "\n"
        )
    }
)



# Other methods ----------------------------------------------------------------

#' Methods for [`muscomic`][muscomic-class] and [`muscadet`][muscadet-class]
#' objects
#'
#' Methods to facilitate access to data within the [`muscomic`][muscomic-class]
#' and [`muscadet`][muscadet-class] objects.
#' - `Cells()`: Get cell identifiers (addition of methods for `muscomic` and
#' `muscadet` to [SeuratObject::Cells()]).
#' - `Features()`: Get feature identifiers (addition of methods for `muscomic`
#' and `muscadet` to [SeuratObject::Features()]).
#' - `coordFeatures()`: Get coordinates of features data frames.
#' - `matCounts()`: Get raw count matrices.
#' - `matLogRatio()`: Get log R ratio matrices.
#'
#' The `Cells()`, `Features()`, and `coordFeatures()` functions return information
#' associated with the coverage slot of `muscomic` objects. Their behavior depends
#' on whether log R ratios have been computed:
#' - **Before log R ratio computation:** they return data from the *counts* layer
#'   of the coverage slot.
#' - **After log R ratio computation:** if the coverage slot contains a *logratio*
#'   layer, they return the corresponding data instead of the raw counts.
#'
#' @name muscadet-methods
#' @rdname muscadet-methods
#' @aliases muscomic-methods
#'
#' @param x A [`muscomic`][muscomic-class] or [`muscadet`][muscadet-class]
#'   object.
#' @param ... Other arguments (ignored).
#'
#' @seealso [muscomic-class], [muscadet-class]
#'
#' @examples
#' library("SeuratObject")
#'
#' # Load example muscadet object
#' # data("muscadet_obj")
#'
#' Cells(muscadet_obj) # list of 2 cell names vectors for the 2 omics
#' Cells(muscadet_obj)$ATAC # cell names vector from the omic ATAC
#' Cells(muscadet_obj$ATAC) # cell names vector from the ATAC muscomic object
#' length(Cells(muscadet_obj))
#' length(Cells(muscadet_obj)$ATAC)
#' length(Cells(muscadet_obj$ATAC))
#'
#' Features(muscadet_obj)
#' Features(muscadet_obj)$ATAC
#' Features(muscadet_obj$ATAC)
#'
#' coordFeatures(muscadet_obj)
#' coordFeatures(muscadet_obj)$ATAC
#' coordFeatures(muscadet_obj$ATAC)
#'
#' matCounts(muscadet_obj)
#' matCounts(muscadet_obj)$ATAC
#' matCounts(muscadet_obj$ATAC)
#'
#' matLogRatio(muscadet_obj)
#' matLogRatio(muscadet_obj)$ATAC
#' matLogRatio(muscadet_obj$ATAC)
#'
NULL

#' @rdname muscadet-methods
#' @export
coordFeatures <- function(x) {
  UseMethod(generic = "coordFeatures", object = x)
}

#' @rdname muscadet-methods
#' @export
matCounts <- function(x) {
  UseMethod(generic = "matCounts", object = x)
}

#' @rdname muscadet-methods
#' @export
matLogRatio <- function(x) {
    UseMethod(generic = "matLogRatio", object = x)
}


#' @rdname muscadet-methods
#'
#' @importFrom SeuratObject Cells
#'
#' @return
#' `Cells`:
#' - if `x` is a [muscadet::muscomic()] object: a vector of cell names.
#' - if `x` is a [muscadet::muscadet()] object: a list of cell names vectors,
#' one list element per omic.
#'
#' @method Cells muscomic
#' @export
Cells.muscomic <- function(x, ...) {
    if (!is.null(slot(x, "coverage")[["logratio"]])) {
        rownames(slot(x, "coverage")[["logratio"]][["mat"]])
    } else if (!is.null(slot(x, "coverage")[["counts"]])) {
        rownames(slot(x, "coverage")[["counts"]][["mat"]])
    }
}

#' @rdname muscadet-methods
#'
#' @importFrom SeuratObject Cells
#'
#' @method Cells muscadet
#' @export
Cells.muscadet <- function(x, ...) {
    lapply(slot(x, "omics"), function(omic) {
        if (!is.null(slot(omic, "coverage")[["log.ratio"]])) {
            colnames(slot(omic, "coverage")[["log.ratio"]])
        } else if (!is.null(slot(omic, "coverage")[["mat.counts"]])) {
            colnames(slot(omic, "coverage")[["mat.counts"]])
        }
    })
}

#' @rdname muscadet-methods
#'
#' @importFrom SeuratObject Features
#'
#' @return
#' `Features`:
#' - if `x` is a [muscadet::muscomic()] object: a vector of feature names.
#' - if `x` is a [muscadet::muscadet()] object: a list of feature names vectors,
#' one list element per omic.
#'
#' @method Features muscomic
#' @export
Features.muscomic <- function(x, ...) {
    if (!is.null(slot(x, "coverage")[["logratio"]])) {
        colnames(slot(x, "coverage")[["logratio"]][["mat"]])
    } else if (!is.null(slot(x, "coverage")[["counts"]])) {
        colnames(slot(x, "coverage")[["counts"]][["mat"]])
    }
}

#' @rdname muscadet-methods
#'
#' @importFrom SeuratObject Features
#'
#' @method Features muscadet
#' @export
Features.muscadet <- function(x, ...) {
    lapply(slot(x, "omics"), function(omic) {
        if (!is.null(slot(omic, "coverage")[["log.ratio"]])) {
            rownames(slot(omic, "coverage")[["log.ratio"]])
        } else if (!is.null(slot(omic, "coverage")[["mat.counts"]])) {
            rownames(slot(omic, "coverage")[["mat.counts"]])
        }
    })
}


#' @rdname muscadet-methods
#'
#' @return
#' `coordFeatures`:
#' - if `x` is a [muscadet::muscomic()] object: a data frame of feature coordinates.
#' - if `x` is a [muscadet::muscadet()] object: a list of feature coordinates data frames,
#' one list element per omic.
#'
setMethod(
  f = "coordFeatures",
  signature = signature(x = "muscomic"),
  definition = function(x) {
      if (!is.null(slot(x, "coverage")[["logratio"]])) {
          slot(x, "coverage")[["logratio"]][["coord.features"]]
      } else if (!is.null(slot(x, "coverage")[["counts"]])) {
          slot(x, "coverage")[["counts"]][["coord.features"]]
      }
  }
)



#' @rdname muscadet-methods
#'
setMethod(
    f = "coordFeatures",
    signature = signature(x = "muscadet"),
    definition = function(x) {
        lapply(slot(x, "omics"), function(omic) {
            slot(omic, "coverage")[["coord.features"]]
        })
    }
)

#' @rdname muscadet-methods
#'
#' @return
#' `matCounts`:
#' - if `x` is a [muscadet::muscomic()] object: a [`dgCMatrix`][Matrix::dgCMatrix-class] *features x cells*.
#' - if `x` is a [muscadet::muscadet()] object: a list of [`dgCMatrices`][Matrix::dgCMatrix-class]
#' *features x cells*, one list element per omic.
#'
setMethod(
  f = "matCounts",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(slot(x, "coverage")[["counts"]][["mat"]])
  }
)

#' @rdname muscadet-methods
#'
setMethod(
    f = "matCounts",
    signature = signature(x = "muscadet"),
    definition = function(x) {
        lapply(slot(x, "omics"), function(omic) {
            slot(omic, "coverage")[["mat.counts"]]
        })
    }
)

#' @rdname muscadet-methods
#'
#' @return
#' `matLogRatio`:
#' - if `x` is a [muscadet::muscomic()] object: a `matrix` *features x cells*.
#' - if `x` is a [muscadet::muscadet()] object: a list of `matrices` *features x cells*,
#' one list element per omic.
#'
setMethod(
    f = "matLogRatio",
    signature = signature(x = "muscomic"),
    definition = function(x) {
        return(slot(x, "coverage")[["logratio"]][["mat"]])
    }
)

#' @rdname muscadet-methods
#'
setMethod(
    f = "matLogRatio",
    signature = signature(x = "muscadet"),
    definition = function(x) {
        lapply(slot(x, "omics"), function(omic) {
            slot(omic, "coverage")[["log.ratio"]]
        })
    }
)
