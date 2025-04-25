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
#' A `muscadet` object encapsulates data from different single-cell omics
#' as \code{\link{muscomic}} objects as well as downstream analysis results after
#' clustering and CNA calling.
#'
#' @slot omics List of \code{\link{muscomic}} objects, one per single-cell omic
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

#' Create a muscomic object
#'
#' Create a \code{\link{muscomic}} object.
#'
#' @param type Type of single cell omic, either "ATAC" or "RNA"
#'   (`character` string).
#' @param mat_counts Matrix of raw counts *features x cells* (`matrix` or
#'   `dgCMatrix`). Rows are features (they must correspond to the id column of
#'   `features`), and columns are cells.
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
#'   heterozygous positions retrieved by bulk sequencing. The data frame format
#'   is based on the Variant Calling Format (VCF), thereby it must contain the
#'   following columns : `cell`, `id`, `CHROM`, `POS`, `REF`, `ALT`, `RD`, `AD`,
#'   `DP`, `GT`. See [allele_counts] for details.
#' @param label.omic Label for the single cell omic (`character` string). By
#'   default "scATAC-seq" is used for "ATAC" type and "scRNA-seq" for "RNA"
#'   type.
#' @param label.features Label for features (`character` string). By default
#'   "peaks" is used for "ATAC" type and "genes" for "RNA" type.
#'
#' @return
#' A \code{\link{muscomic}} object.
#'
#' @import dplyr
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix summary
#' @importFrom stringr str_remove
#' @importFrom gtools mixedsort
#' @importFrom methods new
#' @importFrom rlang .data
#'
#' @export
#'
#' @seealso
#' * [CreateMuscadetObject()]
#' * \code{\link{muscomic-class}}, \code{\link{muscadet-class}}
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' atac
#'
#' rna <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_tumor,
#'   allele_counts = allele_counts_rna_tumor,
#'   features = genes
#' )
#' rna
#'
#' atac_ref <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_ref,
#'   allele_counts = allele_counts_atac_ref,
#'   features = peaks
#' )
#'
#' rna_ref <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_ref,
#'   allele_counts = allele_counts_rna_ref,
#'   features = genes
#' )
#' rna_ref
#'
#' # without allele counts data (not required for clustering step)
#' atac2 <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
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

  # Check count matrix
  stopifnot(
    "`mat_counts` must be either a matrix or a dgCMatrix." =
      any(class(mat_counts) %in% c("matrix", "dgCMatrix"))
  )
  # Convert to dgCMatrix class objects if matrix
  if (any(class(mat_counts) == "matrix")) {
    mat_counts <- Matrix::Matrix(mat_counts)
  }

  # Check coordinates of features
  stopifnot("`features` must be a data frame." = class(features) == "data.frame")
  colnames(features) <- c("CHROM", "start", "end", "id")
  # Modify class if needed
  features$CHROM <- as.character(features$CHROM)
  features$start <- as.integer(features$start)
  features$end <- as.integer(features$end)
  features$id <- as.character(features$id)
  # Remove "chr" if necessary
  features$CHROM <- stringr::str_remove(features$CHROM, "chr")
  # Ordered chromosomes
  features$CHROM <- ordered(features$CHROM, levels = gtools::mixedsort(unique(features$CHROM)))
  # Reorder data
  features <- features[order(features[, "start"]), ]
  features <- features[order(features[, "CHROM"]), ]

  # Check row names of matrix matching features id
  stopifnot("Row names of count matrix `mat_counts` must not be NULL and must match `features` id column." =
                !is.null(rownames(mat_counts)))
  stopifnot("Row names of count matrix `mat_counts` must match `features` id column." =
                any(features$id %in% rownames(mat_counts)))

  # Sort and filter matrix based on provided features
  mat_counts <- mat_counts[features$id[features$id %in% rownames(mat_counts)], ]

  # Filter matrix to remove cells with zero counts
  zero_count_cells <- colnames(mat_counts)[colSums(mat_counts) == 0]
  if (length(zero_count_cells) > 0) {
      warning("The following cells are removed due to zero counts in `mat_counts`: ",
              paste(zero_count_cells, collapse = ", "))
  }
  mat_counts <- mat_counts[, colSums(mat_counts) > 0, drop = FALSE]

  # Get unique index value for features
  coord.df <- dplyr::mutate(features, index = as.numeric(1:nrow(features)))
  # Generate table of raw counts by converting matrix to summary data frame
  table_counts <- as.data.frame(Matrix::summary(mat_counts)) %>%
      dplyr::mutate(cell = colnames(mat_counts)[.data$j], omic = type.omic) %>%
      dplyr::rename(index = "i", DP = "x") %>%
      dplyr::left_join(coord.df, by = "index") %>%
      dplyr::mutate(POS = round((.data$start + .data$end) / 2, 0)) %>% # POS as center of the feature
      dplyr::select(c("omic", "cell", "id", "CHROM", "POS", "DP"))
  # Ordered chromosomes
  table_counts$CHROM <- ordered(table_counts$CHROM, levels = gtools::mixedsort(unique(table_counts$CHROM)))
  # Reorder data
  table_counts <- table_counts[order(table_counts[, "POS"]), ]
  table_counts <- table_counts[order(table_counts[, "CHROM"]), ]

  # Coverage data
  coverage <- list(
    mat.counts = mat_counts,
    table.counts = table_counts,
    coord.features = features,
    label.features = label.features
  )

  # Allelic data
  if (!is.null(allele_counts)) {

      # Remove "chr" if necessary
      allele_counts$CHROM <- stringr::str_remove(allele_counts$CHROM, "chr")
      # Ordered chromosomes
      allele_counts$CHROM <- ordered(allele_counts$CHROM, levels = gtools::mixedsort(unique(allele_counts$CHROM)))
      # Reorder data
      allele_counts <- allele_counts[order(allele_counts[, "POS"]), ]
      allele_counts <- allele_counts[order(allele_counts[, "CHROM"]), ]

      allelic <- list(table.counts = data.frame(omic = type.omic, allele_counts))

  } else {

      allelic <- list(table.counts = NULL)
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


#' Create a muscadet object
#'
#' Create a \code{\link{muscadet}} object.
#'
#' @param omics A list of \code{\link{muscomic}} objects (`list`). The names of
#'   the list will set the names of omics in the final object, if the list is
#'   unamed, the type is taken instead.
#' @param bulk.lrr A data frame containing log R ratio per genomic segments from
#'   bulk sequencing data (`data.frame`). One row per segment and 4 columns
#'   ordered as followed: chromosome (`integer`), start position (`integer`),
#'   end position (`integer`), and Log R ratio value (`numeric`).
#' @param bulk.label Label for bulk data (`character` string).
#' @param genome Reference genome name among: "hg38", "hg19" and "mm10"
#'   (`character` string). "hg38" by default.
#'
#' @return
#' A \code{\link{muscadet}} object.
#'
#' @importFrom stats setNames ave
#' @importFrom methods new
#'
#' @export
#'
#' @note A \code{\link{muscadet}} object can contain several
#'   \code{\link{muscomic}} objects of the same type (`ATAC` or `RNA` for slot
#'   `type`) but they can't have identical labels (slot `label.omic`).
#'
#' @seealso
#' * [CreateMuscomicObject()]
#' * \code{\link{muscomic-class}}, \code{\link{muscadet-class}}
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
#' @description Enable autocompletion for `$` access for \code{\link{muscadet}}
#'   or \code{\link{muscomic}} objects. For `muscadet` objects, it also lists
#'   omic datasets contained inside the `omics` slot.
#'
#' @rdname musc-auto
#'
#' @inheritParams utils::.DollarNames
#' @param x A \code{\link{muscadet}} or \code{\link{muscomic}} object.
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
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' @return The selected slot or the omic dataset (\code{\link{muscomic}} object)
#'   for `muscadet` objects. The selected slot for `muscomic` objects.
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [ muscadet
#'
#' @examples
#' # Load example muscadet object
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' @description Assign new data in a \code{\link{muscadet}} or
#'   \code{\link{muscomic}} object. For `muscadet` objects, the omic datasets in
#'   the `omics` slot can be directly reassigned.
#'
#' @param value The new value to assign.
#'
#' @return The updated \code{\link{muscadet}} or \code{\link{muscomic}} object.
#'
#' @importFrom methods slot
#' @importFrom methods slotNames
#'
#' @export
#' @method [<- muscadet
#'
#' @examples
#' # Load example muscadet object
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' data(muscadet_obj)
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
#' Overview of a \code{\link{muscomic}} object.
#'
#' @param object A \code{\link{muscomic}} object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @keywords internal
#'
#' @seealso \code{\link{muscomic-class}}, [CreateMuscomicObject()]
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' atac
setMethod(
    f = "show",
    signature = signature(object = "muscomic"),
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

        summary <- get_muscomic_summary(object)
        cat(
            "A muscomic object of type", summary$type,
            "labelled", summary$label, "containing:", "\n",
            summary$matrix_used, "coverage data matrix", "\n",
            summary$cells, "cells", "\n",
            summary$features, "features:", summary$feature_labels, "\n",
            summary$vars, "variant positions", "\n"
        )
    }
)

#' muscadet object overview
#'
#' Overview of a \code{\link{muscadet}} object.
#'
#' @param object A \code{\link{muscadet}} object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @importFrom methods slot
#'
#' @keywords internal
#'
#' @seealso \code{\link{muscadet-class}}, [CreateMuscadetObject()]
#'
#' @examples
#' # Load example muscadet object
#' data(muscadet_obj)
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
                "k =",
                length(unique(object@cnacalling$clusters)),
                ";",
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

#' Methods for \code{\link{muscomic}} and \code{\link{muscadet}} objects
#'
#' Methods to facilitate access to data within the \code{\link{muscomic}} and
#' \code{\link{muscadet}} objects.
#' - `Cells()`: Get cell identifiers (addition of methods for `muscomic` and
#'    `muscadet` to [SeuratObject::Cells()]).
#' - `Features()`: Get feature identifiers (addition of methods for `muscomic`
#'    and `muscadet` to [SeuratObject::Features()]).
#' - `coordFeatures()`: Get coordinates of features data frames.
#' - `matCounts()`: Get raw count matrices.
#' - `matLogRatio()`: Get log R ratio matrices.
#'
#' @name muscadet-methods
#' @rdname muscadet-methods
#' @aliases muscomic-methods
#'
#' @param x A \code{\link{muscomic}} or \code{\link{muscadet}} object.
#' @param ... Other arguments (ignored).
#'
#' @seealso \code{\link{muscomic-class}} \code{\link{muscadet-class}}
#'
#' @examples
#' library("SeuratObject")
#' # Load example muscadet object
#' data(muscadet_obj)
#'
#' Cells(muscadet_obj) # list of 2 cell names vectors for the 2 omics
#' Cells(muscadet_obj)$ATAC # cell names vector from the omic ATAC
#' Cells(muscadet_obj$ATAC) # cell names vector from the ATAC muscomic object
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
#' - if `x` is a `muscomic` object: A vector of cell names.
#' - if `x` is a `muscadet` object: A list of vectors of cell names, one list
#' element per omic.
#'
#' @method Cells muscomic
#' @export
Cells.muscomic <- function(x, ...) {
    if (!is.null(slot(x, "coverage")[["log.ratio"]])) {
        cells <- colnames(slot(x, "coverage")[["log.ratio"]])
    } else if (!is.null(slot(x, "coverage")[["mat.counts"]])) {
        cells <- colnames(slot(x, "coverage")[["mat.counts"]])
    }
    return(cells)
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
            cells <- colnames(slot(omic, "coverage")[["log.ratio"]])
        } else if (!is.null(slot(omic, "coverage")[["mat.counts"]])) {
            cells <- colnames(slot(omic, "coverage")[["mat.counts"]])
        }
        return(cells)
    })
}

#' @rdname muscadet-methods
#'
#' @importFrom SeuratObject Features
#'
#' @return
#' `Features`:
#' - if `x` is a `muscomic` object: A vector of feature names.
#' - if `x` is a `muscadet` object: A list of vectors of feature names, one list
#' element per omic.
#'
#' @method Features muscomic
#' @export
Features.muscomic <- function(x, ...) {
    if (!is.null(slot(x, "coverage")[["log.ratio"]])) {
        features <- rownames(slot(x, "coverage")[["log.ratio"]])
    } else if (!is.null(slot(x, "coverage")[["mat.counts"]])) {
        features <- rownames(slot(x, "coverage")[["mat.counts"]])
    }
    return(features)
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
            features <- rownames(slot(omic, "coverage")[["log.ratio"]])
        } else if (!is.null(slot(omic, "coverage")[["mat.counts"]])) {
            features <- rownames(slot(omic, "coverage")[["mat.counts"]])
        }
        return(features)
    })
}


#' @rdname muscadet-methods
#'
#' @return
#' `coordFeatures`:
#' - if `x` is a `muscomic` object: A data frame of feature coordinates.
#' - if `x` is a `muscadet` object: A list of data frames of feature
#' coordinates, one list element per omic.
#'
setMethod(
  f = "coordFeatures",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(slot(x, "coverage")[["coord.features"]])
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "coordFeatures",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(slot(omic, "coverage")[["coord.features"]])
    })
  }
)

#' @rdname muscadet-methods
#'
#' @return
#' `matCounts`:
#' - if `x` is a `muscomic` object: A \code{\link{dgCMatrix-class}} *features x cells*.
#' - if `x` is a `muscadet` object: A list of \code{\link{dgCMatrix-class}} *features x cells*,
#' one list element per omic.
#'
setMethod(
  f = "matCounts",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(slot(x, "coverage")[["mat.counts"]])
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "matCounts",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(slot(omic, "coverage")[["mat.counts"]])
    })
  }
)

#' @rdname muscadet-methods
#'
#' @return
#' `matLogRatio`:
#' - if `x` is a `muscomic` object: A \code{\link{matrix}} *features x cells*.
#' - if `x` is a `muscadet` object: A list of \code{\link{matrix}} *features x cells*,
#' one list element per omic.
#'
setMethod(
    f = "matLogRatio",
    signature = signature(x = "muscomic"),
    definition = function(x) {
        return(slot(x, "coverage")[["log.ratio"]])
    }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "matLogRatio",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(slot(omic, "coverage")[["log.ratio"]])
    })
  }
)
