#####rlang .data prevents R CMD check from giving a NOTE about undefined global variables
#' @import ModClusterR
#' @import RColorBrewer
#' @import ggplot2
#' @import dplyr
#' @import methods
#' @import multiplot
#' @import RNA
#' @importFrom data.table fread rbindlist setnames melt data.table .I .N .SD
#' @importFrom rlang .data
#' @importFrom grDevices dev.off pdf
#' @importFrom stats sd
#' @importFrom graphics plot axis abline text par mtext
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom EnrichedHeatmap EnrichedHeatmap '+.AdditiveUnit'
#' @importFrom circlize colorRamp2
#' @importFrom utils read.table write.table setTxtProgressBar
#' @importFrom tibble as_tibble
#' @importFrom readr read_tsv write_tsv
#' @importFrom plyr mapvalues

#' @title A ChIP-seq bed file of 4 columns of chr, start, end and clus
#' @name chip
#' @aliases chip-class
#' @description Store a ChIP-seq peak file
#' @slot bed A \code{tbl} object
#' @exportClass chip
chip <- setClass(
  "chip",
  slots = c(bed = "tbl"),
  validity = function(object) {
    if(nrow(object@bed) < 1) {
      return("Empty data.table was given.")
    }
    if(!identical(colnames(object@bed), c("chr","start","end","clus"))) {
      return("Colnames for bed must be \"chr\",\"start\",\"end\", and \"clus\"")
    }
    return(TRUE)
  }
)

#' @title Store ChIP-seq with binnarized intersections
#' @name chiploj
#' @aliases chiploj-class
#' @description Store ChIP-seq with binnarized intersections
#' @slot binmat A \code{tbl} object
#' @exportClass chiploj
chiploj <- setClass(
    "chiploj",
    slots = c(binmat = "tbl"),
    contains = "chip",
    validity = function(object) {
      if (nrow(object@binmat)<1) {
        return("Eempty binmat was givin.")
      }
      return(TRUE)
    }
)

#' @title A ChIP-seq bed file (4 columns) with 2 additional columns of gid and dist
#' @name gene2peak
#' @aliases gene2peak-class
#' @description The \code{gene2peak} object may be initialized before/after filtering of gid and/or dist.
#' @slot bed A \code{tbl} object
#' @exportClass gene2peak
gene2peak <- setClass(
  "gene2peak",
  slots=c(bed = "tbl"),
  validity=function(object) {
    if(nrow(object@bed) < 1) {
      return("Empty data.table was given.")
    }
    if(!identical(colnames(object@bed), c("chr","start","end","clus","gid","dist"))) {
      return("Colnames for bed must be \"chr\",\"start\",\"end\",\"clus\", \"gid\", and \"dist\"")
    }
    return(TRUE)
  }
)

#' @title Store TPM and related information for plotting
#' @name tpm4plot
#' @aliases tpm4plot-class
#' @rdname tpm4plot-class
#' @description Store TPM and related information for plotting
#' @slot tpm.value A code{tbl} object, representing raw TPM values of either individual or grouped samples'
#' @slot tpm.ascale A code{tbl} object, representing TPM values standardized of all samples
#' @slot tpm.rscale A code{tbl} object, representing TPM values standardized across rows of all samples, zero variance rows were filtered according to var.idx
#' @slot grp.before A character vector the same data value as the 4th column of the \code{bed} slot of the \code{chip} object. The length of \code{grp.before} is identical to the row numbers of \code{tpm.value} and \code{tpm.ascale}
#' @slot grp.after A character vector of \code{grp.before} filtered by \code{var.idx}.
#' @slot vari.idx A logical vector of the same length as \code{grp.before} or the same row numbers as either \code{tpm.val} or \code{tpm.ascale}.
#' @exportClass tpm4plot
tpm4plot <- setClass(
  Class = "tpm4plot",
  slots = c(tpm.val = "tbl", tpm.ascale = "tbl", tpm.rscale = "tbl", info = "tbl", grp.before = "character", grp.after = "character", var.idx = "logical"),
  validity = function(object){
    if (nrow(object@tpm.val) < 1 || nrow(object@tpm.rscale) < 1 || nrow(object@tpm.ascale) < 1 || nrow(object@info) < 1) {
      return("Empty tpm.val, tpm.rscale, tpm.ascale, or info tbl_df was given.")
    }
    if (sum(sapply(colnames(object@tpm.val),is.null)) != 0 || sum(sapply(colnames(object@tpm.rscale),is.null)) != 0) {
      return("column names of tpm.grp or tpm.val were empty.")
    }
    if (length(object@grp.before) < 1 || length(object@grp.after) < 1) {
      return("Empty grp.before or grp.after was given.")
    }
    if (length(object@var.idx) < 1) {
      return("Eempty var.idx was given.")
    }
    if (!identical(dim(object@tpm.val), dim(object@tpm.ascale))){
      return("Matrix dimensions of tpm.val and tpm.ascale do not match. ")
    }
    if (nrow(object@tpm.val) != length(object@grp.before) || length(object@grp.before) != length(object@var.idx) || nrow(object@tpm.val) != nrow(object@info)) {
      return("tpm.val, info, grp.before, and var.idx were not of the same length.")
    }
    if (nrow(object@tpm.rscale) != length(object@grp.after)) {
      return("tpm.rscale and grp.after were not of the same length.")
    }
    return(TRUE)
  }
)

#' @title Linking peak interval files to limma differential statistics
#' @name limma
#' @aliases limma-class
#' @description Linking peak interval files to limma differential statistics
#' @slot dat A \code{tbl} object, representing the linked peak cluster with differential analysis comparison/direction (cmp/direction), percentage/count (pct/cnt) of overlap, a peak count within each cluster, and a count of differential expressed genes/transcripts for such a comparison and direction. Generated by limmaConst function using limma differential analysis file and gene2peak object.
#' @exportClass limma
limma <- setClass(
  "limma",
  slots = c(dat = "tbl"),
  validity = function(object) {
    if(nrow(object@dat) < 1) {
      return("Empty dat was given.")
    }
    if(!identical(colnames(object@dat), c("clus", "cmp", "direction", "pct", "cnt", "peakcnt", "diffcnt"))) {
      return("Colnames for bed must be \"clus\", \"cmp\", \"direction\", \"pct\", \"cnt\",
        \"peakcnt\", and \"diffcnt\"")
    }
    return(TRUE)
  }
)
