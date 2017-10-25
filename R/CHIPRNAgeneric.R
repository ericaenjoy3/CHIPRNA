#' @include CHIPRNAclass.R CHIPRNAconst.R

#' @title Update the \code{bed} slot of the \code{chip} object based upon a \code{gene2peak} object or a \code{tpm4plot} object.
#' @name updateChip
#' @rdname updateChip-methods
#' @description Update the \code{bed} slot of the \code{chip} object based upon a \code{bed} slot of a \code{gene2peak} object or a \code{var.idx} slot of a \code{tpm4plot} object.
#' @param chip.obj A \code{chip} object.
#' @param other.obj A \code{gene2peak} object or a \code{tpm4plot} object.
#' @return A code{chip} object.
#' @export updateChip
setGeneric(name = "updateChip",
  def = function(chip.obj, other.obj) {
    standardGeneric("updateChip")
  }
)

#' @title Update the \code{bed} slot of the \code{chiploj} object based upon a \code{gene2peak} object or a \code{tpm4plot} object.
#' @name updateChipLOJ
#' @rdname updateChipLOJ-methods
#' @description Update the \code{bed} slot of the \code{chiploj} object based upon a \code{gene2peak} object or a \code{tpm4plot} object.
#' @param chiploj.obj A \code{chiploj} object.
#' @param other.obj A \code{gene2peak} object or a \code{tpm4plot} object.
#' @return A code{chiploj} object.
#' @export updateChipLOJ
setGeneric(name = "updateChipLOJ",
  def = function(chiploj.obj, other.obj) {
    standardGeneric("updateChipLOJ")
  }
)

#' @title Update the \code{bed} slot of the \code{gene2peak} object based upon a \code{chip} object.
#' @name updateG2P
#' @rdname updateG2P-methods
#' @description Update the \code{bed} slot of the \code{gene2peak} object based upon a \code{chip} object.
#' @param chip.obj A \code{chip} object.
#' @param gene2peak.obj A \code{gene2peak} object.
#' @return A code{gene2peak} object.
#' @export updateG2P
setGeneric(name = "updateG2P",
  def = function(chip.obj, gene2peak.obj) {
    standardGeneric("updateG2P")
  }
)

#' @include CHIPRNAclass.R
#' @title find nearest gene to peak
#' @name closeGene2Peak
#' @rdname closeGene2Peak-methods
#' @description find nearest gene to peak
#' @param obj A \code{chip} object
#' @param genef A genef (gene or tss) bed formatted file
#' @param genomef A genome file named *.genome
#' @return An \code{gene2peak} object
#' @export closeGene2Peak
setGeneric(name="closeGene2Peak",
  def=function(obj, genef="~/athena/Gencode/mm10/annotation/transcript_tss.bed", genomef="~/athena/Gencode/mm10/sequence/chrNameLength.genome"){
    standardGeneric("closeGene2Peak")
  }
)

#' @title deduplicate the same gene closest to multiple peaks
#' @name deDupGene
#' @rdname deDupGene-methods
#' @description deduplicate the same gene closest to multiple peaks
#' @param obj A \code{gene2peak} object
#' @param distThresh A distance threshold for assigning genes to peaks.
#' @return An updated \code{gene2peak} object
#' @export deDupGene
setGeneric(name="deDupGene",
  def=function(obj, distThresh){
    standardGeneric("deDupGene")
  }
)

#' @title clustering of expressions of genes that are closest to peaks of interest
#' @name clusing
#' @rdname clusing-methods
#' @description clustering of genes that are closest to peaks of interest
#' @slot dat A \code{data.frame} object
#' @slot pdffout A pdf file output file.
#' @return pr_mb predicted cluster membership in integers
#' @export clusing
setGeneric(name="clusing",
  def=function(dat, pdffout) {
    standardGeneric("clusing")
  }
)

#' @title plot heatmap of genes cloest to peaks
#' @name plotHeat
#' @rdname plotHeat-methods
#' @description plot heatmap of genes cloest to peaks
#' @param tpm.obj A code{tpm4plot} object
#' @param pdffout A pdf output file
#' @param nms A character vector of length one, specifying directory and prefix of plotting
#' @return cluster membership for peaks correspond to row scaled TPM
#' @export plotHeat
setGeneric(name="plotHeat",
  def=function(tpm.obj, pdffout, nms) {
    standardGeneric("plotHeat")
  }
)

#' @title Plot binarized heatmap of \code{chiploj} object and gene expression heatmap in one.
#' @name plotHeatComb
#' @rdname plotHeatComb-methods
#' @description Plot binarized heatmap of \code{chiploj} object and gene expression heatmap in one.
#' @param chip.obj A code{chip} object
#' @param tpm.obj A code{tpm4plot} object.
#' @param pdffout A pdf output file.
#' @param option A numeric value specifying the type of tpm values for plotting. 1: row scaled TPM values; 2: all scaled TPM values; 3: log2 TPM values.
#' @return NULL
#' @export plotHeatComb
setGeneric(name="plotHeatComb",
  def=function(chip.obj, tpm.obj, pdffout, option) {
    standardGeneric("plotHeatComb")
  }
)

#' @title Plot peak intersected percentages of differentially expressed genes as a bar plot.
#' @name plotLimma
#' @rdname plotLimma-methods
#' @description Plot peak intersected percentages of differentially expressed genes as a bar plot.
#' @param limma.obj A code{limma} object
#' @param pdffout A pdf output file.
#' @return NULL
#' @export plotLimma
setGeneric(name="plotLimma",
  def=function(limma.obj, pdffout) {
    standardGeneric("plotLimma")
  }
)

#' @title plot boxplots of gene expressions positioned at gene expression samples separated by the identity of peak cluster membership.
#' @name plotBox
#' @rdname plotBox-methods
#' @description plot boxplots of gene expressions positioned at gene expression samples separated by the identity of peak cluster membership.
#' @param tpm.obj An \code{tpm4plot} object.
#' @param pdffout A pdf output file.
#' @return NULL
#' @export plotBox
setGeneric(name="plotBox",
  def=function(tpm.obj, pdffout, option) {
    standardGeneric("plotBox")
  }
)

#' @title plot stacked barplots of counts of genes expressed less than 5 TPM vs more than 5 TPM.
#' @name plotBar
#' @rdname plotBar-methods
#' @description plot stacked barplots of counts of genes expressed less than 5 TPM vs more than 5 TPM.
#' @param tpm.obj An \code{tpm4plot} object.
#' @param pdffout An output pdf file.
#' @return NULL
#' @export plotBar
setGeneric(name="plotBar",
  def=function(tpm.obj, pdffout) {
    standardGeneric("plotBar")
  }
)
