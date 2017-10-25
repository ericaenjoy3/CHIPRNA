#!/usr/bin/env Rscript

###
# Given a query bed file, the query and target config files for matrix series
#
# EL
# Created on 5 July, 2017
# Last modified on 12 July, 2017
###

suppressPackageStartupMessages(library("argparse"))
library(GOplot)
library(CHIPRNA)

parser <- ArgumentParser()
parser$add_argument("--chip", type = "character", nargs = "+", required = TRUE, help = "one or multiple chipseq bed files")
parser$add_argument("--sms", type = "character", nargs = "*", help = "if multiple chip-seq files supplied, sample names must bed supplied for each bed file")
parser$add_argument("--tpm", type = "character", nargs = "*", help = "TPM file(s)")
parser$add_argument("--limma", type = "character", nargs = "*", help = "Limma file(s)")
parser$add_argument("--plot", type = "character", nargs = "+", required = TRUE, help = "heatmap, boxplot and/or barplots")
parser$add_argument("--distThresh", type = "integer", nargs = 1, required = TRUE, help = "distance of gene to peak or peak to gene")
parser$add_argument("--nms", type = "character", nargs = 1, required = TRUE, help = "output directory with file pattern")
parser$add_argument("--gene", action  = "store_true", help = "cloest peaks to genes")
parser$add_argument("--loj", action = "store_true", help = "LOJ bed input")
args <- parser$parse_args()

if (!args$loj) {
  chip.obj <- chipConst(chipF = args$chip, sms = NULL)
} else {
  chip.obj <- chiplojConst(chipF = args$chip, reverse = TRUE)
}

gene2peak.obj <- closeGene2Peak(chip.obj)
gene2peak.obj <- deDupGene(gene2peak.obj, args$distThresh)

if (args$loj) {
  chip.obj <- updateChipLOJ(chip.obj, gene2peak.obj)
} else {
  chip.obj <- updateChip(chip.obj, gene2peak.obj)
}

gene2peak.obj <- updateG2P(chip.obj, gene2peak.obj)

tpm.obj <- tpmConst(tpm = args$tpm, gene2peak.obj)


plotHeatComb(chip.obj, tpm.obj, pdffout = paste0(args$nms, "_plotHeatComb_rscale.pdf"), option = 1)
plotHeatComb(chip.obj, tpm.obj, pdffout = paste0(args$nms, "_plotHeatComb_ascale.pdf"), option = 2)
plotHeatComb(chip.obj, tpm.obj, pdffout = paste0(args$nms, "_plotHeatComb_logtpm.pdf"), option = 3)

# clus <- plotHeat(tpm.obj, pdffout = paste0(args$nms,"_plotHeat_rscale.pdf"), nms = args$nms, option = 1)
# clus <- plotHeat(tpm.obj, pdffout = paste0(args$nms,"_plotHeat_rscale.pdf"), nms = args$nms, option = 2)
# clus <- plotHeat(tpm.obj, pdffout = paste0(args$nms,"_plotHeat_rscale.pdf"), nms = args$nms, option = 3)


if (!is.null(args$limma)) {
  limma.list <- limmaConst(args$limma, gene2peak.obj)
  limma.obj <- limma.list$limma.obj
  limma_dat <- limma.list$dat
  plotLimma(limma.obj, pdffout = paste0(args$nms,"_plotLimma.pdf"))
}

plotBox(tpm.obj, pdffout = paste0(args$nms, "_plotBox_rscale.pdf"), option = 1)
plotBox(tpm.obj, pdffout = paste0(args$nms, "_plotBox_ascale.pdf"), option = 2)
plotBox(tpm.obj, pdffout = paste0(args$nms, "_plotBox_log2tpm.pdf"), option = 3)

# heatmap

#
gset.obj <- selDB(major = "C2.CP", minor = "Reactome", type= "symbols", species = "mouse")
gene=factor(sapply(strsplit(gene2peak.obj@bed$gid, "|", fixed=TRUE), function(vec)vec[2])[tpm.obj@var.idx])
idx.rm <- duplicated(gene)
if (any(idx.rm)) {
  gene <- gene[!idx.rm]
  clus <- clus[!idx.rm]
}
gclus.obj<-new("gclus", tbl=as_data_frame(data.frame(gene = gene, clus = clus)))

# GO
res.list <- GO(gclus.obj, gset.obj, filterP = 0.05, filterOR = TRUE)
go_set.obj <- res.list$go_set.obj
go_res.obj <- res.list$go_res.obj
write_GO(go_set.obj, go_res.obj, args$nms)
simi(go_set.obj, go_res.obj, args$nms)




plotBar(tpm.obj, pngfout = paste0(args$nms, "_plotBar.png"))

# constitute closest peaks and genes
# to be worked on
#closeGene2Peak<-function(chipobj,genef,genomef){
#}

# barplots of gene expressions by chip-seq groups or subsetted by diffiferential expressed genes

# heatmap of gene expressions separated by chip-seq groups subsetted by diffiferential expressed genes


# or heatmap of gene expressions separated by diffiferential expressed genes subsetted by chip-seq groups


# box plots of gene expression separated by chip-seq groups subsetted by diffiferential gexpressed genes
# or box plots of gene expression separated by diffiferential gexpressed genes
