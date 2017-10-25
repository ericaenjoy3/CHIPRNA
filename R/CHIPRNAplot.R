#' @include CHIPRNAclass.R CHIPRNAconst.R CHIPRNAgeneric.R CHIPRNAmethods.R
#' @rdname plotHeat-methods
setMethod(f = "plotHeat",
  signature = c(tpm.obj = "tpm4plot"),
  definition = function(tpm.obj, pdffout, nms){
    # validation
    dir.create(dirname(pdffout), showWarnings = FALSE)
    # proc
    group <- if(!is.factor(tpm.obj@grp.after)){factor(tpm.obj@grp.after)}
    if (length(levels(group))==1) {
      mat.list <- list(tpm.obj@tpm.rscale);
      names(mat.list) <- levels(group)
    } else {
      mat.list <- split(as.data.frame(tpm.obj@tpm.rscale),group)
    }
    clus.list <- lapply(seq_along(mat.list),
      function(i,nms,mat.list){
        mat <- mat.list[[i]]
        pr_mb <- clusing(mat, pdffout = paste0(nms, "_", names(mat.list)[i], ".pdf"))
        return(paste(names(mat.list)[i], pr_mb, sep = "_"))
      }, nms = nms, mat.list = mat.list)
    clus <- factor(unsplit(clus.list, group))
    fac <- apply(data.frame(clus = clus) %>% mutate(clus=gsub("_","\n", clus)) %>% group_by(clus) %>% mutate(n = n()) %>% ungroup(), 1, paste, collapse = "\n")
    ht_list <- Heatmap(tpm.obj@tpm.rscale, show_row_names = FALSE, show_column_names = TRUE,
      cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, show_column_dend = FALSE,
      heatmap_legend_param = list(title = "", color_bar = "continuous"), clustering_distance_rows = "spearman",
      clustering_method_rows = "average", clustering_distance_columns = "spearman",
      clustering_method_columns = "average", split=fac,gap = unit(3, "mm"))
    # png(pngfout,width=2500*2,height=2500,res=300)
    pdf(pdffout)
    draw(ht_list)
    dev.off()
    return(clus=clus)
  }
)

#' @rdname plotHeatComb-methods
setMethod(f = "plotHeatComb",
  signature = c(chip.obj = "chip", tpm.obj = "tpm4plot"),
  definition = function(chip.obj, tpm.obj, pdffout, option){
    if (option==1) {
      if (class(chip.obj)=="chiploj") {
        chip.obj <- updateChipLOJ(chip.obj, tpm.obj);
      } else {
        chip.obj <- updateChip(chip.obj, tpm.obj);
      }
      mat <- tpm.obj@tpm.rscale
    } else if (option == 2){
      mat <- tpm.obj@tpm.ascale
    } else if (option == 3) {
      mat <- log2(tpm.obj@tpm.val + 0.05)
    } else {
      stop("option is out of bound (1:3)")
    }
    grps <- pull(chip.obj@bed, "clus")
    ngrp <- grps %>% levels() %>% length()
    col1 <- if(ngrp < 3 && ngrp > 1) {
        brewer.pal(3, "Dark2")[1:2]
      } else if ( ngrp > 8 ) {
        brewer.pal(ngrp, "Paired")
      } else {
        brewer.pal(ngrp, "Dark2")
      }
    split <- if (length(unique(as.character(grps)))==1) {NULL} else {grps}
    ht_list <- Heatmap(as.character(grps), col = structure(col1, names = levels(grps)),
      name = "Groups", show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "mm"),
      use_raster = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, split = split,
      combined_name_fun = NULL)
    if (class(chip.obj) == "chiploj") {
      ht_list <- ht_list + Heatmap(chip.obj@binmat,
        col = colorRamp2(range(chip.obj@binmat), c("blue","red")),
        show_row_names = FALSE, show_column_names = TRUE, use_raster = FALSE,
        cluster_columns=FALSE, cluster_rows=FALSE, split = split,
        combined_name_fun = NULL, show_heatmap_legend = FALSE)
    }
    ht_list <- ht_list + Heatmap(mat, show_row_names = FALSE,
      show_column_names = TRUE, cluster_rows = TRUE, show_row_dend = FALSE,
      cluster_columns = FALSE, show_column_dend = FALSE,
      heatmap_legend_param = list(title = "", color_bar = "continuous"),
      clustering_distance_rows = "spearman", clustering_method_rows = "average",
      clustering_distance_columns = "spearman", clustering_method_columns = "average",
      split = split, gap = unit(5, "mm"))
    pdf(pdffout)
    draw(ht_list)
    dev.off()
  }
)

#' @rdname plotBox-methods
setMethod(f = "plotBox",
  signature = c(tpm.obj = "tpm4plot"),
  definition = function(tpm.obj, pdffout, option){
    # validation
    dir.create(dirname(pdffout), showWarnings = FALSE)
    if (option==1) {
      grp <- tpm.obj@grp.after
      tpm <- tpm.obj@tpm.rscale
      ylab <- expression(paste("row standardized ", log[2], "(TPM)"))
    } else if (option==2){
      grp <- tpm.obj@grp.before
      tpm <- tpm.obj@tpm.ascale
      ylab <- expression(paste("all standardized ", log[2], "(TPM)"))
    } else if (option==3) {
      grp <- tpm.obj@grp.before
      tpm <- log2(tpm.obj@tpm.val + 0.05)
      ylab <- expression(paste(log[2], "(TPM)"))
    } else {
      stop("option is out of bound (1:3)")
    }
    ldat <- melt(data.table(grp = grp, tpm), id.vars = "grp")
    browser()
    theme_set(theme_grey(base_size = 15))
    p1 <- ggplot(ldat, aes_(x = ~variable, y = ~value, color = ~grp)) +
      geom_boxplot(na.rm = TRUE, notch = FALSE) +
      labs(x = "",y = ylab)+
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top",
        axis.text.x = element_text(angle = 90))
    p2 <- ggplot(ldat, aes_(x = ~grp, y = ~value, color = ~variable)) +
      geom_boxplot(na.rm = TRUE, notch = FALSE) +
      labs(x = "",y = ylab)+
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top",
        axis.text.x = element_text(angle = 90))
    pdf(pdffout, height = 7*2)
    multiplot(p1,p2,cols = 1)
    dev.off()
  }
)

#' @rdname plotLimma-methods
setMethod(f = "plotLimma",
  signature = c(limma.obj = "limma"),
  definition = function(limma.obj, pdffout) {
    lim <- max(abs(range(limma.obj@dat$pct)))
    p1 <- ggplot(limma.obj@dat, aes_(x = ~cmp, y = ~pct, fill = ~direction))+
    geom_bar(stat="identity") +
    labs(x="",y="Peak With Nearby DiffGenes (%)") +
    theme(legend.position='none') +
    scale_y_continuous(limits = c(-lim, lim),
      breaks = seq(-lim, lim, by = 10),
      labels = abs(seq(-lim, lim, by = 10))) +
    geom_text(aes_(x = ~cmp, y = ~pct, label = abs(~pct)),
      hjust = 0.5, vjust = 0, size = 4, data = limma.obj@dat, inherit.aes = FALSE) +
    coord_flip()
    if (length(unique(as.character(limma.obj@dat$clus)))>1) {
      p1 <- p1 + facet_grid(. ~ clus, scales = "free")
    }
    pdf(pdffout, width = max(7, 3.5 * length(unique(as.character(limma.obj@dat$clus)))), pointsize = 14)
    multiplot(p1, cols = 1)
    dev.off()
  }
)

#' @rdname plotBar-methods
setMethod(f="plotBar",
  signature=c(tpm.obj="tpm4plot"),
  definition=function(tpm.obj,pdffout) {
    tpm.cat <- apply(tpm.obj@tpm.val, 2, function(vec)factor(cut(vec, breaks = c(-Inf,1e-10,5,Inf), right = TRUE, labels = c("0", "<5 TPM", ">5 TPM"))))
    cnt.pct <- melt(data.table(grp = tpm.obj@grp.before, tpm.cat), id.var="grp") %>% group_by(.data$grp, .data$variable, .data$value) %>% tally() %>% group_by(.data$grp) %>% mutate(pct=(100*.data$n)/sum(.data$n))
    theme_set(theme_grey(base_size = 15))
    p1 <- ggplot(cnt.pct, aes_(x = ~variable, fill = ~value,y = ~n))+ geom_bar(stat = "identity") +
      xlab("") + ylab("Count") + facet_grid(. ~ grp) +
      theme(legend.title = element_blank(),panel.spacing = unit(2, "lines"),
        legend.position = "top", axis.text.x = element_text(angle=90))
    p2 <- ggplot(cnt.pct, aes_(x = ~variable,fill = ~value,y = ~pct)) + geom_bar(stat = "identity") +
      xlab("") + ylab("%") + facet_grid(. ~ grp) +
      theme(legend.title = element_blank(), panel.spacing=unit(2, "lines"), legend.position="top",
        axis.text.x = element_text(angle = 90))
    pdf(pdffout)
    multiplot(p1, p2, cols = 2)
    dev.off()
  }
)
