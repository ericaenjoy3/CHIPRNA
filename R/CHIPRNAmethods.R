#' @include CHIPRNAclass.R CHIPRNAgeneric.R
#' @rdname closeGene2Peak-methods
setMethod(f="closeGene2Peak",
  signature="chip",
  definition=function(obj,genef,genomef){
    before<-tempfile()
    after<-tempfile()
    write_tsv(obj@bed,path=before,col_names=FALSE)
    command<-paste("sort -k1,1V -k2,2n", before, "| bedtools closest -a stdin -b", genef, "-g", genomef, "-d | awk 'BEGIN{OF=\"\\t\";OFS=\"\\t\"} {print $1,$2,$3,$4,$(NF-3),$NF}' | sort -k1,1V -k2,2n -k5,5V | uniq >",after,sep=" ")
    cat(command,"\n")
    Sys.setenv(SHELL="/bin/bash")
    try(system(command))
    newObj<-new("gene2peak",bed=read_tsv(after,col_names=c("chr","start","end","clus","gid","dist")) %>% select(1:6))
    unlink(before);unlink(after);
    return(newObj)
  }
)

#' @rdname deDupGene-methods
setMethod(f="deDupGene",
  signature="gene2peak",
  definition=function(obj,distThresh){
    dat<-obj@bed %>% filter(dist!=-1, dist<=distThresh)
    # the same gene dected multiple times
    norep<- dat %>% group_by(.data$gid) %>% filter(n()==1)
    # repeititon of gid within the same cluster
    repwithin <- dat %>% group_by(.data$gid) %>% filter(n()>1) %>% filter(length(unique(.data$clus))==1) %>% sample_n(1)
    # reptition of gid within a different cluster;
    # keep only the gene, if the min distance is only present in 1 cluster;
    # remove the gene, if the min distance is present in multiple clusters.
    repbet <- dat %>% group_by(.data$gid) %>% filter(n()>1) %>% filter(length(unique(.data$clus))>1, near(.data$dist,min(.data$dist))) %>% filter(n()==1)
    # if such gene is present in multiple clusters, remove such genes
    filter.dat<-do.call(rbind,list(norep,repwithin,repbet)) %>% ungroup()
    obj@bed<-filter.dat
    return(obj)
  }
)

#' @rdname clusing-methods
setMethod(f = "clusing",
  signature=c(dat = "data.frame", pdffout="character"),
  definition=function(dat, pdffout){
    if(!is.null(dirname(pdffout)) || dirname(pdffout)==".") dir.create(dirname(pdffout), showWarnings = FALSE)
    pdf(pdffout)
    opt <- Optimal_Clusters_KM(dat, max_clusters = min(10,ncol(dat)), plot_clusters=TRUE, criterion = 'distortion_fK', fK_threshold = 0.85, initializer= 'optimal_init', tol_optimal_init = 0.2)
    dev.off()
    km_mb <- MiniBatchKmeans(dat, clusters = opt, batch_size = 20, num_init = 5, max_iters = 100,
    init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,verbose = F)
    pr_mb <- predict_MBatchKMeans(dat, km_mb$centroids)
    return(pr_mb)
  }
)

#' @rdname plotHeat-methods
setMethod(f = "plotHeat",
  signature=c(tpm.obj = "tpm4plot"),
  definition=function(tpm.obj, pdffout, nms){
    # validation
    dir.create(dirname(pdffout), showWarnings = FALSE)
    # proc
    group <- if(!is.factor(tpm.obj@grp.after)){factor(tpm.obj@grp.after)}
    if (length(levels(group))==1) {
      mat.list<-list(tpm.obj@tpm.rscale);
      names(mat.list)<-levels(group)
    } else {
      mat.list<-split(as.data.frame(tpm.obj@tpm.rscale),group)
    }
    clus.list<-lapply(seq_along(mat.list),
      function(i,nms,mat.list){
        mat <- mat.list[[i]]
        pr_mb <- clusing(mat, pdffout = paste0(nms, "_", names(mat.list)[i], ".pdf"))
        return(paste(names(mat.list)[i], pr_mb, sep = "_"))
      }, nms = nms, mat.list = mat.list)
    clus <- factor(unsplit(clus.list, group))
    fac <- apply(data.frame(clus = clus) %>% mutate(clus=gsub("_","\n", clus)) %>% group_by(clus) %>% mutate(n = n()) %>% ungroup(), 1, paste, collapse = "\n")
    ht_list <- Heatmap(tpm.obj@tpm.rscale, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, show_row_dend = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, heatmap_legend_param = list(title = "", color_bar = "continuous"), clustering_distance_rows = "spearman", clustering_method_rows = "average", clustering_distance_columns = "spearman", clustering_method_columns = "average", split=fac,gap = unit(3, "mm"))
    # png(pngfout,width=2500*2,height=2500,res=300)
    pdf(pdffout)
    draw(ht_list)
    dev.off()
    return(clus=clus)
  }
)

#' @rdname plotBox-methods
setMethod(f = "plotBox",
  signature=c(tpm.obj = "tpm4plot"),
  definition=function(tpm.obj, pdffout){
    # validation
    dir.create(dirname(pdffout), showWarnings = FALSE)
    # proc
    ldat <- melt(data.table(grp = tpm.obj@grp.before,tpm.obj@tpm.ascale),id.vars = "grp")
    theme_set(theme_grey(base_size = 15))
    p1 <- ggplot(ldat, aes(variable, value, color = grp)) +
    #geom_jitter(alpha=I(1/4), aes(color=grp),na.rm=TRUE) +
    # geom_violin(alpha=I(1/4), aes(fill=grp),na.rm=TRUE) +
    geom_boxplot(na.rm = TRUE, notch = FALSE) +
    labs(x = "",y = expression(paste("standardized ", log[2], "(TPM)")))+
    theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top", axis.text.x = element_text(angle = 90))
    # png(pngfout,width=min(3000,length(unique(ldat$variable))*800),height=3000,res=300)
    pdf(pdffout)
    multiplot(p1,cols = 1)
    dev.off()
  }
)

#' @rdname plotBar-methods
setMethod(f="plotBar",
  signature=c(tpm.obj="tpm4plot"),
  definition=function(tpm.obj,pdffout) {
    tpm.cat <- apply(tpm.obj@tpm.val, 2, function(vec)factor(cut(vec, breaks = c(-Inf,1e-10,5,Inf), right = TRUE, labels = c("0", "<5 TPM", ">5 TPM"))))
    cnt.pct <- melt(data.table(grp = tpm.obj@grp.before, tpm.cat), id.var="grp") %>% group_by(grp, variable, value) %>% tally() %>% group_by(grp) %>% mutate(pct=(100*n)/sum(n))
    theme_set(theme_grey(base_size = 15))
    p1 <- ggplot(cnt.pct,aes(x = variable,fill = value,y = n))+ geom_bar(stat = "identity") + xlab("") + ylab("Count") + facet_grid(. ~ grp) + theme(legend.title = element_blank(),panel.spacing = unit(2, "lines"),legend.position = "top", axis.text.x = element_text(angle=90))
    p2 <- ggplot(cnt.pct,aes(x = variable,fill = value,y = pct)) + geom_bar(stat = "identity") + xlab("") + ylab("%") +
    facet_grid(. ~ grp)+theme(legend.title=element_blank(),panel.spacing=unit(2, "lines"),legend.position="top", axis.text.x = element_text(angle = 90))
    pdf(pdffout)
    multiplot(p1, p2, cols = 2)
    dev.off()
  }
)
