#' @include CHIPRNAclass.R CHIPRNAconst.R CHIPRNAgeneric.R
#' @rdname closeGene2Peak-methods
setMethod(f = "updateChip",
  signature = c(chip.obj = "chip", other.obj = "gene2peak"),
  definition = function(chip.obj, other.obj) {
    str_before <- apply(select(chip.obj@bed, chr:clus), 1, paste, collapse = "-")
    str_after <- apply(select(other.obj@bed, chr:clus), 1, paste, collapse = "-")
    idx <- which(str_before %in% str_after)
    return(new("chip", bed = chip.obj@bed[idx,]))
  }
)

#' @rdname updateChip-methods
setMethod(f = "updateChip",
  signature = c(chip.obj = "chip", other.obj = "tpm4plot"),
  definition = function(chip.obj, other.obj) {
    idx <- tpm.obj@var.idx
    return(new("chip", bed = chip.obj@bed[idx,]))
  }
)

#' @rdname updateChipLOJ-methods
setMethod(f = "updateChipLOJ",
  signature = c(chiploj.obj = "chiploj", other.obj = "gene2peak"),
  definition = function(chiploj.obj, other.obj) {
    str_before <- apply(select(chiploj.obj@bed, chr:clus), 1, paste, collapse = "-")
    str_after <- apply(select(other.obj@bed, chr:clus), 1, paste, collapse = "-")
    idx <- which(str_before %in% str_after)
    return(new("chiploj", bed = chiploj.obj@bed[idx,], binmat = chiploj.obj@binmat[idx,]))
  }
)

#' @rdname updateChipLOJ-methods
setMethod(f = "updateChipLOJ",
  signature = c(chiploj.obj = "chiploj", other.obj = "tpm4plot"),
  definition = function(chiploj.obj, other.obj) {
    idx <- tpm.obj@var.idx
    return(new("chiploj", bed = chiploj.obj@bed[idx,], binmat = chiploj.obj@binmat[idx,]))
  }
)

#' @rdname updateG2P-methods
setMethod(f = "updateG2P",
  signature = c(chip.obj = "chip", gene2peak.obj = "gene2peak"),
  definition = function(chip.obj, gene2peak.obj) {
    str_before <- apply(select(chip.obj@bed, chr:clus), 1, paste, collapse = "-")
    str_after <- apply(select(gene2peak.obj@bed, chr:clus), 1, paste, collapse = "-")
    idx <- order(match(str_after, str_before))
    return(new("gene2peak", bed = gene2peak.obj@bed[idx,]))
  }
)

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
    newObj <- new("gene2peak",bed=read_tsv(after,col_names=c("chr","start","end","clus","gid","dist")) %>% select(1:6))
    unlink(before)
    unlink(after)
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
    set.seed(888)
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
