#' @include CHIPRNAclass.R

#' @title create a \code{chip} object from args list
#' @name chipConst
#' @rdname chipConst-methods
#' @description create a \code{chip} object from args list
#' @param chip bed file(s)
#' @param sms sample name(s) for corresponding bed file(s)
#' @return a \code{chip} object
#' @export chipConst
chipConst <- function(chip, sms) {
  if(!is.null(sms) && length(sms)>1) {
    stopifnot(length(args$chip)==length(args$sms))
    dat.list<-lapply(seq_along(chip),function(i){
      dat<-read_tsv(chip[i],col_names=c("chr","start","end")) %>% select(1:3) %>% mutate(clus = sms[i]);
      return(dat)
      })
    dat<-do.call(rbind,dat.list)
  } else {
    stopifnot(length(chip)==1)
    dat<-read_tsv(chip, col_names=c("chr", "start", "end", "clus")) %>% select(.data$chr : .data$clus);
  }
  dat <- dat %>% mutate(clus = gsub(";", "_", clus))
  return(new("chip", bed = dat))
}

#' @title TPM contrast
#' @name tpmConst
#' @rdname tpmConst-methods
#' @description Do something
#' @param tpm a tpm object
#' @param gene2peak a gene2peak class
#' @param small a numeric value for addding to TPM value
#' @param logit a logical value for whethether to log2 transform TPM value before proceeding to tpm.ascale and tpm.rscale
#' @return A \code{tpm4plot} object
#' @export tpmConst
tpmConst<-function(tpm, gene2peak, small = 0.05, logit = TRUE){
  if (!is.null(tpm)) {
    tpm.grp <- SepTPMCnt(tpm)$tpm.grp
    tpm.grp <- data.frame(gid = rownames(tpm.grp), tpm.grp) %>% mutate(gid = as.character(gid))
    # mean expression of equal-distance genes
    dat <- left_join(gene2peak@bed, tpm.grp, by = c('gid' = 'gid')) %>% group_by(.data$chr, .data$start, .data$end) %>% mutate_at(vars(-(.data$clus:.data$dist)),mean) %>% ungroup()
    info <- select(dat,1:6)
    tpm.val <- select(dat,-c(1:6))
    grp.before <- pull(gene2peak@bed, .data$clus)
    # transform tpm for standardization of all samples or standardization across samples of a gene
    tpm.trans <- if (logit) {tpm.val %>% mutate_all(funs(`+`), small) %>% mutate_all(funs(log2))}
    tpm.ascale <- as_tibble(matrix(scale(as.vector(as.matrix(tpm.trans))), ncol=ncol(tpm.trans), dimnames=list(NULL, colnames(tpm.trans))))
    var.idx <- apply(tpm.trans, 1, sd) > 0
    grp.after <- grp.before[var.idx]
    tpm.rscale <- as_tibble(matrix(t(apply(tpm.trans[var.idx,], 1, scale)), ncol = ncol(tpm.trans), dimnames = list(NULL, colnames(tpm.trans))))
    return(new("tpm4plot", tpm.val = tpm.val, tpm.ascale = tpm.ascale, tpm.rscale = tpm.rscale, info = info, grp.before = grp.before, grp.after = grp.after, var.idx = var.idx))
  }
}
