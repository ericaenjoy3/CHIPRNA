#' @include CHIPRNAclass.R

#' @title create a \code{chip} object from args list
#' @name chipConst
#' @rdname chipConst-methods
#' @description create a \code{chip} object from args list
#' @param chipF bed file
#' @param sms sample name(s) for corresponding bed file(s)
#' @return a \code{chip} object
#' @export chipConst
chipConst <- function(chipF, sms) {
  if (!is.null(sms) && length(sms)>1) {
    stopifnot(length(chipF) == length(sms))
    dat.list <- lapply(seq_along(chipF), function(i){
      dat <- read_tsv(chip[i], col_names=c("chr","start","end")) %>% select(1:3) %>% mutate(clus = sms[i]);
      return(dat)
      })
    dat <- do.call(rbind,dat.list)
  } else {
    stopifnot(length(chipF) == 1)
    dat <- read_tsv(chipF, col_names = c("chr", "start", "end", "clus")) %>% select(.data$chr : .data$clus);
  }
  dat <- dat %>% mutate(clus = gsub(";", "_", .data$clus))
  return(new("chip", bed = dat))
}

#' @title Create a \code{chiploj} object from binnarized peak intersection file
#' @name chiplojConst
#' @rdname chiplojConst-methods
#' @description Create a \code{chiploj} object from binnarized peak intersection file
#' @param chipF Bed file with binnarized intersections in additional columns
#' @param reverse Reverse order of binnarized categories, when TRUE.
#' @return a \code{chiploj} object
#' @export chiplojConst
chiplojConst <- function(chipF, reverse = FALSE) {
  firstLine <- system(paste("head -1", chipF), intern = TRUE)
  dat <- read.table(chipF, header = ifelse(!grepl("^chr", firstLine), TRUE, FALSE),
    col.names = c("chr", "start", "end", "clus"), as.is = TRUE, sep="\t")
  if (ncol(dat) > 4) {
    dat[,-c(1:4)] <- apply(dat[,-c(1:4)], 2,
      function(vec){vec[vec>1] = 1;
      return(vec)
    })
    strvec <- apply(dat[,-c(1:4)], 1, paste, collapse = ";")
    sms <- colnames(dat)[-c(1:4)]
    orderdat <- expand.grid(rep(list(0:1), length(sms)))
    orderdat <- orderdat[order(rowSums(orderdat)), ]
    orderlab <- apply(orderdat, 1, function(vec)paste(sms[which(vec == 1)], collapse = ";"))
    orderstr <- apply(orderdat, 1, paste, collapse = ";")
    strvec <- factor(strvec, levels = orderstr, labels = orderlab, ordered = TRUE)
    strvec <- droplevels(strvec)
    idx <- order(as.numeric(dat[,4]), as.numeric(strvec))
    strvec <- strvec[idx]
    dat <- dat[idx,]
    fac <- paste(as.character(dat[,4]), strvec, sep = ";")
    dat[, 4] <- factor(gsub(";$", "", fac))
    return(new("chiploj", bed = as_tibble(dat[, 1:4]), binmat = as_tibble(dat[, -c(1:4)])))
  } else {
    sms <- strsplit(dat[, 4][which.max(length(dat[, 4]))], ";", fixed = TRUE)[[1]][-1]
    strvec <- gsub("^[^\\;]+;?", "", dat[, 4])
    orderdat <- expand.grid(rep(list(0:1), length(sms)))
    orderdat <- orderdat[order(rowSums(orderdat)), ]
    orderlab <- apply(orderdat, 1, function(vec)paste(sms[which(vec == 1)], collapse = ";"))
    strvec <- if (reverse) {
      factor(strvec, levels = rev(orderlab), ordered = TRUE)
    } else {
      factor(strvec, levels = orderlab, ordered = TRUE)
    }
    strvec <- droplevels(strvec)
    idx <- order(as.numeric(factor(gsub("^([^\\;]+);?.*", "\\1", dat[, 4]))),
      as.numeric(strvec))
    strvec <- strvec[idx]
    dat <- dat[idx,]
    dat[, 4] <- factor(dat[, 4])
    col.list <- strsplit(as.character(strvec), ";", fixed = TRUE)
    binmat <- data.frame(matrix(0, ncol = length(sms), nrow = nrow(dat), dimnames = list(NULL,sms)))
    invisible(sapply(seq_along(col.list), function(i){binmat[i, col.list[[i]]] <<- 1}))
    return(new("chiploj", bed = as_tibble(dat[, 1:4]), binmat = as_tibble(binmat)))
  }
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
tpmConst <- function(tpm, gene2peak.obj, small = 0.05, logit = TRUE){
  tpm.grp <- SepTPMCnt(tpm)$tpm.grp
  tpm.grp <- data.frame(gid = rownames(tpm.grp), tpm.grp) %>% mutate(gid = as.character(gid))
  # mean expression of equal-distance genes
  dat <- left_join(gene2peak.obj@bed, tpm.grp, by = c('gid' = 'gid')) %>% select(-c(5:6)) %>% data.table()
  dat <- as_tibble(dat[, lapply(.SD,mean), by = c("chr", "start", "end", "clus")])
  info <- select(dat, 1:4)
  tpm.val <- select(dat, -c(1:4))
  grp.before <- pull(dat, clus)
  # transform tpm for standardization of all samples or standardization across samples of a gene
  tpm.trans <- if (logit) {tpm.val %>% mutate_all(funs(`+`), small) %>% mutate_all(funs(log2))}
  tpm.ascale <- as_tibble(matrix(scale(as.vector(as.matrix(tpm.trans))), ncol=ncol(tpm.trans), dimnames=list(NULL, colnames(tpm.trans))))
  var.idx <- apply(tpm.trans, 1, sd) > 0
  grp.after <- grp.before[var.idx]
  tpm.rscale <- as_tibble(matrix(t(apply(tpm.trans[var.idx,], 1, scale)), ncol = ncol(tpm.trans), dimnames = list(NULL, colnames(tpm.trans))))
  return(new("tpm4plot", tpm.val = tpm.val, tpm.ascale = tpm.ascale, tpm.rscale = tpm.rscale, info = info, grp.before = grp.before, grp.after = grp.after, var.idx = var.idx))
}

limmaConst <- function(limmaF, gene2peak.obj) {
  dat <- read.table(limmaF, header = TRUE, as.is = TRUE, sep="\t")
  dat <- dat[, grep("^DEG|gid|gname|gtype", colnames(dat))]
  dat$gid <- apply(dat[, grep("^gid|gname|gtype", colnames(dat))], 1, paste, collapse = "|")
  dat <- left_join(gene2peak.obj@bed, select(dat, -c(2:3)), by = c('gid' = 'gid')) %>% select(-c(1:3))
  dat <- as.data.frame(dat)
  sub.list <- split(dat, dat$clus)
  col_names <- colnames(dat)[grep("^DEG", colnames(dat))]

  cells <- unique(c(gsub(".+_([^_\\.]+)\\.([^_\\.]+)", "\\1", col_names),
   gsub(".+_([^_\\.]+)\\.([^_\\.]+)", "\\2", col_names)))
  cells <- cells[c(grep("MEF", cells, ignore.case = TRUE),
       grep("MEF|ESC", cells, ignore.case = TRUE, invert = TRUE),
       grep("ESC", cells, ignore.case = TRUE))]
  cells <- factor(cells, levels = unique(cells),
      ordered = TRUE)

  limma.list <- lapply(col_names, function(nm, dat, sub.list){
    c1 <- gsub(".+_([^_\\.]+)\\.([^_\\.]+)", "\\1", nm)
    c2 <- gsub(".+_([^_\\.]+)\\.([^_\\.]+)", "\\2", nm)
    if (as.numeric(cells[as.character(cells) == c1]) < as.numeric(cells[as.character(cells) == c2])) {
      up_word <- "Down"
      dn_word <- "Up"
      trans <- paste(c1, "to", c2)
    } else {
      up_word <- "Up"
      dn_word <- "Down"
      trans <- paste(c2, "to", c1)
    }
    up_total <- sum(grepl(up_word, dat[, nm]), ignore.case = TRUE)
    dn_total <- sum(grepl(dn_word, dat[, nm]), ignore.case = TRUE)
    return(do.call("rbind", lapply(seq_along(sub.list), function(i) {
      subdat <- sub.list[[i]]
      up_cnt <- sum(grepl(up_word, subdat[, nm],
        ignore.case = TRUE))
      dn_cnt <- sum(grepl(dn_word, subdat[, nm],
        ignore.case = TRUE))
      up_pct <-round(100 * up_cnt / nrow(subdat), digits = 2)
      dn_pct <-round(100 * dn_cnt / nrow(subdat), digits = 2)
      return(rbind(
        data.frame(clus = names(sub.list)[i], cmp = trans, direction = "up",
          pct = up_pct, cnt = up_cnt,
          peakcnt = nrow(subdat), diffcnt = up_total),
        data.frame(clus = names(sub.list)[i], cmp = trans, direction = "down",
          pct = -dn_pct, cnt = -dn_cnt,
          peakcnt = nrow(subdat), diffcnt = dn_total)))
    })))
  }, dat = dat, sub.list = sub.list)
  limmaDF <- do.call("rbind", limma.list)
  trans <- c()
  for (i in 1:(length(cells)-1)) {
    c1 <- as.character(cells[i])
    c2 <- as.character(cells[i + 1])
    trans <- c(trans, paste(c1, "to", c2))
  }
  for (i in 3:length(cells)) {
    c1 <- as.character(cells[1])
    c2 <- as.character(cells[i])
    trans <- c(trans, paste(c1, "to", c2))
  }
  limmaDF <- limmaDF[limmaDF[, "cmp"] %in% trans, ]
  limmaDF <- limmaDF[order(match(limmaDF[, "cmp"], trans)), ]
  limmaDF[, "cmp"] <- factor(limmaDF[, "cmp"], levels = trans, ordered = TRUE)
  return(new("limma", dat = as_tibble(limmaDF)))
}
