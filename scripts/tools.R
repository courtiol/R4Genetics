#Note: this was vaguely inspired from the code in https://github.com/jgx65/hierfstat/blob/master/R/pairwise.fst.R.
#TODO: compare to the bootstrap procedure in hierfstat and perhaps do the same.  https://github.com/jgx65/hierfstat/blob/master/R/ppboot.R

pairwise_F <- function(x, stat = "Fst", confint = TRUE, nboot = 100) {
  
  confint_fn <- function(x) {
    boot_res <- replicate(nboot, {
      x_boot <- x[, c(1, sample(attr(x, "loci"), replace = TRUE))]
      mean(pegas::Fst(x_boot)[, stat])
    })
    return(stats::quantile(boot_res, c(0.025, 0.975)))
  }
  
  if ("loci" %in% class(x)) {
    x <- pegas::loci2genind(x)
  }
  if (!"genind" %in% class(x)) stop("x must be a genind object!")
  x <- adegenet::seppop(x)
  pops <- names(x)
  npop <- length(pops)
  if (npop < 2) stop("x needs to contain at least 2 populations!")
  fmat <- matrix(nrow = npop, ncol = npop, dimnames = list(pops, pops))
  fmat_upr <- fmat_lwr <- fmat
  pb <- NULL
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(total = npop*(npop - 1)/2)
  }
  if (!isNamespaceLoaded("progress")) {
    message("install the package progress if you want to see a progress bar.")
  }
  for (a in 2:npop) {
    for (b in 1:(a - 1)) {
      if (!is.null(pb)) pb$tick()
      sub_x <- pegas::as.loci(adegenet::repool(x[pop = c(pops[a], pops[b])]))
      fmat[a, b] <- fmat[b, a] <- mean(pegas::Fst(sub_x)[, stat])
      if (confint) {
        CI <- confint_fn(sub_x)
        fmat_lwr[a, b] <- fmat_lwr[b, a] <- CI[1]
        fmat_upr[a, b] <- fmat_upr[b, a] <- CI[2]
      }
    }
  }
  if (confint) return(list(mean = fmat, upr = fmat_upr, lwr = fmat_lwr))
  return(fmat)
}

similarity <- function(x, id1, id2) {
  idA <- x@tab[id1, ]
  idB <- x@tab[id2, ]
  idA_clean <- names(idA)[!is.na(idA) & idA != 0]
  idB_clean <- names(idB)[!is.na(idB) & idB != 0]
  comp <- c(distinct = length(setdiff(idA_clean, idB_clean)), total = length(unique(c(idA_clean, idB_clean))))
  comp <- c(comp, common = comp[["total"]] - comp[["distinct"]])
  ratio <- comp[[3]] / comp[[2]]
  data.frame(id1 = id1, id2 = id2, distinct = comp[[1]],
             common = comp[[3]], total = comp[[2]], ratio = ratio,
             stringsAsFactors = FALSE)
}

#similarity(myData, id1 = "N7", id2 = "N145")

##TODO: optimise and add pop
pairwise_similarity <- function(x) {
  if ("loci" %in% class(x)) {
    x <- pegas::loci2genind(x)
  }
  if (!"genind" %in% class(x)) stop("x must be a genind object!")
  indivs <- rownames(x@tab)
  nindivs <- length(indivs)
  if (nindivs < 2) stop("x needs to contain at least 2 individuals!")
  ncases <- nindivs*(nindivs - 1) / 2
  res <- data.frame(id1 = rep(NA_character_, ncases),
                    id2 = rep(NA_character_, ncases),
                    distinct = rep(NA_real_, ncases),
                    common = rep(NA_real_, ncases),
                    total = rep(NA_real_, ncases),
                    ratio = rep(NA_real_, ncases),
                    stringsAsFactors = FALSE)
  pb <- NULL
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(total = ncases)
  }
  if (!isNamespaceLoaded("progress")) {
    message("install the package progress if you want to see a progress bar.")
  }
  i <- 1
  for (i_id1 in 2:nindivs) {
    for (i_id2 in 1:(i_id1 - 1)) {
      if (!is.null(pb)) pb$tick()
      res[i, ] <- similarity(x, indivs[i_id1], indivs[i_id2])
      i <- i + 1
    }
  }
  return(res)
}

#pairwise_similarity(myData[pop = 1])
#lapply(seppop(myData), pairwise_similarity)
