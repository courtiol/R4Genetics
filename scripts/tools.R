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

## TODO: FIX BUG? -> different from outcomes of propShared() but not sure why...
similarity <- function(x, id1, id2) {
  idA <- x@tab[id1, ]
  idB <- x@tab[id2, ]
  allelesA <- table(rep(names(idA), idA))
  allelesA <- rep(names(allelesA), allelesA)
  allelesB <- table(rep(names(idB), idB))
  allelesB <- rep(names(allelesB), allelesB)
  comparison <- do.call("rbind", lapply(locNames(x), function(locus) {
    A <- unique(allelesA[grepl(locus, allelesA)])
    B <- unique(allelesB[grepl(locus, allelesB)])
    total <- length(unique(c(A, B)))
    distinct <- length(c(setdiff(A, B), setdiff(B, A)))
    common <- total - length(c(setdiff(A, B), setdiff(B, A)))
    return(data.frame(locus, total, distinct, common))
  }))
  comparison_summary <- apply(comparison[, c("total", "distinct", "common")], 2, sum)
  comparison_summary <- data.frame(id1 = id1, id2 = id2,
                                   distinct = comparison_summary["distinct"],
                                   common = comparison_summary["common"],
                                   total = comparison_summary["total"],
                                   ratio = comparison_summary["common"]/comparison_summary["total"],
                                   stringsAsFactors = FALSE)
  comparison_summary
}

#similarity(myData, id1 = "N7", id2 = "N142")
#propShared(myData[pop = 1])["N7", "N142"]
#data.frame(as.loci(myData[pop = 1, ]))[c(1,3), ]

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

#test <- pairwise_similarity(myData[pop = 1])
#res <- cbind(test, ratio2 = propShared(myData[pop = 1])[cbind(test[, 1], test[, 2])])
#plot(res$ratio, res$ratio2)
