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
  if (requireNamespace("progress", quietly = TRUE) & interactive()) {
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


# TODO: this is ugly and diploid only but I failed to find something smart...
sharedness <- function(id1, id2) {
  a <- id1[1]
  b <- id1[2]
  c <- id2[1]
  d <- id2[2]
  l <- length(unique(c(a, b, c, d)))
  l1 <- length(unique(id1))
  l2 <- length(unique(id2))
  if (l == 1) return(2)
  if (l == 2 & (l1 == 1 & l2 == 1)) return(0)
  if (l == 2 & (l1 != 1 & l2 != 1)) return(2)
  if (l == 2 & ((l1 == 1) | (l2 == 1))) return(1)
  if (l == 3 & ((l1 == 1) | (l2 == 1))) return(0)
  if (l == 3 & (l1 != 1 & l2 != 1)) return(1)
  if (l == 4) return(0)
  print(paste("id1 = ", id1, collapse = " "))
  print(paste("id2 = ", id2, collapse = " "))
  stop("case not accounted for in sharedness -> bug")
  }


similarity <- function(x, id1, id2) {
  
  idAB <- x@tab[c(id1, id2), ]
  cols_with_NA <- apply(apply(idAB, 1, is.na), 1, any)
  idAB <- idAB[, !cols_with_NA]
  
  alleles_count_to_df <- function(allele_counts_row) {
    alleles <- table(rep(names(allele_counts_row), allele_counts_row))
    alleles <- rep(names(alleles), alleles)
    alleles_list <- strsplit(alleles, ".", fixed = TRUE)
    data.frame(locus = unlist(lapply(alleles_list, function(x) x[[1]])),
               allele = unlist(lapply(alleles_list, function(x) x[[2]])))
  }

  allelesA <- alleles_count_to_df(idAB[1, ])
  allelesB <- alleles_count_to_df(idAB[2, ])
  allelesAB <- data.frame(locus = allelesA$locus,
                          idA =  allelesA$allele,
                          idB = allelesB$allele)
  
  common <- sapply(by(allelesAB[, c("idA", "idB")],
                                  INDICES = allelesAB$locus,
                                  function(x) sharedness(x[, 1], x[, 2])),
                   "c")
  
  allelesAB$common <- rep(common, 2) ## for debugging
  
  comparison_summary <- data.frame(id1 = id1, id2 = id2,
                                   common = sum(common),
                                   distinct = length(common) * 2 - sum(common),
                                   total = length(common) * 2,
                                   ratio = sum(common)/(length(common) * 2),
                                   stringsAsFactors = FALSE)
  comparison_summary
}
# data.frame(as.loci(myData[pop = 1, ]))[c(1,3), ]
# similarity(myData, id1 = "AFBIBOR9503", id2 = "AFBIBOR9505")
# propShared(myData[pop = 1])["AFBIBOR9503", "AFBIBOR9505"] 



pairwise_similarity <- function(x, pop = NULL, as_vector = FALSE) {
  if ("loci" %in% class(x)) {
    x <- pegas::loci2genind(x)
  }
  if (!"genind" %in% class(x)) stop("x must be a genind object!")
  
  if (!is.null(pop)) x <- x[pop = pop]
  
  if (as_vector) {
    ratio <- adegenet::propShared(x) # faster than using similarity()
    return(ratio[upper.tri(ratio)])
  }
  
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
  if (requireNamespace("progress", quietly = TRUE) & interactive()) {
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

print("All functions succesfully loaded")
