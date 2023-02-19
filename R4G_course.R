## ----setup, include = FALSE-------------------------------------------------------------------------
options(width = 160) # to prevent premature text wrapping
knitr::opts_chunk$set(cache = FALSE, cache.path = "knitr/cache/", fig.path = "knitr/figures/")


## ----purl, eval = FALSE-----------------------------------------------------------------------------
## # file.remove("R4G_course.R") # uncomment and run to remove existing file
## knitr::purl("R4G_course.Rmd")


## ----one plus one-----------------------------------------------------------------------------------
1 + 1


## ----loading packages, message = FALSE--------------------------------------------------------------
library(adegenet)
library(pegas)
library(poppr)
library(hierfstat)
library(ggplot2)
library(lattice)
library(viridisLite)


## ----installing packages, eval = FALSE--------------------------------------------------------------
## install.packages("adegenet")


## ----sourcing---------------------------------------------------------------------------------------
source("scripts/tools.R")


## ----read genepop files-----------------------------------------------------------------------------
myData <- read.genepop(file = "data/microbov_small.gen", ncode = 3, quiet = TRUE)


## ----formating error while reading, error=TRUE------------------------------------------------------
badFormat <- read.genepop(file = "data/badFormatting.gen", ncode = 3, quiet = TRUE)


## ----input error while reading, error=TRUE----------------------------------------------------------
badInput <- read.genepop(file = "data/microbov_small.gen", ncode = 2, quiet = TRUE)


## ----listing files----------------------------------------------------------------------------------
dir(path = "data/", pattern = "*.gen")


## ----first look at a genind object------------------------------------------------------------------
myData


## ----behind the genind object 1---------------------------------------------------------------------
str(myData)


## ----behind the genind object 2, eval = FALSE-------------------------------------------------------
## print.AsIs(myData)


## ----accessor_nInd----------------------------------------------------------------------------------
nInd(myData) # Number of individuals


## ----accessor_nLoc----------------------------------------------------------------------------------
nLoc(myData) # Number of loci


## ----accessor_nPop----------------------------------------------------------------------------------
nPop(myData) # Number of populations


## ----accessor_locNames------------------------------------------------------------------------------
locNames(myData) # Names of loci


## ----accessor_alleles-------------------------------------------------------------------------------
alleles(myData) # List of all alleles


## ----accessor_nAll----------------------------------------------------------------------------------
nAll(myData, onlyObserved  = TRUE) # Number of alleles for each locus


## ----accessor_indNames------------------------------------------------------------------------------
indNames(myData) # Names of individuals


## ----accessor_popNames------------------------------------------------------------------------------
popNames(myData) # Name of the last individual in each population


## ----accessors for replacement----------------------------------------------------------------------
#let's give the pops some names:

#myPops <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5", "Pop6", "Pop7")
myPops <- paste0("Pop", 1:nPop(myData))
popNames(myData) <- myPops
popNames(myData)


## ----accessors for replacement 2--------------------------------------------------------------------
pop(myData)


## ----genind to df-----------------------------------------------------------------------------------
myData_hierf <- genind2hierfstat(myData)
head(myData_hierf) # Display first 6 rows


## ----genind to loci---------------------------------------------------------------------------------
myData_loci <- genind2loci(myData) # or as.loci(myData)
myData_loci


## ----exploring loci objects 1-----------------------------------------------------------------------
str(myData_loci)


## ----exploring loci objects 2-----------------------------------------------------------------------
head(data.frame(myData_loci), n = 10) # first 10 rows


## ----exploring loci objects 3-----------------------------------------------------------------------
attr(myData_loci, "locicol") # columns that contain loci


## ----info table-------------------------------------------------------------------------------------
missing <- info_table(myData, df = TRUE) # df = TRUE uses a long rather than 
                                         # wide format for the df
missing[order(missing$Missing, decreasing = TRUE), ] # We reorder the table to show the highest missing on top


## ----info table plot--------------------------------------------------------------------------------
info_table(myData, plot = TRUE, low = "white") ## NB: run in the console for best display


## ----ggplot, eval = FALSE---------------------------------------------------------------------------
## theplot <- ggplot_build(last_plot())
## theplot$data[[4]]$size <- 3
## theplot$data[[4]]$angle <- 45
## plot(ggplot_gtable(theplot))


## ----genotype curve---------------------------------------------------------------------------------
gencurv <- genotype_curve(myData)


## ----genotype curve in details----------------------------------------------------------------------
apply(gencurv, 2, range)


## ----remove pop by index, eval = FALSE, message = FALSE---------------------------------------------
## myData[pop = -c(7), drop = TRUE] # drop = TRUE updates the number of remaining alleles!


## ----remove pop-------------------------------------------------------------------------------------
removePop <- c("Pop7")
myData <- myData[pop = !popNames(myData) %in% removePop, drop = TRUE]


## ----remove locus-----------------------------------------------------------------------------------
removeLoc <- c("ILSTS5")
myData <- myData[loc = !locNames(myData) %in% removeLoc, drop = TRUE]


## ----new info data----------------------------------------------------------------------------------
myData


## ----new info table---------------------------------------------------------------------------------
info_table(myData, plot = TRUE, low = "white")


## ----new genotype curve-----------------------------------------------------------------------------
genotype_curve(myData)


## ----mlg--------------------------------------------------------------------------------------------
mlg(myData)


## ----mlg.id-----------------------------------------------------------------------------------------
head(mlg.id(myData))


## ----identifying mlg--------------------------------------------------------------------------------
myMLG <- mlg.id(myData)
myMLG[lengths(myMLG) > 1]


## ----dropping mlg-----------------------------------------------------------------------------------
removeInd <- c("AFBIZEB9483", "AFBIBOR9523","AFBTSOM9374")
myData <- myData[!indNames(myData) %in% removeInd, drop = TRUE]
myData


## ----dropping many loci, eval = FALSE---------------------------------------------------------------
## samplesToKeep <- unlist(lapply(myMLG, function(x) x[1])) # Capture first occurence for each MLG
## myData <- myData[indNames(myData) %in% samplesToKeep, drop = TRUE]


## ----ploidy 1---------------------------------------------------------------------------------------
head(info_table(myData, type = "ploidy"), n = 2)


## ----ploidy 2---------------------------------------------------------------------------------------
table(info_table(myData, type = "ploidy"), useNA = "always") # Counts different values


## ----ploidy 3---------------------------------------------------------------------------------------
nLoc(myData) * nInd(myData)


## ----informloci-------------------------------------------------------------------------------------
informloci(myData)


## ----informloci 2, eval = FALSE---------------------------------------------------------------------
## myData <- informloci(myData)


## ----pairwise_similarity----------------------------------------------------------------------------
pairwise_similarity(myData, pop = "Pop1")


## ----pairwise_similarity2, results = "hide"---------------------------------------------------------
lapply(popNames(myData), function(pop) hist(pairwise_similarity(myData, pop = pop, as_vector = TRUE), main = pop, xlim = c(0, 1))) # don't forget as_vector = TRUE


## ---------------------------------------------------------------------------------------------------
hist(pairwise_similarity(myData, as_vector = TRUE), xlim = c(0, 1)) # don't forget as_vector = TRUE


## ----testing HWE, results='hold'--------------------------------------------------------------------
hw.test(myData, B = 0) ## for Monte Carlo test instead of asymptotic,
                       ## use large value of B (>= 1000)


## ----HWE tested per pop-----------------------------------------------------------------------------
lapply(seppop(myData), hw.test)


## ----pop sizes--------------------------------------------------------------------------------------
unlist(lapply(seppop(myData), nInd))


## ----index of association---------------------------------------------------------------------------
ia(myData, sample = 999)


## ----index of association per pop-------------------------------------------------------------------
lapply(seppop(myData), ia, sample = 999)


## ----pair ia, results='hide'------------------------------------------------------------------------
pair.ia(seppop(myData)[["Pop6"]])
#pair.ia(myData[pop = 6])   # same thing


## ----plot pair ia, results='hide'-------------------------------------------------------------------
pair.ia(myData[pop = "Pop6"], low = "black", high = "green")


## ----plot pair ia all pop, results='hide'-----------------------------------------------------------
lapply(seppop(myData), pair.ia, low = "black", high = "green")


## ----heterozygotes, results='hold'------------------------------------------------------------------
freqP <- seq(from = 0, to = 1, by = 0.01)
freqQ <- 1 - freqP

plot(2*freqP*freqQ ~ freqQ, type = "l", las = 1,
     ylab = "frequency of heterozygote genotypes",
     xlab = "frequency of allele 'a'")
text(0, 0, "AA")
text(1, 0, "aa")
arrows(0.45, 0, 0.05, 0)
text(0.5, 0, "Aa")
arrows(0.55, 0, 0.95, 0)


## ----expected vs observed heterozygozity------------------------------------------------------------
heteroz <- data.frame(Hexp = summary(myData)$Hexp, Hobs = summary(myData)$Hobs)
heteroz$diff <- heteroz$Hexp - heteroz$Hobs
heteroz


## ----plot expected minus observed heterozygozity----------------------------------------------------
barplot(heteroz$diff, names.arg = rownames(heteroz),
        main = "Heterozygosity: expected-observed",
        xlab = "", ylab = "Hexp - Hobs", font.lab = 2, las = 2)


## ----ggplot expected minus observed heterozygozity--------------------------------------------------
heteroz$loci <- rownames(heteroz) ## ggplot needs names stored as a column

ggplot(heteroz, aes(y = Hexp, x = Hobs)) +
  geom_segment(aes(y = Hexp - 0.01, yend = Hobs, xend = Hobs), linetype = "dashed") +
  geom_text(aes(label = loci), size = 3) +
  geom_abline(slope = 1) + 
  scale_x_continuous(limits = range(c(heteroz$Hobs, heteroz$Hexp))) +
  scale_y_continuous(limits = range(c(heteroz$Hobs, heteroz$Hexp))) +
  labs(title = "Heterozygosity: expected vs observed") +
  xlab(expression(bold("Observed heterozygosity"))) +
  ylab(expression(bold("Expected heterozygosity"))) +
  theme_classic() +
  coord_fixed()


## ----expected heterozygosity using Hs---------------------------------------------------------------
Hs(myData)

## ----observed heterozygosity using Ho---------------------------------------------------------------
Ho(myData)


## ----expected and observed heterozygosity per pop---------------------------------------------------
heteroz_per_pop <- data.frame(Hexp = Hs(myData), Hobs = Ho(myData))
heteroz_per_pop$diff <- heteroz_per_pop$Hexp - heteroz_per_pop$Hobs
heteroz_per_pop


## ----heteroz matrix---------------------------------------------------------------------------------
heteroz_matrix <- do.call("cbind", lapply(seppop(myData),
                                          function(x) summary(x)$Hexp - summary(x)$Hobs))
heteroz_matrix


## ----plot heteroz matrix----------------------------------------------------------------------------
levelplot(heteroz_matrix, col.regions = viridis(100),
          main = "Heterozygosity: expected-observed",
          xlab = "Locus", ylab = "Population",
          scales = list(x = list(rot = 90)))


## ----citation hierfstat-----------------------------------------------------------------------------
citation(package = "hierfstat")
packageVersion("hierfstat") ## don't forget to cite the version number!


## ----basic stats hierf------------------------------------------------------------------------------
basic.stats(myData)


## ----basic stats Fis--------------------------------------------------------------------------------
basic.stats(myData)$Fis # Fis for every locus (rows) by population (columns)


## ----basic stats Fis2-------------------------------------------------------------------------------
Fis <- colMeans(basic.stats(myData)$Fis) # close but not equal to (Hs(myData) - Ho(myData))/Hs(myData)
Fis


## ----Fis boot---------------------------------------------------------------------------------------
Fis_CI <- boot.ppfis(myData, nboot = 999)
Fis_CI


## ----Fis all----------------------------------------------------------------------------------------
Fis_table <- data.frame(Fis = Fis, Fis_CI$fis.ci)
colnames(Fis_table)[2:3] <- c("lwr", "upr")
Fis_table$lwr[Fis_table$lwr < 0] <- 0 # lwr than 0 does not make sense (bootstrap artifact)
Fis_table


## ----Fst per pair-----------------------------------------------------------------------------------
pairFst <- pairwise.WCfst(genind2hierfstat(myData)) # does not work directly on genind input :-(
pairFst


## ----run pairwise Fst with uncertainty plot---------------------------------------------------------
levelplot(pairFst, col.regions = rev(grey.colors(30)))


## ----bootstrapping Fst------------------------------------------------------------------------------
boot.ppfst(genind2hierfstat(myData), nboot = 999)


## ----private alleles--------------------------------------------------------------------------------
private <- private_alleles(myData, report = "data.frame") ## report influences the output format
private


## ----private alleles2, message=FALSE----------------------------------------------------------------
private_only <- private[private$count > 0, ]
private_only[order(private_only$population), ]


## ----private alleles3, message=FALSE, error=TRUE----------------------------------------------------
library(dplyr)
private |>
  filter(count > 0) |>
  arrange(population) |>
  rename(nb_indiv = count)


## ----private alleles per pop------------------------------------------------------------------------
private_alleles_per_pop <- rowSums(private_alleles(myData) > 0)
private_alleles_per_pop


## ----nb of alleles vs pop size, results='hold'------------------------------------------------------
ind_per_pop <- sapply(seppop(myData), nInd)
plot(ind_per_pop, private_alleles_per_pop,
     xlab = "Sample size", ylab = "Number of private alleles", las = 1,
     col = NULL)
text(y = private_alleles_per_pop,
     x = ind_per_pop,
     labels = names(ind_per_pop))


## ----proportion of shared alleles-------------------------------------------------------------------
similarity_mat <- propShared(myData[pop = "Pop1"])
similarity_mat[1:5, 1:5]


## ----neighbour joining------------------------------------------------------------------------------
distance_mat <- 1 - similarity_mat
mynj <- nj(distance_mat)
mynj


## ----plot unrooted----------------------------------------------------------------------------------
plot(mynj, type = "unrooted", cex = 0.8)


## ----plot prevosti small----------------------------------------------------------------------------
mynj_prevosti <- nj(prevosti.dist(seppop(myData)[["Pop1"]]))
plot(mynj_prevosti, type = "unrooted", cex = 0.8)


## ----co-phylogenies, fig.width=10-------------------------------------------------------------------
cophyloplot(mynj, mynj_prevosti,
            assoc =  cbind(mynj$tip.label, mynj_prevosti$tip.label), cex = 0.8)


## ----plot prevosti big,  fig.width=10---------------------------------------------------------------
bignj <- nj(prevosti.dist(myData))
plot(bignj, type = "unrooted", cex = 0.8)


## ----plot nj better 1, fig.height = 9, fig.width = 9------------------------------------------------
plot(bignj, type = "fan", show.tip.label = FALSE, x.lim = c(-0.7, 0.7), no.margin = TRUE)
tiplabels(text = rownames(myData@tab),
          frame = "none",
          col = rainbow(nPop(myData))[as.numeric(myData@pop)], cex = 0.8, offset = 0.05)
legend("topleft", fill = rainbow(nPop(myData)),
       legend = popNames(myData), bty = "n",
       title = "Population")


## ----plot nj better 2, fig.height = 12, fig.width = 12----------------------------------------------
plot(bignj, type = "radial", show.tip.label = FALSE, x.lim = c(-0.7, 0.7), no.margin = TRUE)
tiplabels(text = rownames(myData@tab),
          frame = "none",
          col = rainbow(nPop(myData))[as.numeric(myData@pop)], cex = 0.8, offset = 0.05)
legend("topleft", fill = rainbow(nPop(myData)),
       legend = popNames(myData), bty = "n",
       title = "Population")


## ----pca hierf--------------------------------------------------------------------------------------
plot(indpca(myData))


## ----pca--------------------------------------------------------------------------------------------
myData_matrix <- scaleGen(myData, center = FALSE, scale = FALSE, NA.method = "mean")
mypca <- dudi.pca(myData_matrix, center = TRUE, scale = FALSE, scannf = FALSE, nf = Inf)


## ----pca explained 1--------------------------------------------------------------------------------
head(mypca$li[, 1:4]) ## only show head for first 4 axes


## ----pca explained 2--------------------------------------------------------------------------------
head(zapsmall(cor(mypca$li))[, 1:4])  ## only show head for first 4 axes


## ----pca explained 3--------------------------------------------------------------------------------
barplot(mypca$eig,
        names.arg = colnames(mypca$li),
        cex.names = 0.5,
        col = heat.colors(length(mypca$eig)),
        las = 2, ylab = "Inertia")


## ----pca explained 4--------------------------------------------------------------------------------
barplot(cumsum(100*mypca$eig/sum(mypca$eig)),
        names.arg = colnames(mypca$li),
        cex.names = 0.5,
        col = rev(heat.colors(length(mypca$eig))),
        las = 2, log = "y",
        ylab = "Cumulative proportion of variance explained")


## ----plot PCA 1-------------------------------------------------------------------------------------
s.class(mypca$li, fac = pop(myData),
        col = rainbow(nPop(myData)), grid = FALSE, xax = 1, yax = 2, cpoint = 0)
s.label(mypca$li, add.plot = TRUE, boxes = FALSE, clabel = 0.5)
add.scatter.eig(mypca$eig[1:10], xax = 1, yax = 2, ratio = 0.15)


## ----contribution to dimension----------------------------------------------------------------------
s.arrow(mypca$co[, 1:2], boxes = FALSE, cpoint = 0, clabel = 0)
private_alleles_pop6 <- colnames(private_alleles(myData))[private_alleles(myData)["Pop6", ] != 0]
s.arrow(mypca$co[private_alleles_pop6, 1:2], boxes = FALSE, clabel = 0.8, add.plot = TRUE)


## ----trimming missing data--------------------------------------------------------------------------
myData_noNA <- missingno(myData, type = "geno", cutoff = 0) # xvalDapc cannot handle missing values so we remove them


## ----clusters---------------------------------------------------------------------------------------
grp <- find.clusters(myData_noNA, n.pca = Inf, n.clust = 4) ## best not to define n.clust and choose interactivelly, keep n.pca = Inf here


## ----dapc-------------------------------------------------------------------------------------------
dapc1 <- dapc(myData_noNA, pop = grp$grp, n.pca = 10, n.da = Inf) ## for choice of n.pca see below (it is very important!) 
scatter(dapc1)


## ----dapc plot comparing outcome to pop 2-----------------------------------------------------------
compoplot(dapc1, show.lab = TRUE, legend = FALSE, cex.names = 0.4,
          lab = paste(pop(myData_noNA), rownames(dapc1$tab)))


## ----new id-----------------------------------------------------------------------------------------
newID <- data.frame(INRA63.183 = 2, INRA63.181 = 0, INRA63.177 = 0, 
                    INRA63.175 = 0, INRA63.179 = 0, INRA63.185 = 0, INRA63.171 = 1, 
                    INRA5.137 = 0, INRA5.141 = 2, INRA5.143 = 0, INRA5.139 = 0, 
                    INRA5.149 = 0, ETH225.147 = 0, ETH225.157 = 1, ETH225.139 = 0, 
                    ETH225.141 = 0, ETH225.153 = 0, ETH225.149 = 0, ETH225.155 = 0, 
                    ETH225.143 = 0, ETH225.159 = 1, ETH225.151 = 0, ETH225.137 = 0, 
                    ETH225.145 = 0, HE5.149 = 0, HE5.151 = 0, HE5.163 = 0, 
                    HE5.165 = 0, HE5.167 = 2, HE5.155 = 0, HE5.157 = 0, 
                    HE5.153 = 0, HE5.161 = 0, HE5.159 = 0, HE1.103 = 1, 
                    HE1.109 = 0, HE1.105 = 1, HE1.107 = 0, HE1.101 = 0, 
                    HE1.117 = 0, HE1.113 = 0, HE1.111 = 0, INRA35.102 = 0, 
                    INRA35.104 = 2, INRA35.120 = 2, INRA35.110 = 0, INRA35.108 = 0, 
                    INRA35.114 = 0, INRA35.106 = 0, ETH152.191 = 0, ETH152.195 = 1, 
                    ETH152.197 = 1, ETH152.193 = 0, ETH152.201 = 0, ETH152.199 = 0, 
                    ETH152.205 = 0, INRA23.213 = 0, INRA23.215 = 0, INRA23.197 = 0, 
                    INRA23.199 = 1, INRA23.209 = 2, INRA23.205 = 0, INRA23.211 = 1, 
                    INRA23.203 = 0, INRA23.207 = 0, INRA23.217 = 0, INRA23.201 = 0, 
                    INRA23.193 = 0, ETH10.209 = 0, ETH10.211 = 0, ETH10.207 = 1, 
                    ETH10.217 = 0, ETH10.219 = 1, ETH10.221 = 0, ETH10.215 = 0, 
                    ETH10.223 = 0, ETH10.213 = 2, HE9.153 = 0, HE9.159 = 0, 
                    HE9.165 = 1, HE9.155 = 0, HE9.161 = 1, HE9.163 = 0, 
                    HE9.167 = 0, HE9.157 = 0, HE9.149 = 0, HE9.169 = 0, 
                    HE9.151 = 0, CSSM66.181 = 0, CSSM66.189 = 1, CSSM66.185 = 1, 
                    CSSM66.193 = 0, CSSM66.197 = 0, CSSM66.183 = 0, CSSM66.199 = 0, 
                    CSSM66.187 = 0, CSSM66.179 = 0, CSSM66.195 = 0, CSSM66.171 = 0, 
                    CSSM66.191 = 0, INRA32.162 = 1, INRA32.178 = 0, INRA32.180 = 0, 
                    INRA32.184 = 0, INRA32.202 = 0, INRA32.182 = 1, INRA32.176 = 0, 
                    INRA32.174 = 0, INRA32.164 = 0, INRA32.160 = 0, INRA32.204 = 0, 
                    INRA32.168 = 0, ETH3.117 = 0, ETH3.125 = 0, ETH3.115 = 2, 
                    ETH3.127 = 0, ETH3.103 = 0, ETH3.123 = 0, ETH3.129 = 0, 
                    ETH3.131 = 0, ETH3.119 = 0, ETH3.109 = 0, ETH3.113 = 0, 
                    BM2113.140 = 0, BM2113.142 = 0, BM2113.122 = 0, BM2113.134 = 0, 
                    BM2113.136 = 1, BM2113.146 = 1, BM2113.138 = 0, BM2113.130 = 0, 
                    BM2113.124 = 0, BM2113.150 = 0, BM2113.128 = 0, BM2113.132 = 0, 
                    BM2113.126 = 0, BM2113.144 = 0)
predict(dapc1, newdata = newID)


## ----clusters vs pop--------------------------------------------------------------------------------
table.value(table(pop(myData_noNA), grp$grp), col.lab = paste0("cluster_", 1:4))


## ----cross validation pca---------------------------------------------------------------------------
cross_validation <- xvalDapc(myData_noNA, grp = grp$grp, n.da = Inf, n.rep = 30) # increase n.rep for real usage!
cross_validation$`Number of PCs Achieving Highest Mean Success`

