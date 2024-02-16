##################################################-
## Project: autism grmsem 2024
## Description: prepare grm for grmsem analysis
## Date: September 2023
## Author: Beate St Pourcain
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())

# ------------------------------------------- #
# Functions ----------------------------------
# ------------------------------------------- #
grm.input <- function(grm){
  arg <- sqrt(1 + 8*nrow(grm))
  n <- (-1 + arg)/2
  A <- matrix(0, n, n)
  A[upper.tri(A, diag=T)] <- grm[,4]
  A <- A + t(A)
  diag(A) <- .5*diag(A)
  return(A)
}

# ------------------------------------------- #
# 01. Format GRM ----------------------------
# ------------------------------------------- #

# GRM was computed using the following code:
# spark_probands_ibd125 is the cleaned PLINK file after quality control (see Supplementary Methods 1)
# system("./plink --bfile spark_probands_ibd125 --rel-cutoff 0.05 --make-grm-gz --out spark_grm")

# Read in GRM (tab delimited)
grm <- read.table(file = "spark_grm.grm", sep="\t")

# Reformat
G <- grm.input(grm)

rm(grm, grm.input)

# ------------------------------------------- #
# 02. Format GRM ids ------------------------
# ------------------------------------------- #

# Read in ID file (tab delimited)
grm.id <- read.table(file = "spark_grm.grm.id", sep="\t")
colnames(grm.id) <- c("FID", "IID")

# ------------------------------------------- #
# 03. Load data -----------------------------
# ------------------------------------------- #

# Read in your phenotype file
pheno <- as.data.frame(read.table(file = paste0(ph.dir,ph.file), header = TRUE))

# Create an index based on grm
grm.id$seq <- seq(1:nrow(grm.id))

# Merge your phenotype with the grm based on the index
ph.srt <- merge(x = grm.id, y = pheno, by = c("FID","IID"), all.x = TRUE)
ph.srt <- ph.srt[order(ph.srt$seq),]
summary(ph.srt)

rm(pheno, na_string, ph.dir, separator, w.dir, grm.RData, ph.file)
save.image("grm.spark_grm_pheno.RData")
