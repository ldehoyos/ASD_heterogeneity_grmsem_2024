##################################################-
## Project: autism grmsem 2024
## Description: predict number of genetic factors
## Date: September 2023
## Author: Lucia de Hoyos
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())


# ------------------------------------------- #
# 01. Load libraries and functions ----------
# ------------------------------------------- #

library(nFactors)
library(ggplot2)

# ------------------------------------------- #
# 02. Prepare data --------------------------
# ------------------------------------------- #

# Load Cholesky
fname_cholesky <- "Cholesky_SPARK_model.RData"

myEnv <- new.env()
load(file = fname_cholesky, envir = myEnv)

# Names of the phenotypes
ph.nms <- myEnv$out$ph.nms

# Correlation matrix
corrMat <- myEnv$var.out$RG   
colnames(corrMat) <- rownames(corrMat) <- ph.nms


# ------------------------------------------- #
# 03. PCA check -----------------------------
# ------------------------------------------- #

# Run analysis of the number of factors to retain 

# Obtain eigenvalues and eigenvectors
ev <- eigen(corrMat)

# Run analysis of the number of factors to retain 
nS <- nScree(x=ev$values)

# Get vector with eigenvalues and positions
eig <- ev$values
k   <- 1:length(eig) 

# Run linear model
vp.p <- lm(eig[c(nS$Components$noc + 1, length(ev$values))] ~ k[c(nS$Components$noc + 1, length(ev$values))]) 

# get slope and intercept (y= m*x + n)
sl     <- coef(vp.p)[2]
intrct <- coef(vp.p)[1]

cat("\nRetain the number of factors based on Kaiser (eigenvalues > 1) and optimal coordinate criterion\n", 
    "\nKaiser criterion: \t\t",nS$Components$nkaiser, "\n",
    "\nOptimal coordiate criterion: \t",nS$Components$noc, "\n")


var_explained_df <- data.frame(C= paste0(1:length(ev$values)),
                               Eigenvalue=ev$values,
                               KC=c(rep("Y", nS$Components$nkaiser),rep("N", length(ph.nms)-nS$Components$nkaiser)),
                               OC=c(rep("Y", nS$Components$noc),rep("N", length(ph.nms)-nS$Components$noc)))

# Plot
ggplot(var_explained_df, aes(x=C,y=Eigenvalue)) +
  geom_abline(slope = sl, intercept = intrct, size=0.75, linetype="dashed") +
  geom_hline(yintercept = 1, alpha=0.2) +
  geom_line(aes(group=1),size=0.75) + xlab("Number of factors") +
  geom_point(aes(color=KC), size = 4) + scale_color_manual(values = c("grey","black")) +
  theme_classic()  
