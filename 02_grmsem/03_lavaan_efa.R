##################################################-
## Project: autism grmsem 2024
## Description: predict number of genetic factors
## Date: September 2023
## Author: Lucia de Hoyos, Maria T (Mariska) Barendse, Beate St Pourcain
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())


# ------------------------------------------- #
# 01. Load libraries and functions ---------
# ------------------------------------------- #

library(Matrix)
library(lavaan)

# ------------------------------------------- #
# 02. Prepare data --------------------------
# ------------------------------------------- #

# Load Cholesky
fname_cholesky <- "Cholesky_SPARK_model.RData"

myEnv <- new.env()
load(file = fname_cholesky, envir = myEnv)

# Names of the phenotypes
ph.nms <- myEnv$out$ph.nms

# Covariance matrix
covMat <- round(myEnv$var.out$VA, digits = 6)
colnames(covMat) <- rownames(covMat) <- ph.nms

# smooth covariance matrix
covMat.smooth <- as.matrix((nearPD(covMat))$mat)

# weight matrix 
wMat.simple <- 1/((myEnv$var.out$VA.se)^2) 
wMat        <- diag(diag(wMat.simple))
if(length(wMat)==0){wMat <- diag(x = 1, nrow = length(ph.nms))}
colnames(wMat) <- rownames(wMat) <- ph.nms

# ------------------------------------------- #
# 03. Run EFA -------------------------------
# ------------------------------------------- #

numFacts <- 3 # from previous script

# Define model
modelEFA <- paste0(paste0("efa('efa')*f",1:numFacts, collapse = " + "), " =~ ",
                   paste0(ph.nms, collapse = " + "), "\n",  paste0("f",1:numFacts, "~~1*f", 1:numFacts, collapse = "\n"))
cat(modelEFA)


# Fit EFA model
fit.modelEFA <- sem(model = modelEFA, 
                    sample.cov = covMat.smooth, sample.nobs = myEnv$out$n.ind,
                    estimator="DWLS", WLS.V = wMat, 
                    rotation="varimax", rotation.args = list(orthogonal = TRUE))
summary(fit.modelEFA)


# ------------------------------------------- #
# 04. Define factor structure for grmsem ----
# ------------------------------------------- #

# Common genetic factors
lavaanFit.AC  <- lavInspect(fit.modelEFA, what = "est")$lambda
colnames(lavaanFit.AC) <- paste0("F",1:numFacts)

# Specific genetic factors
lavaanFit.AS    <- suppressWarnings(data.frame(AS=sqrt(diag(lavInspect(fit.modelEFA, what = "est")$theta))))
lavaanFit.AS$AS <- ifelse(is.na(lavaanFit.AS$AS),0,lavaanFit.AS$AS)

# Merge common and specific genetic factor loadings
EFA.grmsem <- cbind(lavaanFit.AC,lavaanFit.AS)

# ------------------------------------------- #
# 05. Get starting values for grmsem --------
# ------------------------------------------- #

#### Get starting values, genetic part
aVal <- unlist(EFA.grmsem, use.names = FALSE)
paste0("A.v.startval=c(",paste0(round(aVal,3),collapse = ","),")")

#### Get starting values, residual part
eVal <- myEnv$out$model.out
paste0("E.v.startval=c(",paste0(round(eVal$estimates[grep("^e", eVal$label)],4),collapse = ","),")")

# ------------------------------------------- #
# 06. Define free parameters for IPC --------
# ------------------------------------------- #

# For common genetic factors, retain factor loadings above threshold, 
#         in this 0.1 (1% phenotypic variation)
freepar.AC <- ifelse(c(abs(lavaanFit.AC))>0.1,1,0)

# Specific factors are always estimated, never constrained
freepar.AS <- rep(1, times = nrow(lavaanFit.AS))

# Define free parameters, if 0 it is not estimated
Av.freepar <- c(freepar.AC, freepar.AS)
paste0("A.v.free=c(",paste0(Av.freepar,collapse = ","),")")


# ------------------------------------------- #
# 07. Define free parameters for bifactor ----
# ------------------------------------------- #

# Set the parameters of the first factor to estimated
Av.freepar.bifactor <- replace(Av.freepar, 1:length(ph.nms), 1)
paste0("A.v.free=c(",paste0(Av.freepar.bifactor,collapse = ","),")")


