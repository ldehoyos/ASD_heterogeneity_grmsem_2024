##################################################-
## Project: autism grmsem 2024
## Description: run phenotype SEM in lavaan
## Date: September 2023
## Author: Lucia de Hoyos & Beate St Pourcain
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())


# ------------------------------------------- #
# 01. Load libraries and functions ----------
# ------------------------------------------- #

library(nFactors)
library(psych)
library(lavaan)
library(dplyr)
library(tidyr)

# ------------------------------------------- #
# 02. Prepare data --------------------------
# ------------------------------------------- #

# Load phenotype data
pheno <- as.data.frame(read.table(file = paste0(ph.dir,ph.file), header = TRUE))

# Load grmsem best-fitting model output
grmsem.model <- "IPC_3factor_SPARK_model.RData"
grmsem.env   <- new.env()
load(grmsem.model, envir = grmsem.env)

# Get phenotypes of interest
ph.nms   <- data.frame(ph=grmsem.env$out$ph.nms,stringsAsFactors = F)
ph.nms$names <- c("lang_dis","lang_lvl","bhv_odd","age_crawl","selffeed","contr_mv","selfinj","sameness")
print(ph.nms$ph)

# ------------------------------------------- #
# 03. Run data checks -----------------------
# ------------------------------------------- #

# Get information on missing data
pheno %>% select(ph.nms$ph) %>% ff_glimpse() %>% .[[1]] 

# Get missing plot
pheno %>% select(ph.nms$ph) %>% missing_plot()


# Get information on which phenotypes are NA
pheno.miss <- pheno %>% select(ph.nms$ph) %>% is.na()
pheno.misvar <- data.frame(FID= pheno$FID, IID= pheno$IID, nmis= rowSums(pheno.miss), 
                           misspat= apply(pheno.miss, 1, function(x) paste(as.numeric(x), collapse = ""))) 

# Subset your phenotype data  for individuals that have data for at least 4 phenotypes
pheno <- pheno.misvar %>% filter(nmis<=4) %>% inner_join(y = pheno, by = c("FID","IID"))
cat("\nAfter filtering for missing data\n", nrow(pheno), "invididuals remain")

# ------------------------------------------- #
# 04. PCA check -----------------------------
# ------------------------------------------- #
# Run analysis of the number of factors to retain 

# Select variables
pheno.sel <- pheno[,ph.nms$ph]

## Compute correlation matrix
corrMat <- cor(x = pheno.sel, use = "pairwise.complete.obs", method = "pearson")

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

# Number of factors
numFacts <- nS$Components$noc

# ------------------------------------------- #
# 05. Divide sample into two halves ---------
# ------------------------------------------- #
# Match on missingness and sex

# Get EFA dataset
pheno.efa     <- pheno %>% group_by(nmis,Sex) %>% slice_sample(prop=0.5) %>% ungroup() 
pheno.efa.sel <- pheno.efa %>% select(ph.nms$ph)
pheno.efa.sel  %>% is.na() %>% colSums()

# Get CFA dataset
pheno.cfa     <- pheno %>% anti_join(pheno.efa, by=c("FID","IID"))
pheno.cfa.sel <- pheno.cfa %>% select(ph.nms$ph)
pheno.cfa.sel %>% is.na() %>% colSums()

# Colnames 
colnames(pheno.efa.sel) <- colnames(pheno.cfa.sel) <- ph.nms$names

# ------------------------------------------- #
# 06. Run EFA on training sample ------------
# ------------------------------------------- #
modelEFA <- paste0(paste0("efa('efa')*f",1:numFacts, collapse = " + "), " =~ ",
                   paste0(ph.nms$names, collapse = " + "), "\n",  paste0("f",1:numFacts, "~~1*f", 1:numFacts, collapse = "\n"))

cat(modelEFA)

## 6.1. Correlated factors: Oblimin rotation
efaOBL <- sem(model = modelEFA, data = pheno.efa.sel, estimator = "ML", missing = "ML", rotation = "oblimin")
summary(object = efaOBL, standardized = TRUE)


## 6.2. Uncorrelated latent factors: varimax rotation
efaORT <- sem(model = modelEFA, data = pheno.efa.sel, estimator="ML", missing="ML", rotation = "varimax", rotation.args = list(orthogonal = TRUE))
summary(object = efaORT, standardized = TRUE)

# ------------------------------------------- #
# 07. Run CFA on test sample ----------------
# ------------------------------------------- #
cutoff<- 0.1

## 7.1. Correlated factors: Oblimin rotation
efaOBL.values <- lavInspect(efaOBL, what= "est")$lambda %>% as.data.frame
efaOBL.values$trait <- factor(rownames(efaOBL.values), levels = c("lang_dis","lang_lvl","bhv_odd","age_crawl","selffeed","contr_mv","selfinj","sameness"))


model.cfaOBL <- efaOBL.values %>% 
  gather(key = "fact", value = "est", -trait) %>% 
  # Apply cut-off
  filter(abs(est)>cutoff) %>% 
  # Reformat to get a matrix where each row is a factor
  mutate(est=NULL) %>% group_by(fact) %>%
  summarise(outp = paste0(trait, collapse = " + ")) %>%
  # Output a vector that defines the model
  summarise(outp = paste0(fact, "=~ NA*", outp, collapse = "\n")) %>% as.character()

cat(model.cfaOBL)
# f1=~ NA*lang_dis + lang_lvl + contr_mv
# f2=~ NA*lang_dis + bhv_odd + age_crawl + selffeed + contr_mv
# f3=~ NA*bhv_odd + contr_mv + selfinj + sameness

# Fit in lavaan
cfaOBL <- sem(model = model.cfaOBL, 
              data = pheno.cfa.sel, 
              estimator="ML", missing="ML", 
              rotation = "oblimin")
summary(object = cfaOBL, standardized=TRUE)



## 7.2. Uncorrelated latent factors: varimax rotation
efaORT.values <- lavInspect(efaORT, what= "est")$lambda %>% as.data.frame
efaORT.values$trait <- factor(rownames(efaORT.values), levels = c("lang_dis","lang_lvl","bhv_odd","age_crawl","selffeed","contr_mv","selfinj","sameness"))


model.cfaORT <- efaORT.values %>% 
  gather(key = "fact", value = "est", -trait) %>% 
  # Apply cut-off
  filter(abs(est)>cutoff) %>% 
  # Reformat to get a matrix where each row is a factor
  mutate(est=NULL) %>% group_by(fact) %>%
  summarise(outp = paste0(trait, collapse = " + ")) %>%
  # Output a vector that defines the model
  summarise(outp = paste0(fact, "=~ NA*", outp, collapse = "\n")) %>% as.character()

cat(model.cfaORT)
# f1=~ NA*lang_dis + lang_lvl + selffeed + contr_mv + selfinj
# f2=~ NA*lang_dis + bhv_odd + age_crawl + selffeed + contr_mv + sameness
# f3=~ NA*bhv_odd + contr_mv + selfinj + sameness


# Fit in lavaan
cfaORT <- sem(model = model.cfaORT, 
              data = pheno.cfa.sel, 
              estimator="ML", missing="ML", 
              rotation = "varimax", orthogonal = TRUE)
summary(object = cfaORT, standardized=TRUE)


# ------------------------------------------- #
# 08. Check CFA model fit  ------------------
# ------------------------------------------- #
modelFits <-  data.frame(rbind(lavInspect(cfaOBL, "fit"), 
                               lavInspect(cfaORT, "fit")),
                         row.names = c(paste0("CFA.", numFacts,"f.oblimin"),paste0("CFA.", numFacts,"f.varimax")))

modelFits.summary <- round(modelFits[c("cfi","tli","rmsea","srmr","aic","bic")],2)
kable(modelFits.summary) %>% kable_paper("hover", full_width = T, html_font = "arial")

# ------------------------------------------- #
# 09. Fit best model in full sample ---------
# ------------------------------------------- #
# Run step 7 using full sample instead with the rotation that provided the best-fit





