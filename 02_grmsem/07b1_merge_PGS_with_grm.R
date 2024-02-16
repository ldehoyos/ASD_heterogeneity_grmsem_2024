##################################################-
## Project: autism grmsem 2024
## Description: merge PGS with GRM for grmsem analysis
## Date: September 2023
## Author: Lucia de Hoyos & Fenja Schlag
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())


# ------------------------------------------- #
# 01. Load data -----------------------------
# ------------------------------------------- #

# Phenotype data
pheno <- read.table(file = paste0(ph.dir,ph.file), header = TRUE)
f_PRS <- subset(x = pheno, select = c("family_id","subject_sp_id"))
colnames(f_PRS) <- c("FID","IID")

# Polygenic score data (Z-standardised)
prsData <- read.table(paste0(mainDir,args[i],"/EA_Lee18.prscs.scoring.input.pgs"), header = T, stringsAsFactors = F)
f_PRS <- merge(x = f_PRS, y = prsData, by = c("FID","IID"), all.x = TRUE)

# ------------------------------------------- #
# 02. Merge with GRM ------------------------
# ------------------------------------------- #

load("grm.spark_grm_pheno.RData")
grm.id$seq <- seq(1:nrow(grm.id))

ph.srt <- merge(x = grm.id, y = f_PRS, by = c("FID","IID"), all.x = T)
ph.srt <- ph.srt[order(ph.srt$seq),]
summary(ph.srt)

rm(list = setdiff(ls(), c("grm.id","G", "ph.srt")))
save.image("grm.spark_grm_pheno.pgs.RData")
