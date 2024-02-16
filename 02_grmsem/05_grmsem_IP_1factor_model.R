##################################################-
## Project: autism grmsem 2024
## Description: run grmsem IP 1-factor model
## Date: September 2023
## Author: Lucia de Hoyos & Beate St Pourcain
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())
set.seed(12345)

# ------------------------------------------- #
# 01. Load libraries and functions ----------
# ------------------------------------------- #

require(parallel)
require(numDeriv)
require(msm)
require(stats)
require(grmsem)

# ------------------------------------------- #
# 02. Define variables ----------------------
# ------------------------------------------- #

# Define model to run
model <- "IP"

# Output filename
out_nm <- "1factor_SPARK_model"

# Phenotypes to include in the model
ph.sel=c("dev_lang_dis.drcat","language_level.rnkrq","behav_odd.drcat",
         "crawled_age_mos.rnkrq","fed_self_spoon_age_mos.rnkrq","control_during_movement.rnkrq",
         "ii_self_injurious_score.rnkrq","v_sameness_behavior_score.rnkrq")

# Free parameters, genetic part
A.v.free=NULL

# Starting values, genetic part
A.v.startval=NULL

# Free parameters, residual part
E.v.free=NULL

# Starting values, residual part
E.v.startval=NULL

# ------------------------------------------- #
# 03. Load data -----------------------------
# ------------------------------------------- #

# Load the GRM matrix G
load("grm.spark_grm_pheno.RData")

# Subset your dataframe to the selected phenotypes
ph <- subset(x = ph.srt, select = ph.sel)

# ------------------------------------------- #
# 04. Fit and estimate the parameters -------
# ------------------------------------------- #

# Fit the model in grmsem
out <- grmsem.fit(
  ph = ph,
  G = G,
  LogL = TRUE,
  estSE = TRUE,
  model = model,
  n.AC = 1,
  A.v.free = A.v.free,
  A.v.startval = A.v.startval,
  E.v.startval = E.v.startval,
  E.v.free = E.v.free,
  cores = 4,
  compl.ph = FALSE,
  cluster = "FORK",
  optim = "ucminf"
)

print("GRM-SEM output: unstandardised parameters")
print(out)

# Save output to RData file
outfile <- paste0(model,"_",out_nm,".RData")
save.image(outfile)
