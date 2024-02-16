##################################################-
## Project: autism grmsem 2024
## Description: run grmsem IPC 3-factor model with Asperger subcategory
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
model <- "IPC"

# Output filename
out_nm <- "3factor_SPARK_model_Asp"

# Phenotypes to include in the model
ph.sel=c("Asp.drcat","dev_lang_dis.drcat","language_level.rnkrq","behav_odd.drcat",
         "crawled_age_mos.rnkrq","fed_self_spoon_age_mos.rnkrq","control_during_movement.rnkrq",
         "ii_self_injurious_score.rnkrq","v_sameness_behavior_score.rnkrq")

# Free parameters, genetic part
A.v.free=c(1,1,1,0,1,1,0,1,0,1,1,0,0,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,0)

# Starting values, genetic part
A.v.startval=c(0.134,-0.355,0.451,0,-0.118,-0.404,0,0.172,0,0.034,-0.201,0,0,0.474,0.221,-0.332,0.357,0,0.058,-0.199,0,0.446,-0.154,0,0,0.171,0.382,0,0,0,0,0,0.284,0.385,0,0)

# Free parameters, residual part
E.v.free=NULL

# Starting values, residual part
E.v.startval=c(0.626,-0.281,0.168,0.033,0.002,-0.095,0.054,-0.079,-0.024,0.846,-0.06,0.056,0.104,-0.011,-0.196,0.218,0.092,0.875,0.076,0.043,0.022,0.248,-0.21,0.009,0.89,0.032,-0.062,0.053,0.083,-0.072,0.851,0.141,-0.02,-0.141,0.087,0.823,-0.177,0.119,0.057,0.778,0.117,-0.168,0.812,0.441,0.778)

# ------------------------------------------- #
# 03. Load data -----------------------------
# ------------------------------------------- #

# Load the GRM matrix G
load("grm.spark_grm_pheno.ASDsubcategories.RData")

# Subset your dataframe to the selected phenotypes
ph <- subset(x = ph.srt, select = ph.sel)

# ------------------------------------------- #
# 04. Fit and estimate the parameters -------
# ------------------------------------------- #

# Fit the model in grmsem
out <- grmsem.fit(
  ph = ph,
  G= G,
  model = model,
  LogL = TRUE,
  estSE = TRUE,
  n.AC = 3,
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
