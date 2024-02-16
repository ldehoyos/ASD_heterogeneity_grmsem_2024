##################################################-
## Project: autism grmsem 2024
## Description: run grmsem IPC 3-factor bifactor model
## Date: September 2023
## Author: Lucia de Hoyos, Maria T (Mariska) Barendse & Beate St Pourcain
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
out_nm <- "3factor_bifactor_SPARK_model"

# Phenotypes to include in the model
ph.sel=c("dev_lang_dis.drcat","language_level.rnkrq","behav_odd.drcat",
         "crawled_age_mos.rnkrq","fed_self_spoon_age_mos.rnkrq","control_during_movement.rnkrq",
         "ii_self_injurious_score.rnkrq","v_sameness_behavior_score.rnkrq")

# Free parameters, genetic part
A.v.free=c(1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1)

# Starting values, genetic part
A.v.startval=c(0.358,-0.398,0,0.118,0.404,0,-0.191,0,-0.173,0,0,0.532,0.217,-0.343,0.316,0,-0.197,0,0.352,-0.171,0,0,0.204,0.438,0.098,0.25,0.289,0,0.311,0.374,0.269,0)

# Free parameters, residual part
E.v.free=NULL

# Starting values, residual part
E.v.startval=c(0.891,-0.1188,0.0511,0.0979,-0.0024,-0.2029,0.2225,0.1137,0.8717,0.0672,0.0743,-0.025,0.2482,-0.2473,-0.0762,0.8857,0.0361,-0.0569,0.0649,0.1124,-0.05,0.8463,0.1502,-0.0049,-0.1185,0.1125,0.8168,-0.1458,0.1053,0.0642,0.7828,0.0932,-0.1188,0.7692,0.3517,0.7906)

# ------------------------------------------- #
# 03. Load data -----------------------------
# ------------------------------------------- #

# Load the GRM matrix G
load("grm.spark_grm_pheno.RData")

# Subset your dataframe to the selected phenotypes
ph <- subset(x = ph.srt,select = ph.sel)

# ------------------------------------------- #
# 04. Fit and estimate the parameters -------
# ------------------------------------------- #

# Fit the model in grmsem
out <- grmsem.fit(
  ph = ph,
  G = G,
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
