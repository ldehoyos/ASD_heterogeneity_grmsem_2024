##################################################-
## Project: autism grmsem 2024
## Description: run grmsem IPC 3-factor model with PGS
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
out_nm <- "3factor_SPARK_model_PGS"

# Phenotypes to include in the model
ph.sel=c("EA_Lee18.prscs.pgs.z","dev_lang_dis.drcat","language_level.rnkrq","behav_odd.drcat",
         "crawled_age_mos.rnkrq","fed_self_spoon_age_mos.rnkrq","control_during_movement.rnkrq",
         "ii_self_injurious_score.rnkrq","v_sameness_behavior_score.rnkrq")

# Free parameters, genetic part
A.v.free=c(1,1,1,0,1,1,0,1,0,1,1,0,0,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,0)

# Starting values, genetic part
A.v.startval=c(0.046,-0.346,0.466,0,-0.152,-0.368,0,0.153,0,-0.012,-0.202,0,0,0.476,0.189,-0.334,0.366,0,-0.151,-0.185,0,0.429,-0.121,0,0,0.228,0.426,-0.936,0,0,0,0,0.349,0.338,0,0)

# Free parameters, residual part
E.v.free=NULL

# Starting values, residual part
E.v.startval=c(0.075,-0.163,0.022,-0.234,-0.056,-0.074,-0.284,-0.785,-0.168,0.883,-0.104,-0.01,0.087,0.016,-0.257,0.09,0.068,0.879,0.083,0.067,-0.004,0.244,-0.201,0.01,0.868,0.003,-0.086,-0.039,-0.138,-0.127,0.849,0.141,-0.04,-0.167,0.07,0.823,-0.212,0.025,0.037,0.721,-0.245,-0.265,0.13,0.829,0.072)

# ------------------------------------------- #
# 03. Load data -----------------------------
# ------------------------------------------- #

# Load the GRM matrix G
load("grm.spark_grm_pheno.pgs.RData")

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
