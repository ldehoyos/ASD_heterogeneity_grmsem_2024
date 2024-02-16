##################################################-
## Project: autism grmsem 2024
## Description: compare grmsem models based on model fit parameters
## Date: September 2023
## Author: Lucia de Hoyos
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism

rm(list = ls())


# ------------------------------------------- #
# 01. Load libraries and functions ----------
# ------------------------------------------- #

fnames.models <- c(Cholesky="Cholesky_SPARK_model.RData",
                   IP_1f="IP_1factor_SPARK_model.RData",
                   IPC_1f="IPC_1factor_SPARK_model.RData",
                   IPC_3f_bifactor="IPC_3factor_bifactor_SPARK_model.RData",
                   IPC_3f="IPC_3factor_SPARK_model.RData")

# ------------------------------------------- #
# 02. Get model fit parameters for each model ----
# ------------------------------------------- #

df <- data.frame(model=names(fnames.models),logLik=NA,Npar=NA,AIC=NA,BIC=NA)

for(i in 1:length(fnames.models)){
  
  # Load the model
  modelIn <- new.env()
  load(fnames.models[i], envir = modelIn)
  
  # Parameters
  df$Npar[i] <- sum(modelIn$out$model.in$freepar)
  
  # Model fit estimates
  df$logLik[i] <- modelIn$out$model.fit$LL
  df$AIC[i] <- modelIn$out$model.fit$LL.AIC
  df$BIC[i] <- modelIn$out$model.fit$LL.BIC
  
  rm(modelIn)
  
}

# ------------------------------------------- #
# 03. Run likelihood ratio test (LRT) with Cholesky ----
# ------------------------------------------- #

df2 <- cbind(df, data.frame(LRT_cholesky_Chi=NA,LTR_cholesky_devPAR=NA,LRT_cholesky_P=NA))

for(j in grep("Cholesky", fnames.models, invert = T)){
  
  # Chi difference with Cholesky
  df2$LRT_cholesky_Chi[j] <- devChi <- -2*(df2$logLik[j] - df2$logLik[df2$model=="Cholesky"])
  
  # Deviance in fitted parameters (parsat - parmodel)
  df2$LTR_cholesky_devPAR[j] <- devPAR <- df2$Npar[df2$model=="Cholesky"] - df2$Npar[j]
  
  # Likelihood ratio test (LRT) Cholesky versus IPC
  df2$LRT_cholesky_P[j] <- pchisq(q = devChi, df = devPAR, lower.tail = F)
  
  rm(devPAR,devChi)
}

# ------------------------------------------- #
# 04. Run likelihood ratio test (LRT) with bifactor ----
# ------------------------------------------- #

df3 <- cbind(df2, data.frame(LRT_bifactor_Chi=NA,LTR_bifactor_devPAR=NA,LRT_bifactor_P=NA))

k <- which(df3$model=="IPC_3f")

# Chi difference with bifactor
df3$LRT_bifactor_Chi[k] <- devChi <- -2*(df3$logLik[k] - df3$logLik[df3$model=="IPC_3f_bifactor"])

# Deviance in fitted parameters
df3$LTR_bifactor_devPAR[k] <- devPAR <- df3$Npar[df3$model=="IPC_3f_bifactor"] - df3$Npar[k]

# Likelihood ratio test (LRT) Cholesky versus IPC
df3$LRT_bifactor_P[k] <- pchisq(q = devChi, df = devPAR, lower.tail = F)

rm(devChi, devPAR)

