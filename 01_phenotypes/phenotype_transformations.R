##################################################-
## Project: autism grmsem 2024
## Description: Functions to transform phenotypes
## Date: September 2023
## Author: Lucia de Hoyos, Fenja Schlag, Ellen Verhoef, Beate St Pourcain
## Max Planck Institute for Psycholinguistics, Nijmegen, The Netherlands
##################################################-
## de Hoyos et al. Structural models of genome-wide covariance identify multiple common dimensions in autism


# dset     [dataframe] dataset that includes your trait variable and your covariates
# traitVar [character] trait (or phenotype) to transform
# covVars  [character vector] vector listing the variable name of the covariates to include (in this study is PC, sex, age, age^2)


# Deviance residuals
drcat <- function(dset, traitVar, covVars){
  
  # Get covariates
  cov <- covVars[!is.na(covVars)]
  
  # Run binomial regression 
  F1   <- as.formula(paste(traitVar, paste(cov, collapse = " + "), sep = " ~ "))
  print(F1)
  binom.reg <- glm(formula = F1, data = dset, na.action = na.exclude, family = "binomial")
  
  # Get deviance residuals
  res_binom.reg  <- resid(binom.reg)
  
  # Scale residuals
  res_dich <- scale(res_binom.reg)
  return(res_dich)
}


# Function twice residualised rank-transformed residuals
rqrnkrq <- function(dset, traitVar, covVars){
  
  # Get covariates
  cov <- covVars[!is.na(covVars)]
  
  # Run linear regression
  F1 <- as.formula(paste(traitVar,paste(cov,collapse=" + "),sep=" ~ "))
  print(F1)
  ols.reg1 <- lm(F1, data = dset, na.action = na.exclude)
  
  # Get residuals
  res_ols.reg1 <- resid(ols.reg1)
  
  # Standardize residuals (ranked residuals)/N
  std_cont  <- qnorm((rank(res_ols.reg1,na.last="keep")-0.5)/sum(!is.na(res_ols.reg1)))
  
  # Add column to dataset
  dset   <- as.data.frame(cbind(dset,std_cont))
  rnkVar <- paste0(traitVar[i],".rnk")
  colnames(dset)[ncol(dset)] <- rnkVar
  
  # Run linear regression on standardized residuals
  F2 <- as.formula(paste(rnkVar,paste(cov,collapse=" + "),sep=" ~ "))
  print(F2)
  ols.reg2 <- lm(F2, data = dset, na.action = na.exclude)
  
  # Get residuals
  res_ols.reg2 <- resid(ols.reg2)
  
  # Scale residuals
  dev_std_cont <- scale(res_ols.reg2)
  return(dev_std_cont)
}




# Read in phenotype file
pheno <- read.table(paste0(workDir,ph.file), sep = "\t", na.strings = "NA", header=TRUE)

# Define sex variable
sexVar <- "Sex"

# Define ancestry-informative principal components
pcVar <- paste0("pc", 1:10)


# -------------------------------------------- #
# Run transformations on categorical phenotypes
# -------------------------------------------- #

# Define the phenotypes and age variables
cat_vars <- data.frame(var="SPARK_BMS_behav_adhd",ageVar="age_years_bms",ageSqVar="age_years_bms_sq")

for (i in 1:nrow(cat_vars)){
  # Create variable name
  nameVar <- paste0(cat_vars$var[i],'.drcat.', 'pcsa')
  print(nameVar, quote = FALSE)
  
  # Select covariates
  cov_sel <- c(sexVar, pcVar, cat_vars$ageVar[i], cat_vars$ageSqVar[i])
  
  # Compute deviance residuals
  ph_tr <- drcat(dset = pheno, traitVar = cat_vars$var[i], covVars = cov_sel)
  
  # Add column to dataset
  pheno[nameVar] <- ph_tr
  
  rm(nameVar,ph_tr,cov_sel)
}

# -------------------------------------------- #
# Run transformations on continuous phenotypes
# -------------------------------------------- #

# Define the phenotypes and age variables
cont_vars <- data.frame(var="SPARK_SCQ_scq_final_score",ageVar="age_years_scq",ageSqVar="age_years_scq_sq")

for (i in 1:nrow(cont_vars)){
  # Create variable name
  nameVar <- paste0(cont_vars$var[i],'.rqrnkrq.', 'pcsa')
  print(nameVar, quote = FALSE)
  
  # Select covariates
  cov_sel <- c(sexVar, pcVar, cont_vars$ageVar[i], cont_vars$ageSqVar[i])
  
  # Compute deviance residuals
  ph_tr <- rqrnkrq(dset = pheno, traitVar = cont_vars$var[i], covVars = cov_sel)
  
  # Add column to dataset
  pheno[nameVar] <- ph_tr
  
  rm(nameVar,ph_tr,cov_sel)
}
