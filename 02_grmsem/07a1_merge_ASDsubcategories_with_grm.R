##################################################-
## Project: autism grmsem 2024
## Description: merge information on ASD subcategories with GRM for grmsem analysis
## Date: September 2023
## Author: Lucia de Hoyos
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

# ------------------------------------------- #
# 02. Create dummy variable -----------------
# ------------------------------------------- #

# Dummy code the diagnosis variable 
pheno$diagnosis <- as.factor(pheno$diagnosis)
table(pheno$diagnosis)

# Asperger's
pheno$Asp <- ifelse(pheno$diagnosis == "Asperger's Disorder", 1, 0)

# Autism
pheno$Aut <- ifelse(pheno$diagnosis == "Autism or Autistic Disorder", 1, 0)

# PDDNOS
pheno$PDDNOS <- ifelse(pheno$diagnosis == "Pervasive Developmental Disorder - Not Otherwise Specified (PDD-NOS)", 1, 0)

# Set ASD as NAs
na.pos <- which(pheno$diagnosis== 'Autism Spectrum Disorder')
pheno$Asp[na.pos] <- pheno$Aut[na.pos] <- pheno$PDDNOS[na.pos] <- NA

f_dummy <- subset(x = pheno, select=c("family_id","subject_sp_id","sex", paste0("PC",1:10), c("age_at_registration_years","age_at_registration_years_sq"),"Asp","Aut","PDDNOS"))
colnames(f_dummy) <- c("FID", "IID","Sex", paste0("PC",1:10), c("age_at_registration_years","age_at_registration_years_sq"),"Asp", "Aut", "PDDNOS")

# -------------------------------------------- #
# 03. Transform the variables ---------------
# ------------------------------------------- #
# Usign drcat function from 01_phenotypes/phenotype_transformations.R

# Define sex variable
sexVar <- "Sex"

# Define ancestry-informative principal components
pcVar <- paste0("PC", 1:10)

# Define the phenotypes and age variables
cat_vars <- data.frame(var= c("Asp", "Aut", "PDDNOS"), ageVar="age_at_registration_years", ageSqVar="age_at_registration_years_sq")

for (i in 1:nrow(cat_vars)){
  # Create variable name
  nameVar <- paste0(cat_vars$var[i],'.drcat.', 'pcsa')
  print(nameVar, quote = FALSE)
  
  # Select covariates
  cov_sel <- c(sexVar, pcVar, cat_vars$ageVar[i], cat_vars$ageSqVar[i])
  
  # Compute deviance residuals
  ph_tr <- drcat(dset = f_dummy, traitVar = cat_vars$var[i], covVars = cov_sel)
  
  # Add column to dataset
  f_dummy[nameVar] <- ph_tr
  
  rm(nameVar,ph_tr,cov_sel)
}

# ------------------------------------------- #
# 04. Merge with GRM ------------------------
# ------------------------------------------- #

load("grm.spark_grm_pheno.RData")
grm.id$seq <- seq(1:nrow(grm.id))

ph.srt <- merge(x = grm.id, y = f_dummy, by = c("FID","IID"), all.x = T)
ph.srt <- ph.srt[order(ph.srt$seq),]
summary(ph.srt)

rm(list = setdiff(ls(), c("grm.id","G", "ph.srt")))
save.image("grm.spark_grm_pheno.ASDsubcategories.RData")
