# 03. GRM-SEM analysis

Code use to run GRM-SEM (genetic-relationship-matrix structural equation modelling). 

Within this study we implemented a novel data-driven approach to model the multivariate architecture of several phenotypes. This approach consists of 7 steps.


## STEP 1. DESCRIBE THE FULL GENETIC ARCHITECTURE

To describe the architecture of the phenotypes modelled, we fit a saturated (Cholesky) model to the data. See `01_grmsem_Cholesky_model.R`.

## STEP 2. PREDICT NUMBER OF SHARED GENETIC FACTORS

To identify how many genetic factors are needed to best describe the genetic architecture of these phenotypes, we run an eigenvalue decomposition of the genetic correlation matrix from the Cholesky model. See `02_predict_nAC.R`.

## STEP 3. APPROXIMATE GENETIC FACTOR STRUCTURE

To uncover the optimal structure between the genetic factors and the phenotypes, we fit an exploratory factor analysis (EFA) in lavaan to the estimates Cholesky genetic covariance matrix. See `03_lavaan_efa.R`.

## STEP 4. DEFINE MULTI-FACTOR MODELS

Using the information from steps 2 and 3, in GRM-SEM, we fit an IPC (hybrid independent pathway / Cholesky) model to the data. See `04_grmsem_IPC_3factor_model.R`.

To confirm the independence of genetic factors, we also fit a bifactor model. See `04_grmsem_IPC_bifactor_model.R`.

## STEP 5. DEFINE ONE-FACTOR MODELS

We also fit a one-factor IP (independent pathway) and IPC model to the data, using GRM-SEM. See `05_grmsem_IP_1factor_model.R` and `05_grmsem_IPC_1factor_model.R`.

## STEP 6. IDENTIFY BEST-FITTING MODEL

We compare all the fitted models in terms of LRT, AIC, BIC and SRMR. See `06_model_comparison.R`.

## STEP 7. CHARACTERISE IDENTIFIED SHARED GENETIC FACTORS

Once the best fitting model has been identified, we mapped the factor structure to:
   - liability to Asperger. See `07a1_merge_ASDsubcategories_with_grm.R` and `07a2_grmsem_IPC_3factor_model_with_ASDsubcategories.R`.
   - PGS for educational attainment (EA3, Lee et al. 2018). See `07b1_merge_PGS_with_grm.R` and `07b2_grmsem_IPC_3factor_model_with_pgs.R`.


In this analysis, we used grmsem version 1.1.2

References and links
- Developmental Changes Within the Genetic Architecture of Social Communication Behavior: A Multivariate Study of Genetic Variance in Unrelated Individuals  [link to paper](https://doi.org/10.1016/j.biopsych.2017.09.020)
- https://gitlab.gwdg.de/beate.stpourcain/grmsem