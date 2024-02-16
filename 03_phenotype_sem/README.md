# 05. Phenotype SEM analysis

Code use to run structural equation modelling (SEM) analysis across final phenotypes in SPARK model, see `phenotype_sem_analysis.R`.

The pipeline is analogous to the data-driven modelling approach used at the genetic level, and has the following steps:

1. Data checks, filtering missingness on more than half of the phenotypes
2. PCA check to determine number of factors in the full sample
3. Split sample in two halves matching on missigness and sex
4. Run exploratory factor analysis (EFA) using varimax and oblimin rotation
5. Run confirmatory factor analysis (CFA) using varimax and oblimin rotation
6. Check which CFA model fits the data best
7. Run best-fitting model in the full sample


nFactors version: 2.4.1

lavaan version: 0.6-10

References
- lavaan: An R Package for Structural Equation Modeling [link to paper](https://doi.org/10.18637/jss.v048.i02)
- Non-Graphical Solutions for Cattell's Scree Test [link to paper](https://doi.org/10.1027/1614-2241/a000051)
