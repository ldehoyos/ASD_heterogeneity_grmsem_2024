# 01. Phenotype transformations

Code use to run phenotype transformations: `phenotype_transformations.R` file



## Categorical phenotypes

Categorical phenotypes were transformed using deviance residuals. These are the residuals of the logistic regression of the phenotype on the covariates.

```
# F1 <- trait ~ sex + age + age^2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10

binom.reg <- glm(F1, data = dset, na.action = na.exclude, family = "binomial")
out_score <- scale(resid(binom.reg))
```

## Continuous phenotypes

Continuous phenotypes were transformed using residualised rank-transformed residuals.

First, we run a linear regression of the phenotype on the covariates and the residuals of that regression are rank-transformed.
```
# F1 <- trait ~ sex + age + age^2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10

ols.reg1 <- lm(F1, data = dset, na.action = na.exclude)
res_ols.reg1 <- resid(ols.reg1)
std_cont  <- qnorm((rank(res_ols.reg1,na.last="keep")-0.5)/sum(!is.na(res_ols.reg1)))
```

Then, we run a linear regression on the rank transformed residuals. And get the residuals of that regression.  
```
# F2 <- trait.rnk  ~ sex + age + age^2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10

ols.reg2  <- lm(F2, data = dset, na.action = na.exclude)
out_score <- scale(resid(ols.reg2))
```