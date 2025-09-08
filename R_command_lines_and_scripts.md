# **R command lines and script**

We analyzed data from 175 wild sloths captured between 1994 and 1995 during the flooding of the Petit Saut Dam (5°03′43″ N, 53°03′00″ O) on the Sinnamary River (French Guiana, South America). The clinical data include the following variables for each examined sloth: 
- `species` : Sloth species (Bt: *Bradypus tridactylus*; Cd: *Choloepus didactylus*)
- `sex` : Sex of the sloth (F: Female; M: Male)
- `age_class` : Age category (A: Adult; J: Juvenile)
- `season` : Season of capture (W: Wet; D: Dry)
- `weight` : Body weight (quantitative variable, in kg)
- `total_length` : Total body length (quantitative variable, in cm)
- `wither_height` : Height at the withers (quantitative variable, in cm)
- `neck_size` : Neck circumference (quantitative variable, in cm)
- `temperature` : Body temperature (quantitative variable, in °C)
- `hematocrit` : Hematocrit level (quantitative variable, in %)
- `health_condition` : Overall health status (G: Good; D: Deteriorated)
- `anaplasma` : Infection status with *Anaplasma* (0: Uninfected; 1: Infected)
- `tick` : Presence of ticks in the fur (0: Absent; 1: Present)
- `microfilaria` : Infection status with microfilariae (0: Uninfected; 1: Infected)
- `trypanosome` : Infection status with trypanosomes (0: Uninfected; 1: Infected)
- `babesia` : Infection status with _Babesia_ (0: Uninfected; 1: Infected)
- `blood_parasite` : Combined infection status for blood parasites (microfilariae + trypanosome + _Babesia_, but excluding _Anaplasma_; 0: Uninfected; 1: Infected)
Details about all the experimental methods and measures are available in the related manuscript.


## Step 1. Retrieving the data

All veterinary clinical data for the two sloth species are available here: https://github.com/olivierduron/Anaplasma_sloth_infections/blob/main/data_sloth.csv

This database will be referred to as `data_sloth` throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R:
```
data_sloth <- read.csv("https://raw.githubusercontent.com/olivierduron/Anaplasma_sloth_infections/main/data_sloth.csv", sep="\t")
```


## Step 2. Prepare the data for analysis

Convert categorical variables into factors:
```
data_sloth$anaplasma      <- as.factor(data_sloth$anaplasma)
data_sloth$species        <- as.factor(data_sloth$species)
data_sloth$season         <- as.factor(data_sloth$season)
data_sloth$sex            <- as.factor(data_sloth$sex)
data_sloth$age            <- as.factor(data_sloth$age)
data_sloth$tick           <- as.factor(data_sloth$tick)
data_sloth$microfilaria   <- as.factor(data_sloth$microfilaria)
data_sloth$trypanosome    <- as.factor(data_sloth$trypanosome)
data_sloth$babesia        <- as.factor(data_sloth$babesia)
data_sloth$bloodparasite  <- as.factor(data_sloth$bloodparasite)
```

Load libraries for analysis: 
```
library(binom)
library(dplyr)
library(MASS)
library(ggplot2)
library(patchwork)
library(smatr)
library(lmtest)
library(akima)
library(pwr)
library(survival)
library(RColorBrewer)
```

## Step 3. Calculate *Anaplasma* infection prevalence
Calculate _Anaplasma_ infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd):

```
prevalence_results <- data_sloth %>% group_by(species) %>% summarise(n = n(), positives = sum(anaplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results are:
```
# A tibble: 2 × 6
  species     n positives prevalence conf_low conf_high
  <fct>   <int>     <int>      <dbl>    <dbl>     <dbl>
1 Bt         92        58      0.630    0.523     0.729
2 Cd         83        40      0.482    0.371     0.594
```

Test if `anaplasma` is influenced by sloth `species`:
```
chisq.test(table(data_sloth$anaplasma, data_sloth$species))
```

Results are:
```
Pearson's Chi-squared test with Yates' continuity correction
data:  table(data_sloth$anaplasma, data_sloth$species)
X-squared = 3.3261, df = 1, p-value = 0.06819
```

## Step 4. Test whether _Anaplasma_ infection prevalence is influenced in _Bradypus tridactylus_ (Bt) by sex, age, season, ticks and blood parasites
Create a subset `data_Bt` containing only records for _Bradypus tridactylus_ (Bt):

```
data_Bt <- subset(data_sloth, species == "Bt")
```

Fit a GLM to test whether `anaplasma` is influenced by interactions among `sex`, `age`, `season`, `tick`, and `blood_parasite` in Bt:
```
model_1 <- glm(anaplasma ~ sex * age * season * tick * bloodparasite, data = data_Bt, family = binomial)
summary(model_1)
```

Results are:
```
glm(formula = anaplasma ~ sex * age * season * tick * bloodparasite, 
    family = binomial, data = data_Bt)
Coefficients: (14 not defined because of singularities)
                                         Estimate Std. Error z value Pr(>|z|)
(Intercept)                             2.877e-01  7.638e-01   0.377    0.706
sexM                                    2.231e-01  1.057e+00   0.211    0.833
ageJ                                   -1.757e+01  3.956e+03  -0.004    0.996
seasonW                                 5.878e-01  9.309e-01   0.631    0.528
tick1                                  -1.785e+01  3.956e+03  -0.005    0.996
bloodparasite1                         -2.877e-01  1.118e+00  -0.257    0.797
sexM:ageJ                               1.706e+01  3.956e+03   0.004    0.997
sexM:seasonW                            4.529e-15  1.317e+00   0.000    1.000
ageJ:seasonW                            1.808e+01  3.956e+03   0.005    0.996
sexM:tick1                              3.491e+01  5.595e+03   0.006    0.995
ageJ:tick1                             -3.462e+01  4.845e+03  -0.007    0.994
seasonW:tick1                          -5.878e-01  4.568e+03   0.000    1.000
sexM:bloodparasite1                     1.163e+00  1.742e+00   0.668    0.504
ageJ:bloodparasite1                            NA         NA      NA       NA
seasonW:bloodparasite1                 -7.696e-02  1.438e+00  -0.054    0.957
tick1:bloodparasite1                   -8.755e-01  4.845e+03   0.000    1.000
sexM:ageJ:seasonW                              NA         NA      NA       NA
sexM:ageJ:tick1                                NA         NA      NA       NA
sexM:seasonW:tick1                     -1.687e+01  6.043e+03  -0.003    0.998
ageJ:seasonW:tick1                             NA         NA      NA       NA
sexM:ageJ:bloodparasite1                       NA         NA      NA       NA
sexM:seasonW:bloodparasite1            -1.897e+00  2.164e+00  -0.877    0.381
ageJ:seasonW:bloodparasite1                    NA         NA      NA       NA
sexM:tick1:bloodparasite1                      NA         NA      NA       NA
ageJ:tick1:bloodparasite1                      NA         NA      NA       NA
seasonW:tick1:bloodparasite1            1.885e+01  5.595e+03   0.003    0.997
sexM:ageJ:seasonW:tick1                        NA         NA      NA       NA
sexM:ageJ:seasonW:bloodparasite1               NA         NA      NA       NA
sexM:ageJ:tick1:bloodparasite1                 NA         NA      NA       NA
sexM:seasonW:tick1:bloodparasite1              NA         NA      NA       NA
ageJ:seasonW:tick1:bloodparasite1              NA         NA      NA       NA
sexM:ageJ:seasonW:tick1:bloodparasite1         NA         NA      NA       NA
(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 121.21  on 91  degrees of freedom
Residual deviance: 100.33  on 74  degrees of freedom
AIC: 136.33
Number of Fisher Scoring iterations: 16
```

Fit a GLM to test whether `anaplasma` infection prevalence is influenced by additive effects of `sex`, `age`, `season`, `tick`, and `blood_parasite` in Bt:
```
model_1a <- glm(anaplasma ~ sex + age + season + tick + bloodparasite, data = data_Bt, family = binomial)
summary(model_1a)
```

Results are:
```
glm(formula = anaplasma ~ sex + age + season + tick + bloodparasite, 
    family = binomial, data = data_Bt)
Coefficients:
               Estimate Std. Error z value Pr(>|z|)
(Intercept)     0.19121    0.47792   0.400    0.689
sexM            0.73449    0.45561   1.612    0.107
ageJ           -1.10598    0.99102  -1.116    0.264
seasonW         0.20467    0.46811   0.437    0.662
tick1          -0.57129    0.62584  -0.913    0.361
bloodparasite1 -0.03394    0.46292  -0.073    0.942
(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 121.21  on 91  degrees of freedom
Residual deviance: 116.74  on 86  degrees of freedom
AIC: 128.74
Number of Fisher Scoring iterations: 4
```

Compare the additive model to the interaction model using a likelihood ratio test:
```
anova(model_1a, model_1, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: anaplasma ~ sex + age + season + tick + bloodparasite
Model 2: anaplasma ~ sex * age * season * tick * bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        86     116.73                     
2        74     100.33 12   16.409   0.1732
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_1, model_1a)
```

Results are:
```
         df      AIC
model_1  18 136.3264
model_1a  6 128.7350
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_1a, test = "Chisq")
res
```

Results are:
```
Single term deletions
Model: anaplasma ~ sex + age + season + tick + bloodparasite
              Df Deviance    AIC     LRT Pr(>Chi)
<none>             116.73 128.74                 
sex            1   119.39 129.40 2.65956   0.1029
age            1   118.01 128.01 1.27259   0.2593
season         1   116.93 126.93 0.19048   0.6625
tick           1   117.56 127.56 0.82359   0.3641
bloodparasite  1   116.74 126.74 0.00537   0.9416
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_1a)
res$delta_AIC <- res$AIC - aic_full
res$delta_AIC
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        128.74   0.00000
sex           129.40   0.65956
age           128.01  -0.72741
season        126.93  -1.80952
tick          127.56  -1.17641
bloodparasite 126.74  -1.99463
```

## Step 4. Test whether _Anaplasma_ infection prevalence is influenced in Choloepus didactylus (Cd) by sex, age, season, ticks and blood parasites

Create a subset `data_Cd` containing only records for Choloepus didactylus (Cd):
```
data_Cd <- subset(data_sloth, species == "Cd")
```



