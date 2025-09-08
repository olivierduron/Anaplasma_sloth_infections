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
- `bloodparasite` : Combined infection status for blood parasites (microfilariae + trypanosome + _Babesia_, but excluding _Anaplasma_; 0: Uninfected; 1: Infected)
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

Fit a GLM to test whether `anaplasma` is influenced by interactions among `sex`, `age`, `season`, `tick`, and `bloodparasite` in Bt:
```
model_1 <- glm(anaplasma ~ sex * age * season * tick * bloodparasite, data = data_Bt, family = binomial)
```

Fit a GLM to test whether `anaplasma` infection prevalence is influenced by additive effects of `sex`, `age`, `season`, `tick`, and `bloodparasite` in Bt:
```
model_1a <- glm(anaplasma ~ sex + age + season + tick + bloodparasite, data = data_Bt, family = binomial)
```

Compare the additive model (model_1a) to the interaction model (model_1) using a likelihood ratio test:
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

Tests for associations between `anaplasma` and the presence of blood parasites (`microfilaria`, `trypanosome`, `babesia`) considered separately in Bt:
```
fisher.test(table(data_Bt$anaplasma, data_Bt$microfilaria))  
fisher.test(table(data_Bt$anaplasma, data_Bt$trypanosome))  
fisher.test(table(data_Bt$anaplasma, data_Bt$babesia))
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Bt$anaplasma, data_Bt$microfilaria)
p-value = 0.8227
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3362994 2.4250879
sample estimates:
odds ratio 
0.894268 

data:  table(data_Bt$anaplasma, data_Bt$trypanosome)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.0591107 71.5824945
sample estimates:
odds ratio 
1.176551 

data:  table(data_Bt$anaplasma, data_Bt$babesia)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0 Inf
sample estimates:
odds ratio 
0
```

## Step 5. Test whether _Anaplasma_ infection prevalence is influenced in _Choloepus didactylus_ (Cd) by sex, age, season, ticks and blood parasites

Create a subset `data_Cd` containing only records for _Choloepus didactylus_ (Cd):
```
data_Cd <- subset(data_sloth, species == "Cd")
```

Fit a GLM to test whether `anaplasma` is influenced by interactions among `sex`, `age`, `season`, `tick`, and `bloodparasite` in Cd:
```
model_2 <- glm(anaplasma ~ sex * age * season * tick * bloodparasite, data = data_Cd, family = binomial)
```

Fit a GLM to test whether `anaplasma` infection prevalence is influenced by additive effects of `sex`, `age`, `season`, `tick`, and `bloodparasite` in Cd:
```
model_2a <- glm(anaplasma ~ sex + age + season + tick + bloodparasite, data = data_Cd, family = binomial)
```

Compare the additive model (model_2a) to the interaction model (model_2) using a likelihood ratio test:
```
anova(model_2a, model_2, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: anaplasma ~ sex + age + season + tick + bloodparasite
Model 2: anaplasma ~ sex * age * season * tick * bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        77    106.525                     
2        61     85.459 16   21.066    0.176
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_2, model_2a)
```

Results are:
```
         df      AIC
model_2  22 129.4592
model_2a  6 118.5250
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_2a, test = "Chisq")
res
```

Results are:
```
Single term deletions
Model:
anaplasma ~ sex + age + season + tick + bloodparasite
              Df Deviance    AIC     LRT Pr(>Chi)  
<none>             106.53 118.53                   
sex            1   107.12 117.12 0.59812  0.43930  
age            1   108.07 118.07 1.54781  0.21346  
season         1   107.82 117.82 1.29158  0.25576  
tick           1   109.35 119.35 2.82844  0.09261
bloodparasite  1   107.03 117.03 0.50223  0.47852  
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_2a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        118.53   0.00000
sex           117.12  -1.40188
age           118.07  -0.45219
season        117.82  -0.70842
tick          119.35   0.82844
bloodparasite 117.03  -1.49777
```

Tests for associations between `anaplasma` and the presence of blood parasites (`microfilaria`, `trypanosome`, `babesia`) considered separately in Cd:
```
fisher.test(table(data_Cd$anaplasma, data_Cd$microfilaria))  
fisher.test(table(data_Cd$anaplasma, data_Cd$trypanosome))  
fisher.test(table(data_Cd$anaplasma, data_Cd$babesia))
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Cd$anaplasma, data_Cd$microfilaria)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.2272016 3.5201430
sample estimates:
odds ratio 
0.9086257 

data:  table(data_Cd$anaplasma, data_Cd$trypanosome)
p-value = 0.2292
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.2028662       Inf
sample estimates:
odds ratio 
Inf 

data:  table(data_Cd$anaplasma, data_Cd$babesia)
p-value = 0.2538
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.4399427 32.1559681
sample estimates:
odds ratio 
2.892291 
```

## Step 6. Test xxxx

xxx
