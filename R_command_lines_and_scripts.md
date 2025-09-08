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
```
# Load the dateset

data_sloth <- read.csv("https://raw.githubusercontent.com/olivierduron/Anaplasma_sloth_infections/main/data_sloth.csv", sep="\t")

# Complete overview

data_sloth
```


## Step 2. Prepare the data for analysis
```
# Convert categorical variables into factors
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

# Load libraries for analysis
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
```
# Calculate Anaplasma infection prevalence and 95% confidence interval for Bradypus tridactylus (Bt) and Choloepus didactylus (Cd)
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

```
#Test if Anaplasma infection prevalence is influenced by sloth species
chisq.test(table(data_sloth$anaplasma, data_sloth$species))
```

Results are:
```
Pearson's Chi-squared test with Yates' continuity correction
data:  table(data_sloth$anaplasma, data_sloth$species)
X-squared = 3.3261, df = 1, p-value = 0.06819
```

## Step 3. 
