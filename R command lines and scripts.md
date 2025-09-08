We analyzed data from 175 wild sloths captured between 1994 and 1995 during the flooding of the Petit Saut Dam. The clinical data include the following variables for each examined sloth: 
- species: Sloth species (Bt: *Bradypus tridactylus*; Cd: *Choloepus didactylus*)
- sex: Sex of the sloth (F: Female; M: Male)
- age_class: Age category (A: Adult; J: Juvenile)
- season: Season of capture (W: Wet; D: Dry)
- weight: Body weight (quantitative variable, in kg)
- total_length: Total body length (quantitative variable, in cm)
- wither_height: Height at the withers (quantitative variable, in cm)
- neck_size: Neck circumference (quantitative variable, in cm)
- temperature: Body temperature (quantitative variable, in °C)
- hematocrit: Hematocrit level (quantitative variable, in %)
- health_condition: Overall health status (G: Good; D: Deteriorated)
- anaplasma: Infection status with *Anaplasma* (0: Uninfected; 1: Infected)
- tick: Presence of ticks in the fur (0: Absent; 1: Present)
- microfilaria: Infection status with microfilariae (0: Uninfected; 1: Infected)
- trypanosome: Infection status with trypanosomes (0: Uninfected; 1: Infected)
- babesia: Infection status with _Babesia_ (0: Uninfected; 1: Infected)
- blood_parasites: Combined infection status for blood parasites (microfilariae + trypanosome + _Babesia_, but excluding _Anaplasma_; 0: Uninfected; 1: Infected)
Details about all the experimental methods and measures are available in the related manuscript.

Step 1. Retrieving the data

All veterinary clinical data for the two sloth species are available here: https://github.com/olivierduron/Anaplasma_sloth_infections/blob/main/data_sloth.csv
This database will be referred to as data_sloth throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

#Load the dateset
data_sloth <- read.csv("https://raw.githubusercontent.com/olivierduron/Anaplasma_sloth_infections/main/data_sloth.csv", sep="\t")
#Complete overview
data_sloth
