# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: trial numbers for LPP by conditioning group

################
### packages ###
################

# loading required packages
library(psych) # ver. 2.3.9
library(here) # ver. 1.0.1
library(tidyr)
library(effectsize)
library(ez)


##################################
### loading and reporting data ###
##################################

# load rating data from text file
pathname <- here()
dfTrials <- read.csv(paste0(pathname, "/experimentData/imaCond3_trialNumbers.txt"), sep=",")
dfTrialsLong <- pivot_longer(data = dfTrials, cols = c(csPlusAversive,csPlusNeutral,csMinus), names_to = "csType", values_to = "trials") 
dfTrialsLong$partCode <- factor(dfTrialsLong$partCode)
dfTrialsLong$group <- factor(dfTrialsLong$group, levels = c("ima","real"))
dfTrialsLong$csType <- factor(dfTrialsLong$csType, levels = c("csPlusAversive","csPlusNeutral","csMinus"))

# descriptive statistics for number of trials across and separate for conditioning groups
describe(dfTrials)
describeBy(dfTrials, group = dfTrials$group)

# ANOVA
anovaTrials <- ezANOVA(data = dfTrialsLong,
                       dv = trials,
                       wid = partCode,
                       within = csType,
                       between = group,
                       type = 3,
                       detailed = TRUE); anovaTrials$ANOVA$pEtaSq <- 
  c(anovaTrials$ANOVA$SSn[1] / (anovaTrials$ANOVA$SSd[1]+anovaTrials$ANOVA$SSn[1]),
    anovaTrials$ANOVA$SSn[2] / (anovaTrials$ANOVA$SSd[2]+anovaTrials$ANOVA$SSn[2]),
    anovaTrials$ANOVA$SSn[3] / (anovaTrials$ANOVA$SSd[3]+anovaTrials$ANOVA$SSn[3]),
    anovaTrials$ANOVA$SSn[4] / (anovaTrials$ANOVA$SSd[4]+anovaTrials$ANOVA$SSn[4])
  ); print(anovaTrials)

anovaTrialsIma <- ezANOVA(data = dfTrialsLong[dfTrialsLong$group == "ima",],
                          dv = trials,
                          wid = partCode,
                          within = csType,
                          type = 3,
                          detailed = TRUE); anovaTrialsIma$ANOVA$pEtaSq <- 
  c(anovaTrialsIma$ANOVA$SSn[1] / (anovaTrialsIma$ANOVA$SSd[1]+anovaTrialsIma$ANOVA$SSn[1]),
    anovaTrialsIma$ANOVA$SSn[2] / (anovaTrialsIma$ANOVA$SSd[2]+anovaTrialsIma$ANOVA$SSn[2])
  ); print(anovaTrialsIma)

anovaTrialsReal <- ezANOVA(data = dfTrialsLong[dfTrialsLong$group == "real",],
                           dv = trials,
                           wid = partCode,
                           within = csType,
                           type = 3,
                           detailed = TRUE); anovaTrialsReal$ANOVA$pEtaSq <- 
  c(anovaTrialsReal$ANOVA$SSn[1] / (anovaTrialsReal$ANOVA$SSd[1]+anovaTrialsReal$ANOVA$SSn[1]),
    anovaTrialsReal$ANOVA$SSn[2] / (anovaTrialsReal$ANOVA$SSd[2]+anovaTrialsReal$ANOVA$SSn[2])
  ); print(anovaTrialsReal)
