# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: Main analyses on LPP with exploratory time window and EEG channels (ANOVAs, pairwise comparisons, plotting)

###################
### preparing R ###
###################

# random number generator seed (necessary for Bayesian ANOVAs)
# created with [set.seed(NULL)] and [sample(2^31 - 1, 1)] 
rngSeed <- 814677222

# loading required packages
library(tidyr) # ver. 1.3.0
library(psych) # ver. 2.3.9
library(effectsize) # ver. 0.8.6
library(ez) # ver. 4.4-0
library(BayesFactor) # ver. 0.9.12-4.5
library(bayestestR) # ver. 0.13.1
library(ggplot2) # ver. 3.4.2
library(scico) # ver. 1.5.0
library(scales) #  ver. 1.2.1
library(ggpubr) #  ver. 0.6.0
library(flextable) # ver. 0.9.4
library(stringr) #  ver. 1.5.0
library(here) # ver. 1.0.1
library(ggbeeswarm) # ver. 0.7.2
library(eegUtils) # ver. 0.7.0
library(OneR) # ver. 2.2



###################################
### Setting time window for LPP ###
###################################

# set by user
sRate <- 1024 # sampling rate of EEG data
startSeg <- -200 # time of first data point relative to CS
TWOI <- c(400,1000) # Time Window Of Interest (in ms)
chanNames <- c("Pz", "POz", "P1", "P2", "P3", "P4", "PO3", "PO4")

# computed
# Sample Window Of Interest
SWOI <- c(round((TWOI[1]-startSeg)*sRate/1000),round((TWOI[2]-startSeg)*sRate/1000)) 



########################
### data preparation ###
########################

# Rating data contains group membership
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# Load channel locations and transform from 3D theta + radius into 2D x & y
# for plotting purposes
loadname <- paste0(pathname,"/channelLocations/chanLocs_biosemi64.txt")
chanLocs <- read.csv(loadname, sep = ";")
chanLocs$thetaRadian <- pi/180*chanLocs$theta
chanLocs$x <- chanLocs$radius*sin(chanLocs$thetaRadian)*200
chanLocs$y <- chanLocs$radius*cos(chanLocs$thetaRadian)*200

chanInd <- match(chanNames, chanLocs$name)

# read LPP data
# create emtpy matrices, rows = participants, columns = sample points (entire ERP)
tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[1], "_IMAKON03_Average_Akquges_EKP_P9P10_S51.dat"), 
                     header = FALSE, sep = " ")
erpMatAV <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
erpMatNEU <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
erpMatMIN <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
topoMatAV <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[2])
topoMatNEU <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[2])
topoMatMIN <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[2])

# read single-subject data and extract Pz
for (partI in 1:length(importRatings$partCode)) {
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S51.dat"), 
                       header = FALSE, sep = " ")
  erpMatAV[partI,] <- rowMeans(tempData[,chanInd])
  topoMatAV[partI,] <- colMeans(tempData[SWOI[1]:SWOI[2],])
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S52.dat"), 
                       header = FALSE, sep = " ")
  erpMatNEU[partI,] <- rowMeans(tempData[,chanInd])
  topoMatNEU[partI,] <- colMeans(tempData[SWOI[1]:SWOI[2],])
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S53.dat"), 
                       header = FALSE, sep = " ")
  erpMatMIN[partI,] <- rowMeans(tempData[,chanInd])
  topoMatMIN[partI,] <- colMeans(tempData[SWOI[1]:SWOI[2],])
}

# create LPP data frame; LPP values are means from 400 to 1000 s post-CS
dataLPP <- data.frame(
  partCode = factor(importRatings$partCode),
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_allTr = rowMeans(erpMatAV[,SWOI[1]:SWOI[2]]),
  Neu_allTr = rowMeans(erpMatNEU[,SWOI[1]:SWOI[2]]),
  Min_allTr = rowMeans(erpMatMIN[,SWOI[1]:SWOI[2]])
)  

# transform into long format
dataLPPLong <- gather(data = dataLPP, key = "cond", value = "LPP", Av_allTr:Min_allTr)
dataLPPLong <- separate(data = dataLPPLong, col = cond, into = c("CS","time"), sep = "_")
dataLPPLong$CS <- factor(dataLPPLong$CS, levels = c("Av","Neu","Min"))
dataLPPLong$time <- factor(dataLPPLong$time)



########################################
### Across groups - primary analyses ###
########################################

# descriptive statistics for LPP ratings across conditioning groups
describe(dataLPP)

# frequentist ANOVA on LPP across conditioning groups
anovaLPP <- ezANOVA(
  data = dataLPPLong[dataLPPLong$time == "allTr",],
  dv = LPP,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaLPP$ANOVA$pEtaSq <- c(anovaLPP$ANOVA$SSn[1] /
                                (anovaLPP$ANOVA$SSd[1]+anovaLPP$ANOVA$SSn[1]),
                              anovaLPP$ANOVA$SSn[2] /
                                (anovaLPP$ANOVA$SSd[2]+anovaLPP$ANOVA$SSn[2]),
                              anovaLPP$ANOVA$SSn[3] /
                                (anovaLPP$ANOVA$SSd[3]+anovaLPP$ANOVA$SSn[3]),
                              anovaLPP$ANOVA$SSn[4] /
                                (anovaLPP$ANOVA$SSd[4]+anovaLPP$ANOVA$SSn[4])
); print(anovaLPP)

# bayesian ANOVA on LPP across conditioning groups
set.seed(rngSeed); anovaBFLPP <- anovaBF(
  formula = LPP ~ usGroup*CS + partInd,
  data = dataLPPLong[dataLPPLong$time == "allTr",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFLPP)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFLPP)

# quick graph of group x CS ANOVA on LPP
ezPlot(
  data = dataLPPLong[dataLPPLong$time == "allTr",],
  dv = LPP,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on LPP across conditioning groups
# CS+av vs CS+neu
lppAcrossAvNeu_t <- t.test(x = dataLPP$Av_allTr,
                           y = dataLPP$Neu_allTr,
                           alternative = "greater", paired = TRUE) # one-sided
lppAcrossAvNeu_d <- cohens_d(x = dataLPP$Av_allTr,
                             y = dataLPP$Neu_allTr,
                             paired = TRUE)
lppAcrossAvNeu_BF <- ttestBF(x = dataLPP$Av_allTr,
                             y = dataLPP$Neu_allTr,
                             nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppAcrossAvMin_t <- t.test(x = dataLPP$Av_allTr,
                           y = dataLPP$Min_allTr,
                           alternative = "greater", paired = TRUE) # one-sided
lppAcrossAvMin_d <- cohens_d(x = dataLPP$Av_allTr,
                             y = dataLPP$Min_allTr,
                             paired = TRUE)
lppAcrossAvMin_BF <- ttestBF(x = dataLPP$Av_allTr,
                             y = dataLPP$Min_allTr,
                             nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppAcrossNeuMin_t <- t.test(x = dataLPP$Neu_allTr,
                            y = dataLPP$Min_allTr,
                            alternative = "two.sided", paired = TRUE) # two-sided
lppAcrossNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr,
                              y = dataLPP$Min_allTr,
                              paired = TRUE)
lppAcrossNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr,
                              y = dataLPP$Min_allTr,
                              nullInterval = NULL, paired = TRUE) # two-sided

# frequentist & bayesian t-tests on LPP (difference scores) between groups
# delta [CS+av - CS+neu]
lppBetweenAvNeu_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                              dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                              dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                            alternative = "two.sided", paired = FALSE) # two-sided
lppBetweenAvNeu_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                              y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                              paired = FALSE)
lppBetweenAvNeu_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                              y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                              nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
lppBetweenAvMin_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                              dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                              dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                            alternative = "two.sided", paired = FALSE) # two-sided
lppBetweenAvMin_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                              y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                              paired = FALSE)
lppBetweenAvMin_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                              y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                              nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
lppBetweenNeuMin_t <- t.test(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                             y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
lppBetweenNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                                 dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                               y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                                 dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                               paired = FALSE)
lppBetweenNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                                 dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                               y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                                 dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided



#####################################################
### Imagery-based conditioning - primary analyses ###
#####################################################

# descriptive statistics for LPP in imagery-based conditioning group
describe(dataLPP[dataLPP$usGroup == "ima",])

# frequentist ANOVA on LPP in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = LPP
anovaLPPIma <- ezANOVA(
  data = dataLPPLong[dataLPPLong$usGroup == "ima" & dataLPPLong$time == "allTr",],
  dv = LPP,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaLPPIma$ANOVA$pEtaSq <- 
  c(anovaLPPIma$ANOVA$SSn[1] / (anovaLPPIma$ANOVA$SSd[1]+anovaLPPIma$ANOVA$SSn[1]),
    anovaLPPIma$ANOVA$SSn[2] / (anovaLPPIma$ANOVA$SSd[2]+anovaLPPIma$ANOVA$SSn[2])
  ); print(anovaLPPIma)

# bayesian ANOVA on LPP in imagery-based conditioning group
set.seed(rngSeed); anovaBFLPPIma <- anovaBF(
  formula = LPP ~ CS + partInd,
  data = dataLPPLong[dataLPPLong$usGroup == "ima" & dataLPPLong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFLPPIma)

# frequentist & bayesian t-tests on LPP in imagery-based conditioning group
# CS+av vs CS+neu
lppImaAvNeu_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                        y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
lppImaAvNeu_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                          y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                          paired = TRUE)
lppImaAvNeu_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                          y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppImaAvMin_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                        y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
lppImaAvMin_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                          y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                          paired = TRUE)
lppImaAvMin_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "ima"],
                          y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppImaNeuMin_t <- t.test(x = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                         y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
lppImaNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           paired = TRUE)
lppImaNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided



#################################################
### Classical conditioning - primary analyses ###
#################################################

# descriptive statistics for LPP in classical conditioning group
describe(dataLPP[dataLPP$usGroup == "real",])

# frequentist ANOVA on LPP in classical conditioning group, including p. eta^2
# IV = CS; DV = LPP
anovaLPPReal <- ezANOVA(
  data = dataLPPLong[dataLPPLong$usGroup == "real" & dataLPPLong$time == "allTr",],
  dv = LPP,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaLPPReal$ANOVA$pEtaSq <- 
  c(anovaLPPReal$ANOVA$SSn[1] / (anovaLPPReal$ANOVA$SSd[1]+anovaLPPReal$ANOVA$SSn[1]),
    anovaLPPReal$ANOVA$SSn[2] / (anovaLPPReal$ANOVA$SSd[2]+anovaLPPReal$ANOVA$SSn[2])
  ); print(anovaLPPReal)

# bayesian ANOVA on LPP in classical conditioning group
set.seed(rngSeed); anovaBFLPPReal <- anovaBF(
  formula = LPP ~ CS + partInd,
  data = dataLPPLong[dataLPPLong$usGroup == "real" & dataLPPLong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFLPPReal)

# frequentist & bayesian t-tests on LPP in classical conditioning group
# CS+av vs CS+neu
lppRealAvNeu_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                         y = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
lppRealAvNeu_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                           paired = TRUE)
lppRealAvNeu_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppRealAvMin_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                         y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
lppRealAvMin_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                           paired = TRUE)
lppRealAvMin_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppRealNeuMin_t <- t.test(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                          y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                          alternative = "two.sided", paired = TRUE) # two-sided
lppRealNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            paired = TRUE)
lppRealNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # two-sided



#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("across groups: CS+av vs CS+neu", "across groups: CS+av vs CS-", "across groups: CSneu vs CS-",
                     "imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-",
                     "between groups: delta CS+av / CS+neu", "between groups: delta CS+av / CS-", "between groups: delta CSneu / CS-"), 3),
  t = c(lppAcrossAvNeu_t$statistic, lppAcrossAvMin_t$statistic, lppAcrossNeuMin_t$statistic,
        lppImaAvNeu_t$statistic, lppImaAvMin_t$statistic, lppImaNeuMin_t$statistic,
        lppRealAvNeu_t$statistic, lppRealAvMin_t$statistic, lppRealNeuMin_t$statistic,
        lppBetweenAvNeu_t$statistic, lppBetweenAvMin_t$statistic, lppBetweenNeuMin_t$statistic),
  df = c(lppAcrossAvNeu_t$parameter, lppAcrossAvMin_t$parameter, lppAcrossNeuMin_t$parameter,
         lppImaAvNeu_t$parameter, lppImaAvMin_t$parameter, lppImaNeuMin_t$parameter,
         lppRealAvNeu_t$parameter, lppRealAvMin_t$parameter, lppRealNeuMin_t$parameter,
         lppBetweenAvNeu_t$parameter, lppBetweenAvMin_t$parameter, lppBetweenNeuMin_t$parameter), 
  p = c(lppAcrossAvNeu_t$p.value, lppAcrossAvMin_t$p.value, lppAcrossNeuMin_t$p.value,
        lppImaAvNeu_t$p.value, lppImaAvMin_t$p.value, lppImaNeuMin_t$p.value,
        lppRealAvNeu_t$p.value, lppRealAvMin_t$p.value, lppRealNeuMin_t$p.value,
        lppBetweenAvNeu_t$p.value*3, lppBetweenAvMin_t$p.value*3, lppBetweenNeuMin_t$p.value*3),  # Bonferroni
  d = c(lppAcrossAvNeu_d$Cohens_d, lppAcrossAvMin_d$Cohens_d, lppAcrossNeuMin_d$Cohens_d,
        lppImaAvNeu_d$Cohens_d, lppImaAvMin_d$Cohens_d, lppImaNeuMin_d$Cohens_d,
        lppRealAvNeu_d$Cohens_d, lppRealAvMin_d$Cohens_d, lppRealNeuMin_d$Cohens_d,
        lppBetweenAvNeu_d$Cohens_d, lppBetweenAvMin_d$Cohens_d, lppBetweenNeuMin_d$Cohens_d),
  BF = c(exp(lppAcrossAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppAcrossAvMin_BF@bayesFactor[["bf"]][1]), exp(lppAcrossNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(lppImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppImaAvMin_BF@bayesFactor[["bf"]][1]), exp(lppImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(lppRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppRealAvMin_BF@bayesFactor[["bf"]][1]), exp(lppRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(lppBetweenAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppBetweenAvMin_BF@bayesFactor[["bf"]][1]), exp(lppBetweenNeuMin_BF@bayesFactor[["bf"]][1]))
)
# round the numbers
tableData$t <-round(tableData$t, 2)
tableData$df <-round(tableData$df, 0)
tableData$p[tableData$p < .001] <- 0
tableData$p[tableData$p > 1] <- 1
tableData$p <-round(tableData$p, 3)
tableData$p <- as.character(tableData$p)
tableData$p[tableData$p == "0"] <- "< .001"
tableData$p <- str_replace(tableData$p, "0\\.", "\\.")
tableData$d <-round(tableData$d, 2)
tableData$BF <- format(tableData$BF, digits = 3)

tableLPP <- flextable(tableData[1:12,])
tableLPP <- add_header_lines(tableLPP, top = TRUE, values = "LPP")
tableLPP <- align(tableLPP, align = "center")

save_as_docx(tableLPP, path = paste0(pathname, "/supplement/04s_tableLPP_exploratory_raw.docx"))



####################
### Plotting LPP ###
####################

csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")
lppAv4GA <- erpMatAV
  lppAv4GAima <- erpMatAV[importRatings$group == "ima",]
  lppAv4GAreal <- erpMatAV[importRatings$group == "real",]
lppNeu4GA <- erpMatNEU
  lppNeu4GAima <- erpMatNEU[importRatings$group == "ima",]
  lppNeu4GAreal <- erpMatNEU[importRatings$group == "real",]
lppMin4GA <- erpMatMIN
  lppMin4GAima <- erpMatMIN[importRatings$group == "ima",]
  lppMin4GAreal <- erpMatMIN[importRatings$group == "real",]
  
lppGA <- data.frame(
  csType = factor(c(rep(1,dim(lppAv4GA)[2]), rep(2,dim(lppAv4GA)[2]), rep(3,dim(lppAv4GA)[2])), labels = csLabels),
  time = rep(seq(startSeg,startSeg+(dim(lppAv4GA)[2]-1)*1000/sRate,1000/sRate),3),
  LPP = c(colMeans(lppAv4GA),
          colMeans(lppNeu4GA),
          colMeans(lppMin4GA)),
  row.names = NULL)
lppGAima <- data.frame(
  csType = factor(c(rep(1,dim(lppAv4GAima)[2]), rep(2,dim(lppAv4GAima)[2]), rep(3,dim(lppAv4GAima)[2])), labels = csLabels),
  time = rep(seq(startSeg,startSeg+(dim(lppAv4GAima)[2]-1)*1000/sRate,1000/sRate),3),
  LPP = c(colMeans(lppAv4GAima),
          colMeans(lppNeu4GAima),
          colMeans(lppMin4GAima)),
  row.names = NULL)
lppGAreal <- data.frame(
  csType = factor(c(rep(1,dim(lppAv4GAreal)[2]), rep(2,dim(lppAv4GAreal)[2]), rep(3,dim(lppAv4GAreal)[2])), labels = csLabels),
  time = rep(seq(startSeg,startSeg+(dim(lppAv4GAreal)[2]-1)*1000/sRate,1000/sRate),3),
  LPP = c(colMeans(lppAv4GAreal),
          colMeans(lppNeu4GAreal),
          colMeans(lppMin4GAreal)),
  row.names = NULL)




# compute difference values
dataLPPdiff <- dataLPPLong
dataLPPdiffAvNeu <- aggregate(LPP ~ partCode + partInd + usGroup, data = dataLPPLong[dataLPPLong$CS == "Av" | dataLPPLong$CS == "Neu",], FUN = "diff")
dataLPPdiffAvMin <- aggregate(LPP ~ partCode + partInd + usGroup, data = dataLPPLong[dataLPPLong$CS == "Av" | dataLPPLong$CS == "Min",], FUN = "diff")
#dataLPPdiffNeuMin <- aggregate(LPP ~ partCode + partInd + usGroup, data = dataLPPLong[dataLPPLong$CS == "Neu" | dataLPPLong$CS == "Min",], FUN = "diff")

#dataLPPdiff <- rbind(dataLPPdiffAvNeu, dataLPPdiffAvMin, dataLPPdiffNeuMin)
dataLPPdiff <- rbind(dataLPPdiffAvNeu, dataLPPdiffAvMin)
dataLPPdiff$LPP <- -dataLPPdiff$LPP

#dataLPPdiff$comp <- factor(c(rep("AVvsNEU", 48), rep("AVvsMIN", 48), rep("NeuvsMIN", 48)), levels = c("AVvsNEU", "AVvsMIN", "NEUvsMIN"))
dataLPPdiff$comp <- factor(c(rep("AVvsNEU", 48), rep("AVvsMIN", 48)), levels = c("AVvsNEU", "AVvsMIN"))

#compLabels = c(expression(paste("CS+"[av], " - CS+"[neu])), expression(paste("\n","CS+"[av], " - CS-")), expression(paste("CS+"[neu], " - CS-")))
compLabels = c(expression(paste("CS+"[av], " - CS+"[neu])), expression(paste("\n","CS+"[av], " - CS-")))


# settings for plotting
lineSize = 1
yMin = -5
yMax = 12
plotFS <- 8
showSig <- TRUE


# ERP at PZ for imagery-based conditioning group
graphLPPima <- ggplot(data = lppGAima, aes(x = time, y = LPP, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = TWOI[1], xmax = TWOI[2], ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), linewidth = lineSize) + 
  scale_x_continuous(breaks = seq(-200,1000,200)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7), labels = csLabels) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Imagery-Based Conditioning", x = "Time (ms)", y = expression(paste("ERP amplitude at Pz ", (µV))), fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = "none") +
  theme(
    legend.position = c(.85,.2),
    legend.direction = "vertical",
    legend.text.align = 0,
    legend.text = element_text(size = plotFS-2),
    legend.key.size = unit(.5, "lines"),
    legend.background = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
    axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

# ERP at PZ for classical conditioning group
graphLPPreal <- ggplot(data = lppGAreal, aes(x = time, y = LPP, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = TWOI[1], xmax = TWOI[2], ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), linewidth = lineSize) + 
  scale_x_continuous(breaks = seq(-200,1000,200)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7), labels = csLabels) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Classical Conditioning", x = "Time (ms)", y = expression(paste("ERP amplitude at Pz ", (µV))), fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = "none") +
  theme(
    legend.position = c(.85,.2),
    legend.direction = "vertical",
    legend.text.align = 0,
    legend.text = element_text(size = plotFS-2),
    legend.key.size = unit(.5, "lines"),
    legend.background = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
    axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

# difference plots for imagery-based conditioning group
graphLPPdiffIma <- ggplot(data = dataLPPdiff[dataLPPdiff$usGroup == "ima",], aes(x = comp, y = LPP)) +
  theme_classic() +
  geom_violin(alpha = .2, color = NA, fill = "gray50", bw = 1, width = 0.75) +
  geom_quasirandom(size = .5, width = .10, color = "gray70") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .10) +
  stat_summary(aes(x = comp, y = LPP), fun = "mean", geom = "crossbar", linewidth = .3, width = .25, color = "black") +
  geom_hline(yintercept = 0) +
  scale_x_discrete(name = " ", labels = compLabels, position = "bottom") +
  scale_y_continuous(name = "Difference in mean LPP amplitude (400-1000 ms)", limits = c(-8,16), breaks = seq(-5,15,5), expand = c(0,0)) +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS-2, color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

# difference plots for classical conditioning group
graphLPPdiffReal <- ggplot(data = dataLPPdiff[dataLPPdiff$usGroup == "real",], aes(x = comp, y = LPP)) +
  theme_classic() +
  geom_violin(alpha = .2, color = NA, fill = "gray50", bw = 1, width = 0.75) +
  geom_quasirandom(size = .5, width = .10, color = "gray70") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .10) +
  stat_summary(aes(x = comp, y = LPP), fun = "mean", geom = "crossbar", linewidth = .3, width = .25, color = "black") +
  geom_hline(yintercept = 0) +
  scale_x_discrete(name = " ", labels = compLabels, position = "bottom") +
  scale_y_continuous(name = "Difference in mean LPP amplitude (400-1000 ms)", limits = c(-8,16), breaks = seq(-5,15,5), expand = c(0,0)) +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS-2, color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

if (showSig == TRUE){
  graphLPPdiffIma <- graphLPPdiffIma +
    geom_text(aes(label = "†", x = 1, y = 15.0), size = plotFS/2, color = "darkred", family = "Helvetica") +
    geom_text(aes(label = "***", x = 2, y = 14.25), size = plotFS/1.5, color = "darkred")
  graphLPPdiffReal <- graphLPPdiffReal +  
    geom_text(aes(label = "**", x = 1, y = 14.25), size = plotFS/1.5, color = "darkred") +
    geom_text(aes(label = "†", x = 2, y = 15.0), size = plotFS/2, color = "darkred", family = "Helvetica")
}



# Average individual amplitudes across whole sample
topoAvgAV <- colMeans(topoMatAV)
topoAvgNEU <- colMeans(topoMatNEU)
topoAvgMIN <- colMeans(topoMatMIN)

topoAvgAVima <- colMeans(topoMatAV[importRatings$group == "ima",])
topoAvgNEUima <- colMeans(topoMatNEU[importRatings$group == "ima",])
topoAvgMINima <- colMeans(topoMatMIN[importRatings$group == "ima",])

topoAvgAVreal <- colMeans(topoMatAV[importRatings$group == "real",])
topoAvgNEUreal <- colMeans(topoMatNEU[importRatings$group == "real",])
topoAvgMINreal <- colMeans(topoMatMIN[importRatings$group == "real",])

# number of electrodes
nrChans = length(topoAvgAV)

# Create data frame with factors lab, driving frequency & modulation function, 
# electrode name, x & y coordinates for plot & LPP amplitude

# across groups
dfTopos <- data.frame(
  electrode = chanLocs$name,
  x = chanLocs$x,
  y = chanLocs$y,
  csp_av = topoAvgAV,
  csp_neu = topoAvgNEU,
  csm = topoAvgMIN
)

# imagery-based conditioning group
dfToposIma <- data.frame(
  electrode = chanLocs$name,
  x = chanLocs$x,
  y = chanLocs$y,
  csp_av = topoAvgAVima,
  csp_neu = topoAvgNEUima,
  csm = topoAvgMINima
)

# classical conditioning group
dfToposReal <- data.frame(
  electrode = chanLocs$name,
  x = chanLocs$x,
  y = chanLocs$y,
  csp_av = topoAvgAVreal,
  csp_neu = topoAvgNEUreal,
  csm = topoAvgMINreal
)

# compute difference values comparing CS+av to the mean of CS+neu and CS-
dfTopos$contrast <- 1*dfTopos$csp_av - 0.5*dfTopos$csp_neu - 0.5*dfTopos$csm
dfToposIma$contrast <- 1*dfToposIma$csp_av - 0.5*dfToposIma$csp_neu - 0.5*dfToposIma$csm
dfToposReal$contrast <- 1*dfToposReal$csp_av - 0.5*dfToposReal$csp_neu - 0.5*dfToposReal$csm



# settings for topography plots
topoRes <- 400
chanCol <- "black"
nrColors <- 8

minLim <- min(dfTopos$contrast)
maxLim <- max(dfTopos$contrast)
absLim <- max(abs(c(minLim, maxLim)))

### create topography plots

# across groups
topoAll <- topoplot(data = dfTopos,
                    contour = FALSE, interp_limit = "head", highlights = chanNames,
                    grid_res = topoRes, quantity = "contrast", scaling = .5)
# transform color scale into discrete steps
topoAll$data$fill <- as.numeric(bin(data = topoAll$data$fill, nbins = nrColors))


# imagery-based conditioning group  
topoIma <- topoplot(data = dfToposIma,
                    contour = FALSE, interp_limit = "head", highlights = chanNames,
                    grid_res = topoRes, quantity = "contrast", scaling = .33, 
                    limits = c(-absLim, absLim)) +
  labs(title = expression(paste("CS+"[av], " vs [CS+"[neu], " & CS-]"))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = plotFS),
        plot.title = element_text(hjust = 0.5, size = plotFS))
# format the color bar
topoIma$guides$fill$barwidth <- unit(6, "lines")
topoIma$guides$fill$barheight <- unit(.5, "lines")
topoIma$guides$fill$title <- expression(paste(Delta, " amplitude ", (µV)))
topoIma$guides$fill$title.theme$size <- plotFS
topoIma$guides$fill$title.theme$angle <- 0
topoIma$guides$fill$title.position <- "bottom"
topoIma$guides$fill$title.hjust <- 0.5
topoIma$guides$fill$nbin <- nrColors


# classical conditioning group  
topoReal <- topoplot(data = dfToposReal,
                     contour = FALSE, interp_limit = "head", highlights = chanNames,
                     grid_res = topoRes, quantity = "contrast", scaling = .33,
                     limits = c(-absLim, absLim)) +
  labs(title = expression(paste("CS+"[av], " vs [CS+"[neu], " & CS-]"))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = plotFS, family = "Arial"),
        plot.title = element_text(hjust = 0.5, size = plotFS, family = "Helvetica"))
# format the color bar
topoReal$guides$fill$barwidth <- unit(6, "lines")
topoReal$guides$fill$barheight <- unit(.5, "lines")
topoReal$guides$fill$title <- expression(paste(Delta, " amplitude ", (µV)))
topoReal$guides$fill$title.theme$size <- plotFS
topoReal$guides$fill$title.theme$angle <- 0
topoReal$guides$fill$title.position <- "bottom"
topoReal$guides$fill$title.hjust <- 0.5
topoReal$guides$fill$nbin <- nrColors

# use the same, discrete color scales for the two groups
topoFillVec <- c(topoIma$data$fill, topoReal$data$fill)
totMin <- min(topoFillVec)
totMax <- max(topoFillVec)

topoFillVec <- as.numeric(bin(data = topoFillVec, nbins = nrColors))

topoFillVec <- (topoFillVec-1) / (nrColors-1)
topoFillVec <- topoFillVec*(totMax-totMin) + totMin

topoIma$data$fill <- topoFillVec[1 : (length(topoFillVec)/2)]
topoReal$data$fill <- topoFillVec[(length(topoFillVec)/2+1) : length(topoFillVec)]



### combining graphs into one figure
# add margins
graphLPPima <- graphLPPima + theme(plot.margin = unit(c(10,5,5,5), "mm"))
graphLPPreal <- graphLPPreal + theme(plot.margin = unit(c(10,5,5,5), "mm"))
graphLPPdiffIma <- graphLPPdiffIma + theme(plot.margin = unit(c(10,2.5,5,7.5), "mm"))
graphLPPdiffReal <- graphLPPdiffReal + theme(plot.margin = unit(c(10,2.5,5,7.5), "mm"))
topoIma <- topoIma + theme(plot.margin = unit(c(10,5,2.6,5), "mm"))
topoReal <- topoReal + theme(plot.margin = unit(c(10,5,2.6,5), "mm"))

# create panels and merge them
graphLPProw1 <- ggarrange(graphLPPima, graphLPPdiffIma, topoIma,
                          ncol = 3, nrow = 1, 
                          widths = c(3,2,2))
graphLPProw2 <- ggarrange(graphLPPreal, graphLPPdiffReal, topoReal,
                          ncol = 3, nrow = 1, 
                          widths = c(3,2,2))
graphLPP <- ggarrange(graphLPProw1,graphLPProw2,
                      ncol = 1, nrow = 2,
                      labels = c("A) Imagery-Based Conditioning", "B) Classical Conditioning"),
                      hjust = -.05
)
# plot
graphLPP

# saving it
ggsave(filename = paste0(pathname, "/supplement/Figure5b_timeCourse_diffPlot_LPP_exploratory.png"),
       plot = graphLPP,
       width = 200,
       height = 150,
       units = "mm",
       dpi = 1200
)

ggsave(filename = paste0(pathname, "/supplement/Figure5b_timeCourse_diffPlot_LPP_exploratory.pdf"),
       plot = graphLPP,
       width = 200,
       height = 150,
       units = "mm"
)
