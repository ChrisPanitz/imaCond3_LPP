# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Feb 2022
# --- content: LPP analyses of imagery-based conditioning data in Panitz & Mueller (2022)

###################
### preparing R ###
###################

# random number generator seed (necessary for Bayesian ANOVAs)
# created with [set.seed(NULL)] and [sample(2^31 - 1, 1)] 
rngSeed <- 814677222

# (installing and) loading required packages
# will install newest package version, not necessarily the version originally used!
if(!is.element("tidyr",installed.packages()[,1])) {install.packages("tidyr")}
library(tidyr) # ver. 1.1.2
if(!is.element("psych",installed.packages()[,1])) {install.packages("psych")}
library(psych) # ver. 2.0.9
if(!is.element("effectsize",installed.packages()[,1])) {install.packages("effectsize")}
library(effectsize) # ver. 2.0.9
if(!is.element("ez",installed.packages()[,1])) {install.packages("ez")}
library(ez) # ver. 4.4-0
if(!is.element("BayesFactor",installed.packages()[,1])) {install.packages("BayesFactor")}
library(BayesFactor) # ver. 2.0.9
if(!is.element("bayestestR",installed.packages()[,1])) {install.packages("BayesFactor")}
library(bayestestR) #
if(!is.element("ggplot2",installed.packages()[,1])) {install.packages("ggplot2")}
library(ggplot2) # ver. 3.3.2
if(!is.element("scico",installed.packages()[,1])) {install.packages("scico")}
library(scico) # ver. 1.2.0
if(!is.element("scales",installed.packages()[,1])) {install.packages("scales")}
library(scales) # 
if(!is.element("ggpubr",installed.packages()[,1])) {install.packages("ggpubr")}
library(ggpubr) # 
if(!is.element("flextable",installed.packages()[,1])) {install.packages("flextable")}
library(flextable) # 
if(!is.element("stringr",installed.packages()[,1])) {install.packages("stringr")}
library(stringr) # 
if(!is.element("here",installed.packages()[,1])) {install.packages("here")}
library(here) # 



###################################
### Setting time window for LPP ###
###################################

# set by user
sRate <- 1024 # sampling rate of EEG data
startSeg <- -200 # time of first data point relative to CS
TWOI <- c(300,700) # Time Window Of Interest (in ms)
chanInd <- 31 # 31 = index of Pz in EEG array

# computed
SWOI <- c(round(TWOI[1]-startSeg*sRate/1000),round(TWOI[2]-startSeg*sRate/1000)) # Sample Window Of Interest



########################
### data preparation ###
########################

# Rating data contains group membership
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname,"/experimentData/imaCond3_demographicsAndRatings.txt"), sep = ",")


# read LPP data
# create emtpy matrices, rows = participants, columns = sample points (entire ERP)
tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[1], "_IMAKON03_Average_Akqu1_EKP_P9P10_S51.dat"), 
                     header = FALSE, sep = " ")
lppMatAV1 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatNEU1 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatMIN1 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatAV2 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatNEU2 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatMIN2 <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])

# read single-subject data and extract Pz
for (partI in 1:length(importRatings$partCode)) {
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu1_EKP_P9P10_S51.dat"), 
                       header = FALSE, sep = " ")
  lppMatAV1[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu1_EKP_P9P10_S52.dat"), 
                       header = FALSE, sep = " ")
  lppMatNEU1[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu1_EKP_P9P10_S53.dat"), 
                       header = FALSE, sep = " ")
  lppMatMIN1[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu2_EKP_P9P10_S51.dat"), 
                       header = FALSE, sep = " ")
  lppMatAV2[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu2_EKP_P9P10_S52.dat"), 
                       header = FALSE, sep = " ")
  lppMatNEU2[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akqu2_EKP_P9P10_S53.dat"), 
                       header = FALSE, sep = " ")
  lppMatMIN2[partI,] <- tempData[,chanInd]
}

# create LPP data frame; LPP values are means from 300 to 700 s post-CS
dataLPP <- data.frame(
  partCode = factor(importRatings$partCode),
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
  Av_1stBl = rowMeans(lppMatAV1[,SWOI[1]:SWOI[2]]),
  Neu_1stBl = rowMeans(lppMatNEU1[,SWOI[1]:SWOI[2]]),
  Min_1stBl = rowMeans(lppMatMIN1[,SWOI[1]:SWOI[2]]),
  Av_2ndBl = rowMeans(lppMatAV2[,SWOI[1]:SWOI[2]]),
  Neu_2ndBl = rowMeans(lppMatNEU2[,SWOI[1]:SWOI[2]]),
  Min_2ndBl = rowMeans(lppMatMIN2[,SWOI[1]:SWOI[2]])  
)  

# transform into long format
dataLPPLong <- gather(data = dataLPP, key = "cond", value = "LPP", Av_1stBl:Min_2ndBl)
dataLPPLong <- separate(data = dataLPPLong, col = cond, into = c("CS","time"), sep = "_")
dataLPPLong$CS <- factor(dataLPPLong$CS)
dataLPPLong$time <- factor(dataLPPLong$time)



############################################################################
### Imagery-based conditioning - supplementary analyses with time factor ###
############################################################################

# descriptive statistics for LPP in imagery-based conditioning group
describe(dataLPP[dataLPP$usGroup == "ima",])

# frequentist CS x Time ANOVA on LPP in imagery-based conditioning group, including p. eta^2
# IV = CS, Time; DV = LPP
anovaLPPIma <- ezANOVA(
  data = dataLPPLong[dataLPPLong$usGroup == "ima",],
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaLPPIma$ANOVA$pEtaSq <- 
  c(anovaLPPIma$ANOVA$SSn[1] / (anovaLPPIma$ANOVA$SSd[1]+anovaLPPIma$ANOVA$SSn[1]),
    anovaLPPIma$ANOVA$SSn[2] / (anovaLPPIma$ANOVA$SSd[2]+anovaLPPIma$ANOVA$SSn[2]),
    anovaLPPIma$ANOVA$SSn[3] / (anovaLPPIma$ANOVA$SSd[3]+anovaLPPIma$ANOVA$SSn[3]),
    anovaLPPIma$ANOVA$SSn[4] / (anovaLPPIma$ANOVA$SSd[4]+anovaLPPIma$ANOVA$SSn[4])
  ); print(anovaLPPIma)
capture.output(print(anovaLPPIma), file = paste0(pathname, "/supplement/04s_LPP_timeFactor_ima_anovaFreq.doc"))

# bayesian CS x Time ANOVA on LPP in imagery-based conditioning group
set.seed(rngSeed); anovaBFLPPIma <- generalTestBF(
  formula = LPP ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataLPPLong[dataLPPLong$usGroup == "ima",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
); print(anovaBFLPPIma)
capture.output(print(anovaBFLPPIma), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_ima_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFLPPIma)
capture.output(bf_inclusion(anovaBFLPPIma), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_ima_BFinclusion.doc"))

# quick graph of CS Type x Time ANOVA for LPP in imagery-based conditioning group
plotLPPIma <- ezPlot(
  data = dataLPPLong[dataLPPLong$usGroup == "ima",],
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotLPPIma 
ggsave(plot = plotLPPIma, filename = paste0(pathname, "/supplement/04s_lpp_timeFactor_ima_plot.jpg"),
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on LPP in imagery-based conditioning group
### 1stBl
# CS+av vs CS+neu
lppImaAvNeu1stBl_t <- t.test(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                             y = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
lppImaAvNeu1stBl_d <- cohens_d(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                               paired = TRUE)
lppImaAvNeu1stBl_BF <- ttestBF(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppImaAvMin1stBl_t <- t.test(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                             y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
lppImaAvMin1stBl_d <- cohens_d(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                               paired = TRUE)
lppImaAvMin1stBl_BF <- ttestBF(x = dataLPP$Av_1stBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppImaNeuMin1stBl_t <- t.test(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                              y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
lppImaNeuMin1stBl_d <- cohens_d(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                                y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                                paired = TRUE)
lppImaNeuMin1stBl_BF <- ttestBF(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "ima"],
                                y = dataLPP$Min_1stBl[dataLPP$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
lppImaAvNeu2ndBl_t <- t.test(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                             y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
lppImaAvNeu2ndBl_d <- cohens_d(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                               paired = TRUE)
lppImaAvNeu2ndBl_BF <- ttestBF(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppImaAvMin2ndBl_t <- t.test(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                             y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
lppImaAvMin2ndBl_d <- cohens_d(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                               paired = TRUE)
lppImaAvMin2ndBl_BF <- ttestBF(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "ima"],
                               y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppImaNeuMin2ndBl_t <- t.test(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                              y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
lppImaNeuMin2ndBl_d <- cohens_d(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                                y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                                paired = TRUE)
lppImaNeuMin2ndBl_BF <- ttestBF(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "ima"],
                                y = dataLPP$Min_2ndBl[dataLPP$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided

tableLPPIma <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(lppImaAvNeu1stBl_t$statistic, lppImaAvMin1stBl_t$statistic, lppImaNeuMin1stBl_t$statistic,
        lppImaAvNeu2ndBl_t$statistic, lppImaAvMin2ndBl_t$statistic, lppImaNeuMin2ndBl_t$statistic),
  df = c(lppImaAvNeu1stBl_t$parameter, lppImaAvMin1stBl_t$parameter, lppImaNeuMin1stBl_t$parameter,
         lppImaAvNeu2ndBl_t$parameter, lppImaAvMin2ndBl_t$parameter, lppImaNeuMin2ndBl_t$parameter), 
  p = c(lppImaAvNeu1stBl_t$p.value, lppImaAvMin1stBl_t$p.value, lppImaNeuMin1stBl_t$p.value,
        lppImaAvNeu2ndBl_t$p.value, lppImaAvMin2ndBl_t$p.value, lppImaNeuMin2ndBl_t$p.value),
  d = c(lppImaAvNeu1stBl_d$Cohens_d, lppImaAvMin1stBl_d$Cohens_d, lppImaNeuMin1stBl_d$Cohens_d,
        lppImaAvNeu2ndBl_d$Cohens_d, lppImaAvMin2ndBl_d$Cohens_d, lppImaNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(lppImaAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(lppImaAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(lppImaNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(lppImaAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(lppImaAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(lppImaNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableLPPIma, file = paste0(pathname, "/supplement/04s_lpp_timeFactor_ima_tTable.doc"))



########################################################################
### Classical conditioning - supplementary analyses with time factor ###
########################################################################

# descriptive statistics for LPP in classical conditioning group
describe(dataLPP[dataLPP$usGroup == "real",])

# frequentist CS x Time ANOVA on LPP in classical conditioning group, including p. eta^2
# IV = CS, Time; DV = LPP
anovaLPPReal <- ezANOVA(
  data = dataLPPLong[dataLPPLong$usGroup == "real",],
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaLPPReal$ANOVA$pEtaSq <- 
  c(anovaLPPReal$ANOVA$SSn[1] / (anovaLPPReal$ANOVA$SSd[1]+anovaLPPReal$ANOVA$SSn[1]),
    anovaLPPReal$ANOVA$SSn[2] / (anovaLPPReal$ANOVA$SSd[2]+anovaLPPReal$ANOVA$SSn[2]),
    anovaLPPReal$ANOVA$SSn[3] / (anovaLPPReal$ANOVA$SSd[3]+anovaLPPReal$ANOVA$SSn[3]),
    anovaLPPReal$ANOVA$SSn[4] / (anovaLPPReal$ANOVA$SSd[4]+anovaLPPReal$ANOVA$SSn[4])
  ); print(anovaLPPReal)
capture.output(print(anovaLPPReal), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_real_anovaFreq.doc"))

# bayesian CS x Time ANOVA on LPP in classical conditioning group
set.seed(rngSeed); anovaBFLPPReal <- generalTestBF(
  formula = LPP ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataLPPLong[dataLPPLong$usGroup == "real",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
); print(anovaBFLPPReal)
capture.output(print(anovaBFLPPReal), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_real_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFLPPReal)
capture.output(bf_inclusion(anovaBFLPPReal), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_real_BFinclusion.doc"))

# quick graph of CS Type x Time ANOVA for LPP in imagery-based conditioning group
plotLPPReal <- ezPlot(
  data = dataLPPLong[dataLPPLong$usGroup == "real",],
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotLPPReal
ggsave(plot = plotLPPReal, filename = paste0(pathname, "/supplement/04s_lpp_timeFactor_real_plot.jpg"),
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on LPP in classical conditioning group
### 1stBl
# CS+av vs CS+neu
lppRealAvNeu1stBl_t <- t.test(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                              y = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
lppRealAvNeu1stBl_d <- cohens_d(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                                paired = TRUE)
lppRealAvNeu1stBl_BF <- ttestBF(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppRealAvMin1stBl_t <- t.test(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                              y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
lppRealAvMin1stBl_d <- cohens_d(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                                paired = TRUE)
lppRealAvMin1stBl_BF <- ttestBF(x = dataLPP$Av_1stBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppRealNeuMin1stBl_t <- t.test(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                               y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
lppRealNeuMin1stBl_d <- cohens_d(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                                 y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                                 paired = TRUE)
lppRealNeuMin1stBl_BF <- ttestBF(x = dataLPP$Neu_1stBl[dataLPP$usGroup == "real"],
                                 y = dataLPP$Min_1stBl[dataLPP$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
lppRealAvNeu2ndBl_t <- t.test(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                              y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
lppRealAvNeu2ndBl_d <- cohens_d(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                                paired = TRUE)
lppRealAvNeu2ndBl_BF <- ttestBF(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
lppRealAvMin2ndBl_t <- t.test(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                              y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
lppRealAvMin2ndBl_d <- cohens_d(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                                paired = TRUE)
lppRealAvMin2ndBl_BF <- ttestBF(x = dataLPP$Av_2ndBl[dataLPP$usGroup == "real"],
                                y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
lppRealNeuMin2ndBl_t <- t.test(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                               y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
lppRealNeuMin2ndBl_d <- cohens_d(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                                 y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                                 paired = TRUE)
lppRealNeuMin2ndBl_BF <- ttestBF(x = dataLPP$Neu_2ndBl[dataLPP$usGroup == "real"],
                                 y = dataLPP$Min_2ndBl[dataLPP$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided

tableLPPReal <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(lppRealAvNeu1stBl_t$statistic, lppRealAvMin1stBl_t$statistic, lppRealNeuMin1stBl_t$statistic,
        lppRealAvNeu2ndBl_t$statistic, lppRealAvMin2ndBl_t$statistic, lppRealNeuMin2ndBl_t$statistic),
  df = c(lppRealAvNeu1stBl_t$parameter, lppRealAvMin1stBl_t$parameter, lppRealNeuMin1stBl_t$parameter,
         lppRealAvNeu2ndBl_t$parameter, lppRealAvMin2ndBl_t$parameter, lppRealNeuMin2ndBl_t$parameter), 
  p = c(lppRealAvNeu1stBl_t$p.value, lppRealAvMin1stBl_t$p.value, lppRealNeuMin1stBl_t$p.value,
        lppRealAvNeu2ndBl_t$p.value, lppRealAvMin2ndBl_t$p.value, lppRealNeuMin2ndBl_t$p.value),
  d = c(lppRealAvNeu1stBl_d$Cohens_d, lppRealAvMin1stBl_d$Cohens_d, lppRealNeuMin1stBl_d$Cohens_d,
        lppRealAvNeu2ndBl_d$Cohens_d, lppRealAvMin2ndBl_d$Cohens_d, lppRealNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(lppRealAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(lppRealAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(lppRealNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(lppRealAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(lppRealAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(lppRealNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableLPPReal, file = paste0(pathname, "/supplement/04s_lpp_timeFactor_real_tTable.doc"))



###############################################################
### Across groups - supplementary analyses with time factor ###
###############################################################

# descriptive statistics for LPP ratings across conditioning groups
describe(dataLPP)

# frequentist Group x CS x Time ANOVA on LPP across conditioning groups
anovaLPP <- ezANOVA(
  data = dataLPPLong,
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaLPP$ANOVA$pEtaSq <- c(
  anovaLPP$ANOVA$SSn[1] / (anovaLPP$ANOVA$SSd[1]+anovaLPP$ANOVA$SSn[1]),
  anovaLPP$ANOVA$SSn[2] / (anovaLPP$ANOVA$SSd[2]+anovaLPP$ANOVA$SSn[2]),
  anovaLPP$ANOVA$SSn[3] / (anovaLPP$ANOVA$SSd[3]+anovaLPP$ANOVA$SSn[3]),
  anovaLPP$ANOVA$SSn[4] / (anovaLPP$ANOVA$SSd[4]+anovaLPP$ANOVA$SSn[4]),
  anovaLPP$ANOVA$SSn[5] / (anovaLPP$ANOVA$SSd[5]+anovaLPP$ANOVA$SSn[5]),
  anovaLPP$ANOVA$SSn[6] / (anovaLPP$ANOVA$SSd[6]+anovaLPP$ANOVA$SSn[6]),
  anovaLPP$ANOVA$SSn[7] / (anovaLPP$ANOVA$SSd[7]+anovaLPP$ANOVA$SSn[7]),
  anovaLPP$ANOVA$SSn[8] / (anovaLPP$ANOVA$SSd[8]+anovaLPP$ANOVA$SSn[8])
); print(anovaLPP)
capture.output(print(anovaLPP), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_acrossGroups_anovaFreq.doc"))

# bayesian ANOVA on LPP across conditioning groups
set.seed(rngSeed); anovaBFLPP <- generalTestBF(
  formula = LPP ~ usGroup*CS*time + partInd + partInd:CS + partInd:time,
  data = dataLPPLong,
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 10000 # only 10,000 iterations because it has to compute 128 models
); print(anovaBFLPP)
capture.output(print(anovaBFLPP), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_acrossGroups_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFLPP)
capture.output(bf_inclusion(anovaBFLPP), file = paste0(pathname, "/supplement/04s_lpp_timeFactor_acrossGroups_BFinclusion.doc"))


# quick graph of US Group x CS Type x Time ANOVA for LPP across groups
plotLPP <- ezPlot(
  data = dataLPPLong,
  dv = LPP,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,
  col = usGroup
) ; plotLPP
ggsave(plot = plotLPP, filename = paste0(pathname, "/supplement/04s_lpp_timeFactor_acrossGroups_plot.jpg"),
       width = 20, height = 10, units = "cm")
