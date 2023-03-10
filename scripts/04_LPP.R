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
if(!is.element("eegUtils",installed.packages()[,1])) {install.packages("eegUtils")}
library(eegUtils) # 
if(!is.element("OneR",installed.packages()[,1])) {install.packages("OneR")}
library(OneR) # 


###################################
### Setting time window for LPP ###
###################################

# set by user
sRate <- 1024 # sampling rate of EEG data
startSeg <- -200 # time of first data point relative to CS
TWOI <- c(300,700) # Time Window Of Interest (in ms)
chanInd <- 31 # 31 = index of Pz in EEG array

# computed
# Sample Window Of Interest
SWOI <- c(round(TWOI[1]-startSeg*sRate/1000),round(TWOI[2]-startSeg*sRate/1000)) 



########################
### data preparation ###
########################

# Rating data contains group membership
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# read LPP data
# create emtpy matrices, rows = participants, columns = sample points (entire ERP)
tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[1], "_IMAKON03_Average_Akquges_EKP_P9P10_S51.dat"), 
                     header = FALSE, sep = " ")
lppMatAV <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatNEU <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])
lppMatMIN <- matrix(data = NA, nrow = dim(importRatings)[1], ncol = dim(tempData)[1])

# read single-subject data and extract Pz
for (partI in 1:length(importRatings$partCode)) {
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S51.dat"), 
                       header = FALSE, sep = " ")
  lppMatAV[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S52.dat"), 
                       header = FALSE, sep = " ")
  lppMatNEU[partI,] <- tempData[,chanInd]
  tempData <- read.csv(paste0(pathname, "/experimentData/erpData/", importRatings$partCode[partI], "_IMAKON03_Average_Akquges_EKP_P9P10_S53.dat"), 
                       header = FALSE, sep = " ")
  lppMatMIN[partI,] <- tempData[,chanInd]
}

# create LPP data frame; LPP values are means from 300 to 700 s post-CS
dataLPP <- data.frame(
  partCode = factor(importRatings$partCode),
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_allTr = rowMeans(lppMatAV[,SWOI[1]:SWOI[2]]),
  Neu_allTr = rowMeans(lppMatNEU[,SWOI[1]:SWOI[2]]),
  Min_allTr = rowMeans(lppMatMIN[,SWOI[1]:SWOI[2]])
)  

# transform into long format
dataLPPLong <- gather(data = dataLPP, key = "cond", value = "LPP", Av_allTr:Min_allTr)
dataLPPLong <- separate(data = dataLPPLong, col = cond, into = c("CS","time"), sep = "_")
dataLPPLong$CS <- factor(dataLPPLong$CS, levels = c("Av","Neu","Min"))
dataLPPLong$time <- factor(dataLPPLong$time)

# create dataframes for topographies


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
                         alternative = "greater", paired = TRUE) # one-sided
lppImaNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           paired = TRUE)
lppImaNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           y = dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y



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
                          alternative = "two.sided", paired = TRUE) # one-sided
lppRealNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            paired = TRUE)
lppRealNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # one-sided x > y



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
bf_nullModel <- 1
bf_usGroup <- exp(anovaBFLPP@bayesFactor$bf[1])
bf_cs <- exp(anovaBFLPP@bayesFactor$bf[2])
bf_interact <- exp(anovaBFLPP@bayesFactor$bf[3])
bf_usGroup_cs <- exp(anovaBFLPP@bayesFactor$bf[4])
bf_usGroup_interact <- exp(anovaBFLPP@bayesFactor$bf[5])
bf_cs_interact <- exp(anovaBFLPP@bayesFactor$bf[6])
bf_fullModel <- exp(anovaBFLPP@bayesFactor$bf[7])

# main effect US group: models [1] and [3] vs. null model and model [2]
bfIncGroupLPP <- (bf_usGroup + bf_usGroup_cs + bf_usGroup_interact + bf_fullModel) / 
  (bf_nullModel + bf_cs + bf_interact + bf_cs_interact); bfIncGroupLPP
# main effect CS type: models "main effect CS" & "main effects CS & group" vs.
#                      null model and "main effect group"
bfIncCsLPP <- (bf_cs + bf_usGroup_cs + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_interact + bf_usGroup_interact); bfIncCsLPP
# interaction: Full model vs. main-effects-only model
bfIncInteractLPP <- (bf_interact + bf_usGroup_interact + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_cs + bf_usGroup_cs); bfIncInteractLPP

# quick & dirty graph of group x CS ANOVA on LPP
ezPlot(
  data = dataLPPLong[dataLPPLong$time == "allTr",],
  dv = LPP,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on LPP (difference scores) across groups
# delta [CS+av - CS+neu]
lppBothAvNeu_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                             dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                         y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                             dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
lppBothAvNeu_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                               dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                               dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           paired = FALSE)
lppBothAvNeu_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                               dataLPP$Neu_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                               dataLPP$Neu_allTr[dataLPP$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
lppBothAvMin_t <- t.test(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                             dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                         y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                             dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
lppBothAvMin_d <- cohens_d(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           paired = FALSE)
lppBothAvMin_BF <- ttestBF(x = dataLPP$Av_allTr[dataLPP$usGroup == "real"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                           y = dataLPP$Av_allTr[dataLPP$usGroup == "ima"] -
                               dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
lppBothNeuMin_t <- t.test(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                              dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                          y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                              dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
lppBothNeuMin_d <- cohens_d(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                            paired = FALSE)
lppBothNeuMin_BF <- ttestBF(x = dataLPP$Neu_allTr[dataLPP$usGroup == "real"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "real"],
                            y = dataLPP$Neu_allTr[dataLPP$usGroup == "ima"] -
                                dataLPP$Min_allTr[dataLPP$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided



# quick & dirty graph of group x CS ANOVA on valence ratings
ezPlot(
  data = dataLPPLong,
  dv = LPP,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  




#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-",
                     "groups: delta CS+av / CS+neu", "groups: delta CS+av / CS-", "groups: delta CSneu / CS-"), 3),
  t = c(lppImaAvNeu_t$statistic, lppImaAvMin_t$statistic, lppImaNeuMin_t$statistic,
        lppRealAvNeu_t$statistic, lppRealAvMin_t$statistic, lppRealNeuMin_t$statistic,
        lppBothAvNeu_t$statistic, lppBothAvMin_t$statistic, lppBothNeuMin_t$statistic),
  df = c(lppImaAvNeu_t$parameter, lppImaAvMin_t$parameter, lppImaNeuMin_t$parameter,
         lppRealAvNeu_t$parameter, lppRealAvMin_t$parameter, lppRealNeuMin_t$parameter,
         lppBothAvNeu_t$parameter, lppBothAvMin_t$parameter, lppBothNeuMin_t$parameter), 
  p = c(lppImaAvNeu_t$p.value, lppImaAvMin_t$p.value, lppImaNeuMin_t$p.value,
        lppRealAvNeu_t$p.value, lppRealAvMin_t$p.value, lppRealNeuMin_t$p.value,
        lppBothAvNeu_t$p.value*3, lppBothAvMin_t$p.value*3, lppBothNeuMin_t$p.value*3),  # Bonferroni
  d = c(lppImaAvNeu_d$Cohens_d, lppImaAvMin_d$Cohens_d, lppImaNeuMin_d$Cohens_d,
        lppRealAvNeu_d$Cohens_d, lppRealAvMin_d$Cohens_d, lppRealNeuMin_d$Cohens_d,
        lppBothAvNeu_d$Cohens_d, lppBothAvMin_d$Cohens_d, lppBothNeuMin_d$Cohens_d),
  BF = c(exp(lppImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppImaAvMin_BF@bayesFactor[["bf"]][1]), exp(lppImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(lppRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppRealAvMin_BF@bayesFactor[["bf"]][1]), exp(lppRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(lppBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(lppBothAvMin_BF@bayesFactor[["bf"]][1]), exp(lppBothNeuMin_BF@bayesFactor[["bf"]][1]))
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
tableData$BF <- format(tableData$BF, digits = 2)

tableLPP <- flextable(tableData[1:9,])
tableLPP <- add_header_lines(tableLPP, top = TRUE, values = "LPP")
tableLPP <- align(tableLPP, align = "center")

save_as_docx(tableLPP, path = paste0(pathname, "/Tables/tableLPP_raw.docx"))



####################
### Plotting LPP ###
####################

csLabels = c("CS+av", "CS+neu", "CS-")
lppAv4GA <- lppMatAV
  lppAv4GAima <- lppMatAV[importRatings$group == "ima",]
  lppAv4GAreal <- lppMatAV[importRatings$group == "real",]
lppNeu4GA <- lppMatNEU
  lppNeu4GAima <- lppMatNEU[importRatings$group == "ima",]
  lppNeu4GAreal <- lppMatNEU[importRatings$group == "real",]
lppMin4GA <- lppMatMIN
  lppMin4GAima <- lppMatMIN[importRatings$group == "ima",]
  lppMin4GAreal <- lppMatMIN[importRatings$group == "real",]
  
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


dataLPPWithin <- dataLPP[,c("partInd","usGroup","Av_allTr","Neu_allTr","Min_allTr")]
# remove each participant's average from each single value
dataLPPWithin[,3:5] <- as.matrix(dataLPPWithin[,3:5]) -
  rowMeans(as.matrix(dataLPPWithin[,3:5])) 
# prepare data frame for bar plot with means from standard dataset and SE from
# dataset without betweem-subject variance
meanLPP <- data.frame(
  usGroup = factor(c(rep(1,3),rep(2,3)),
                   labels = c("Imagery-Based","Classical")),
  CS = factor(c(1,2,3,1,2,3),
              labels = c("CS+\nav","CS+\nneu","CS-\n")),
  #mean = c(describe(dataLPP[dataLPP$usGroup == "ima", c(3,6,9)])$mean,
  #         describe(dataLPP[dataLPP$usGroup == "real", c(3,6,9)])$mean),
  mean = c(describe(dataLPP[dataLPP$usGroup == "ima", c(4,5,6)])$mean,
           describe(dataLPP[dataLPP$usGroup == "real", c(4,5,6)])$mean),
  se = c(describe(dataLPPWithin[dataLPPWithin$usGroup == "ima", 3:5])$se,
         describe(dataLPPWithin[dataLPPWithin$usGroup == "real", 3:5])$se)
)




lineSize = 1
yMin = -5
yMax = 12
plotFS <- 9
showSig <- TRUE

graphLPPima <- ggplot(data = lppGAima, aes(x = time, y = LPP, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = TWOI[1], xmax = TWOI[2], ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), size = lineSize) + 
  scale_x_continuous(breaks = seq(-200,1000,200)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Imagery-Based Conditioning", x = "Time (ms)", y = "ERP amplitude at Pz (??V)", fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = FALSE) +
  theme(
    legend.position = "none",
    plot.title = element_blank(), #element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
    axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
    axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

graphLPPreal <- ggplot(data = lppGAreal, aes(x = time, y = LPP, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = TWOI[1], xmax = TWOI[2], ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), size = lineSize) + 
  scale_x_continuous(breaks = seq(-200,1000,200)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Classical Conditioning", x = "Time (ms)", y = "ERP amplitude at Pz (??V)", fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = FALSE) +
  theme(
    legend.position = "none",
    plot.title = element_blank(), #element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
    axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
    axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
    axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

graphLPPmeansIma <- ggplot(data = meanLPP[meanLPP$usGroup == "Imagery-Based",], aes(x = CS, y = mean, fill = CS)) +
  theme_classic() +
  geom_col(aes(fill = CS), position = position_dodge(width = .9)) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = .1), position = position_dodge(width = .9)) +
  scale_x_discrete(name = " ", labels = csLabels, position = "bottom") +
  scale_y_continuous(name = "Mean LPP amplitude (300-700 ms)", limits = c(0,14), breaks = c(0,5,10), expand = c(0,0)) +
  coord_cartesian(ylim = c(0,12), clip = 'off') +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black")); 

graphLPPmeansReal <- ggplot(data = meanLPP[meanLPP$usGroup == "Classical",], aes(x = CS, y = mean, fill = CS)) +
  theme_classic() +
  geom_col(aes(fill = CS), position = position_dodge(width = .9)) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = .1), position = position_dodge(width = .9)) +
  scale_x_discrete(name = " ", labels = csLabels, position = "bottom") +
  scale_y_continuous(name = "Mean LPP amplitude (300-700 ms)", limits = c(0,14), breaks = c(0,5,10), expand = c(0,0)) +
  coord_cartesian(ylim = c(0,12), clip = 'off') +
  theme(legend.position = "none",
        plot.title = element_blank(), #element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black")); 

if (showSig == TRUE){
  graphLPPmeansIma <- graphLPPmeansIma +
    geom_segment(aes(x = 1, y = mean+se+.5, xend = 2, yend = mean+se+.5), data = meanLPP[1,]) +
    geom_text(aes(label = "???", x = 1.5, y = mean+se+1), size = plotFS/4, data = meanLPP[1,]) +
    geom_segment(aes(x = 1, y = mean+se+1.5, xend = 3, yend = mean+se+1.5), data = meanLPP[1,]) +
    geom_text(aes(label = "*", x = 2, y = mean+se+1.6), size = plotFS/2, data = meanLPP[1,])
  graphLPPmeansReal <- graphLPPmeansReal +  
    geom_segment(aes(x = 1, y = mean+se+.5, xend = 2, yend = mean+se+.5), data = meanLPP[4,]) +
    geom_text(aes(label = "*", x = 1.5, y = mean+se+.6), size = plotFS/2, data = meanLPP[4,])
}

# combining graphs into one figure
graphLPPima <- graphLPPima + theme(plot.margin = unit(c(10,5,5,5), "mm"))
graphLPPreal <- graphLPPreal + theme(plot.margin = unit(c(10,5,5,5), "mm"))
graphLPPmeansIma <- graphLPPmeansIma + theme(plot.margin = unit(c(10,5,5,5), "mm"))
graphLPPmeansReal <- graphLPPmeansReal + theme(plot.margin = unit(c(10,5,5,5), "mm"))

graphLPProw1 <- ggarrange(graphLPPima, graphLPPmeansIma,
                          ncol = 2, nrow = 1, 
                          widths = c(3,2))
graphLPProw2 <- ggarrange(graphLPPreal, graphLPPmeansReal,
                          ncol = 2, nrow = 1, 
                          widths = c(3,2))
graphLPP <- ggarrange(graphLPProw1,graphLPProw2,
                      ncol = 1, nrow = 2,
                      labels = c("Imagery-Based Conditioning", "Classical Conditioning"),
                      label.x = c(0.09,0.18)
                      )
graphLPP

# saving it
ggsave(filename = paste0(pathname, "/Figures/Figure5_timeCourses_barPlot_LPP.eps"),
       plot = graphLPP,
       width = 150,
       height = 150,
       units = "mm",
       dpi = 300
)

ggsave(filename = paste0(pathname, "/Figures/Figure5_timeCourses_barPlot_LPP.pdf"),
       plot = graphLPP,
       width = 150,
       height = 150,
       units = "mm",
       dpi = 300
)

