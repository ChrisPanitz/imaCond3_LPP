# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Mar 2024
# --- content: Supplementary analyses on frontomedial Theta, parallel to main analyses (ANOVAs, pairwise comparisons, plotting)

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
library(bayestestR) # 0.13.1
library(stringr) # ver. 1.5.0
library(here) # ver. 1.0.1
library(ggplot2) # ver. 3.4.2


########################
### data preparation ###
########################

# load rating data from text file
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# load Theta data from text file
# (see imaCond3_theta_readme.txt for more details)
importTheta <- read.csv(paste0(pathname, "/experimentData/imaCond3_theta.txt"), sep = ",")
importTheta <- merge(importRatings[,c("partCode","group")], importTheta, by = "partCode")
importTheta <- importTheta[order(importTheta$partCode),]

# create data frames in wide & long format for fear ratings
dataTheta <- data.frame(
  partInd = factor(1:dim(importTheta)[1]),
  usGroup = factor(importTheta$group, labels = c("ima", "real")),
  Av_allTr = importTheta$Fz_acqTotal_csplusav,
  Av_1stBl = importTheta$Fz_acq1_csplusav,
  Av_2ndBl = importTheta$Fz_acq2_csplusav,
  Neu_allTr = importTheta$Fz_acqTotal_csplusneu,
  Neu_1stBl = importTheta$Fz_acq1_csplusneu,
  Neu_2ndBl = importTheta$Fz_acq2_csplusneu,
  Min_allTr = importTheta$Fz_acqTotal_csminus,
  Min_1stBl = importTheta$Fz_acq1_csminus,
  Min_2ndBl = importTheta$Fz_acq2_csminus
)  
dataThetaLong <- gather(data = dataTheta, key = "cond", value = "theta",
                        Av_allTr:Min_2ndBl)
dataThetaLong <- separate(data = dataThetaLong, col = cond,
                        into = c("CS", "time"), sep = "_")
dataThetaLong$CS <- factor(dataThetaLong$CS)



##########################################
### Imagery-based conditioning - theta ###
##########################################

# descriptive statistics for theta at Fz in imagery-based conditioning group
describe(dataTheta[dataTheta$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = theta at Fz
anovaThetaIma <- ezANOVA(
  data = dataThetaLong[dataThetaLong$usGroup == "ima" & dataThetaLong$time == "allTr",],
  dv = theta,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaThetaIma$ANOVA$pEtaSq <- 
  c(anovaThetaIma$ANOVA$SSn[1] / (anovaThetaIma$ANOVA$SSd[1]+anovaThetaIma$ANOVA$SSn[1]),
    anovaThetaIma$ANOVA$SSn[2] / (anovaThetaIma$ANOVA$SSd[2]+anovaThetaIma$ANOVA$SSn[2])
  ); print(anovaThetaIma)
capture.output(print(anovaThetaIma), file = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_anovaFreq.doc"))

# bayesian CS x Time ANOVA on theta in imagery-based conditioning group
set.seed(rngSeed); anovaBFThetaIma <- anovaBF(
  formula = theta ~ CS + partInd,
  data = dataThetaLong[dataThetaLong$usGroup == "ima" & dataThetaLong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFThetaIma)
capture.output(print(anovaBFThetaIma), file = paste0(pathname, "/supplement/05s_theta_ima_anovaBayes.doc"))

# frequentist & bayesian t-tests on Theta at Fz in imagery-based conditioning group
# CS+av vs CS+neu
thetaImaAvNeuallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                             y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvNeuallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvNeuallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaImaAvMinallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                             y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvMinallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvMinallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaImaNeuMinallTr_t <- t.test(x = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                              y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
thetaImaNeuMinallTr_d <- cohens_d(x = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                paired = TRUE)
thetaImaNeuMinallTr_BF <- ttestBF(x = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided


tableThetaIma <- data.frame(
  comparison = c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"),
  t = c(thetaImaAvNeuallTr_t$statistic, thetaImaAvMinallTr_t$statistic, thetaImaNeuMinallTr_t$statistic),
  df = c(thetaImaAvNeuallTr_t$parameter, thetaImaAvMinallTr_t$parameter, thetaImaNeuMinallTr_t$parameter), 
  p = c(thetaImaAvNeuallTr_t$p.value, thetaImaAvMinallTr_t$p.value, thetaImaNeuMinallTr_t$p.value),
  d = c(thetaImaAvNeuallTr_d$Cohens_d, thetaImaAvMinallTr_d$Cohens_d, thetaImaNeuMinallTr_d$Cohens_d),
  BF = c(exp(thetaImaAvNeuallTr_BF@bayesFactor[["bf"]][1]), exp(thetaImaAvMinallTr_BF@bayesFactor[["bf"]][1]), exp(thetaImaNeuMinallTr_BF@bayesFactor[["bf"]][1])),
  testDir = c("one.sided","one.sided","two.sided")
)
capture.output(tableThetaIma, file = paste0(pathname, "/supplement/05s_theta_ima_tTable.doc"))



#######################################
### Classical conditioning - theta ###
#######################################

# descriptive statistics for theta at Fz in imagery-based conditioning group
describe(dataTheta[dataTheta$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = theta at Fz
anovaThetaReal <- ezANOVA(
  data = dataThetaLong[dataThetaLong$usGroup == "real" & dataThetaLong$time == "allTr",],
  dv = theta,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaThetaReal$ANOVA$pEtaSq <- 
  c(anovaThetaReal$ANOVA$SSn[1] / (anovaThetaReal$ANOVA$SSd[1]+anovaThetaReal$ANOVA$SSn[1]),
    anovaThetaReal$ANOVA$SSn[2] / (anovaThetaReal$ANOVA$SSd[2]+anovaThetaReal$ANOVA$SSn[2])
  ); print(anovaThetaReal)
capture.output(print(anovaThetaReal), file = paste0(pathname, "/supplement/05s_theta_timeFactor_real_anovaFreq.doc"))

# bayesian CS x Time ANOVA on theta in classical conditioning group
set.seed(rngSeed); anovaBFThetaReal <- anovaBF(
  formula = theta ~ CS + partInd,
  data = dataThetaLong[dataThetaLong$usGroup == "real" & dataThetaLong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFThetaReal)
capture.output(print(anovaBFThetaReal), file = paste0(pathname, "/supplement/05s_theta_real_anovaBayes.doc"))

# frequentist & bayesian t-tests on Theta at Fz in classical conditioning group
# CS+av vs CS+neu
thetaRealAvNeuallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                               y = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvNeuallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvNeuallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaRealAvMinallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                               y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvMinallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvMinallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaRealNeuMinallTr_t <- t.test(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
thetaRealNeuMinallTr_d <- cohens_d(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                  paired = TRUE)
thetaRealNeuMinallTr_BF <- ttestBF(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                  nullInterval = NULL, paired = TRUE) # two-sided


tableThetaReal <- data.frame(
  comparison = c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"),
  t = c(thetaRealAvNeuallTr_t$statistic, thetaRealAvMinallTr_t$statistic, thetaRealNeuMinallTr_t$statistic),
  df = c(thetaRealAvNeuallTr_t$parameter, thetaRealAvMinallTr_t$parameter, thetaRealNeuMinallTr_t$parameter), 
  p = c(thetaRealAvNeuallTr_t$p.value, thetaRealAvMinallTr_t$p.value, thetaRealNeuMinallTr_t$p.value),
  d = c(thetaRealAvNeuallTr_d$Cohens_d, thetaRealAvMinallTr_d$Cohens_d, thetaRealNeuMinallTr_d$Cohens_d),
  BF = c(exp(thetaRealAvNeuallTr_BF@bayesFactor[["bf"]][1]), exp(thetaRealAvMinallTr_BF@bayesFactor[["bf"]][1]), exp(thetaRealNeuMinallTr_BF@bayesFactor[["bf"]][1])),
  testDir = c("one.sided","one.sided","two.sided")
)
capture.output(tableThetaReal, file = paste0(pathname, "/supplement/05s_theta_real_tTable.doc"))



###################################
### Across groups - theta at Fz ###
###################################

# descriptive statistics for Theta at Fz across conditioning groups
describe(dataTheta)

# frequentist CS x Time ANOVA on Theta at Fz across conditioning groups
anovaTheta <- ezANOVA(
  data = dataThetaLong[dataThetaLong$time == "allTr",],
  dv = theta,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaTheta$ANOVA$pEtaSq <- c(
  anovaTheta$ANOVA$SSn[1] / (anovaTheta$ANOVA$SSd[1]+anovaTheta$ANOVA$SSn[1]),
  anovaTheta$ANOVA$SSn[2] / (anovaTheta$ANOVA$SSd[2]+anovaTheta$ANOVA$SSn[2]),
  anovaTheta$ANOVA$SSn[3] / (anovaTheta$ANOVA$SSd[3]+anovaTheta$ANOVA$SSn[3]),
  anovaTheta$ANOVA$SSn[4] / (anovaTheta$ANOVA$SSd[4]+anovaTheta$ANOVA$SSn[4])
); print(anovaTheta)
capture.output(print(anovaTheta), file = paste0(pathname, "/supplement/05s_theta_acrossGroups_anovaFreq.doc"))

# bayesian ANOVA on Theta at Fz across conditioning groups
set.seed(rngSeed); anovaBFTheta <- anovaBF(
  formula = theta ~ usGroup*CS + partInd,
  data = dataThetaLong[dataThetaLong$time == "allTr",],
  whichRandom = c("partInd"),
  whichModels = "all",
  iterations = 100000
); print(anovaBFTheta)
capture.output(print(anovaBFTheta), file = paste0(pathname, "/supplement/05s_theta_acrossGroups_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFTheta)
capture.output(bf_inclusion(anovaBFTheta), file = paste0(pathname, "/supplement/05s_theta_acrossGroups_BFinclusion.doc"))

# quick graph of group x CS ANOVA on Theta at Fz
plotTheta <- ezPlot(
  data = dataThetaLong[dataThetaLong$time == "allTr",],
  dv = theta,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
); plotTheta
ggsave(plot = plotTheta, filename = paste0(pathname, "/supplement/05s_theta_acrossGroups_plot.jpg"),
       width = 20, height = 10, units = "cm")

# frequentist & bayesian t-tests on Theta (difference scores) across groups
# delta [CS+av - CS+neu]
thetaBothAvNeuallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                              y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
thetaBothAvNeuallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                  dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                  dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                                paired = FALSE)
thetaBothAvNeuallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                  dataTheta$Neu_allTr[dataTheta$usGroup == "real"],
                                y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                  dataTheta$Neu_allTr[dataTheta$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
thetaBothAvMinallTr_t <- t.test(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                              y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
thetaBothAvMinallTr_d <- cohens_d(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                  dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                  dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                paired = FALSE)
thetaBothAvMinallTr_BF <- ttestBF(x = dataTheta$Av_allTr[dataTheta$usGroup == "real"] -
                                  dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                y = dataTheta$Av_allTr[dataTheta$usGroup == "ima"] -
                                  dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
thetaBothNeuMinallTr_t <- t.test(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"] -
                                 dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                               y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"] - 
                                 dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                               alternative = "two.sided", paired = FALSE) # two-sided
thetaBothNeuMinallTr_d <- cohens_d(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"] -
                                   dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"] -
                                   dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                 paired = FALSE)
thetaBothNeuMinallTr_BF <- ttestBF(x = dataTheta$Neu_allTr[dataTheta$usGroup == "real"] -
                                   dataTheta$Min_allTr[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_allTr[dataTheta$usGroup == "ima"] - 
                                   dataTheta$Min_allTr[dataTheta$usGroup == "ima"],
                                 nullInterval = NULL, paired = FALSE) # two-sided

tableThetaBoth <- data.frame(
  comparison = c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"),
  t = c(thetaBothAvNeuallTr_t$statistic, thetaBothAvMinallTr_t$statistic, thetaBothNeuMinallTr_t$statistic),
  df = c(thetaBothAvNeuallTr_t$parameter, thetaBothAvMinallTr_t$parameter, thetaBothNeuMinallTr_t$parameter), 
  p = c(thetaBothAvNeuallTr_t$p.value*3, thetaBothAvMinallTr_t$p.value*3, thetaBothNeuMinallTr_t$p.value*3), # Bonferroni
  d = c(thetaBothAvNeuallTr_d$Cohens_d, thetaBothAvMinallTr_d$Cohens_d, thetaBothNeuMinallTr_d$Cohens_d),
  BF = c(exp(thetaBothAvNeuallTr_BF@bayesFactor[["bf"]][1]), exp(thetaBothAvMinallTr_BF@bayesFactor[["bf"]][1]), exp(thetaBothNeuMinallTr_BF@bayesFactor[["bf"]][1])),
  testDir = rep("two.sided",3)
)
tableThetaBoth$p[tableThetaBoth$p > 1] <- 1
capture.output(tableThetaBoth, file = paste0(pathname, "/supplement/05s_theta_acrossGroups_tTable.doc"))