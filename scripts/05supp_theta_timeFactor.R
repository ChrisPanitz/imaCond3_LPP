# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Mar 2024
# --- content: Supplementary analyses on frontomedial Theta, including factor "Time"

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
library(flextable) # ver. 0.9.4
library(stringr) #  ver. 1.5.0
library(here) # ver. 1.0.1


########################
### data preparation ###
########################

# load rating data from text file
# (see imaCond3_allratings_readme.txt for more details)
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
  #Av_allTr = importTheta$Fz_acqTotal_csplusav,
  Av_1stBl = importTheta$Fz_acq1_csplusav,
  Av_2ndBl = importTheta$Fz_acq2_csplusav,
  #Neu_allTr = importTheta$Fz_acqTotal_csplusneu,
  Neu_1stBl = importTheta$Fz_acq1_csplusneu,
  Neu_2ndBl = importTheta$Fz_acq2_csplusneu,
  #Min_allTr = importTheta$Fz_acqTotal_csminus,
  Min_1stBl = importTheta$Fz_acq1_csminus,
  Min_2ndBl = importTheta$Fz_acq2_csminus
)  
dataThetaLong <- gather(data = dataTheta, key = "cond", value = "theta",
                        Av_1stBl:Min_2ndBl)
dataThetaLong <- separate(data = dataThetaLong, col = cond,
                        into = c("CS","time"), sep = "_")
dataThetaLong$CS <- factor(dataThetaLong$CS)
dataThetaLong$time <- factor(dataThetaLong$time)



##########################################
### Imagery-based conditioning - theta ###
##########################################

# descriptive statistics for theta at Fz in imagery-based conditioning group
describe(dataTheta[dataTheta$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = theta at Fz
anovaThetaIma <- ezANOVA(
  data = dataThetaLong[dataThetaLong$usGroup == "ima",],
  dv = theta,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaThetaIma$ANOVA$pEtaSq <-
    c(anovaThetaIma$ANOVA$SSn[1] / (anovaThetaIma$ANOVA$SSd[1]+anovaThetaIma$ANOVA$SSn[1]),
      anovaThetaIma$ANOVA$SSn[2] / (anovaThetaIma$ANOVA$SSd[2]+anovaThetaIma$ANOVA$SSn[2]),
      anovaThetaIma$ANOVA$SSn[3] / (anovaThetaIma$ANOVA$SSd[3]+anovaThetaIma$ANOVA$SSn[3]),
      anovaThetaIma$ANOVA$SSn[4] / (anovaThetaIma$ANOVA$SSd[4]+anovaThetaIma$ANOVA$SSn[4])
); print(anovaThetaIma)
capture.output(print(anovaThetaIma), file = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_anovaFreq.doc"))

# bayesian CS x Time ANOVA on LPP in imagery-based conditioning group
set.seed(rngSeed); anovaBFThetaIma <- generalTestBF(
  formula = theta ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataThetaLong[dataThetaLong$usGroup == "ima",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
)

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFThetaIma@bayesFactor$bf <- log(exp(anovaBFThetaIma@bayesFactor$bf) / 
                                      exp(anovaBFThetaIma@bayesFactor$bf[length(anovaBFThetaIma@bayesFactor$bf)]))
anovaBFThetaIma@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFThetaIma)
capture.output(print(anovaBFThetaIma), file = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFThetaIma)
capture.output(bf_inclusion(anovaBFThetaIma), file = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_BFinclusion.doc"))

# quick graph of CS Type x time ANOVA for Theta at Fz in imagery-based conditioning group
plotThetaIma <- ezPlot(
  data = dataThetaLong[dataThetaLong$usGroup == "ima",],
  dv = theta,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotThetaIma 
ggsave(plot = plotThetaIma, filename = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_plot.jpg"),
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on Theta at Fz in imagery-based conditioning group
### 1stBl
# CS+av vs CS+neu
thetaImaAvNeu1stBl_t <- t.test(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                             y = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvNeu1stBl_d <- cohens_d(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvNeu1stBl_BF <- ttestBF(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaImaAvMin1stBl_t <- t.test(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                             y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvMin1stBl_d <- cohens_d(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvMin1stBl_BF <- ttestBF(x = dataTheta$Av_1stBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaImaNeuMin1stBl_t <- t.test(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                              y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
thetaImaNeuMin1stBl_d <- cohens_d(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                                paired = TRUE)
thetaImaNeuMin1stBl_BF <- ttestBF(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_1stBl[dataTheta$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
thetaImaAvNeu2ndBl_t <- t.test(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                             y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvNeu2ndBl_d <- cohens_d(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvNeu2ndBl_BF <- ttestBF(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaImaAvMin2ndBl_t <- t.test(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                             y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
thetaImaAvMin2ndBl_d <- cohens_d(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                               paired = TRUE)
thetaImaAvMin2ndBl_BF <- ttestBF(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "ima"],
                               y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaImaNeuMin2ndBl_t <- t.test(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                              y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
thetaImaNeuMin2ndBl_d <- cohens_d(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                                paired = TRUE)
thetaImaNeuMin2ndBl_BF <- ttestBF(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "ima"],
                                y = dataTheta$Min_2ndBl[dataTheta$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided

tableThetaIma <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(thetaImaAvNeu1stBl_t$statistic, thetaImaAvMin1stBl_t$statistic, thetaImaNeuMin1stBl_t$statistic,
        thetaImaAvNeu2ndBl_t$statistic, thetaImaAvMin2ndBl_t$statistic, thetaImaNeuMin2ndBl_t$statistic),
  df = c(thetaImaAvNeu1stBl_t$parameter, thetaImaAvMin1stBl_t$parameter, thetaImaNeuMin1stBl_t$parameter,
         thetaImaAvNeu2ndBl_t$parameter, thetaImaAvMin2ndBl_t$parameter, thetaImaNeuMin2ndBl_t$parameter), 
  p = c(thetaImaAvNeu1stBl_t$p.value, thetaImaAvMin1stBl_t$p.value, thetaImaNeuMin1stBl_t$p.value,
        thetaImaAvNeu2ndBl_t$p.value, thetaImaAvMin2ndBl_t$p.value, thetaImaNeuMin2ndBl_t$p.value),
  d = c(thetaImaAvNeu1stBl_d$Cohens_d, thetaImaAvMin1stBl_d$Cohens_d, thetaImaNeuMin1stBl_d$Cohens_d,
        thetaImaAvNeu2ndBl_d$Cohens_d, thetaImaAvMin2ndBl_d$Cohens_d, thetaImaNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(thetaImaAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(thetaImaAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(thetaImaNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(thetaImaAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(thetaImaAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(thetaImaNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableThetaIma, file = paste0(pathname, "/supplement/05s_theta_timeFactor_ima_tTable.doc"))



#######################################
### Classical conditioning - theta ###
#######################################

# descriptive statistics for theta at Fz in imagery-based conditioning group
describe(dataTheta[dataTheta$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = theta at Fz
anovaThetaReal <- ezANOVA(
  data = dataThetaLong[dataThetaLong$usGroup == "real",],
  dv = theta,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaThetaReal$ANOVA$pEtaSq <-
  c(anovaThetaReal$ANOVA$SSn[1] / (anovaThetaReal$ANOVA$SSd[1]+anovaThetaReal$ANOVA$SSn[1]),
    anovaThetaReal$ANOVA$SSn[2] / (anovaThetaReal$ANOVA$SSd[2]+anovaThetaReal$ANOVA$SSn[2]),
    anovaThetaReal$ANOVA$SSn[3] / (anovaThetaReal$ANOVA$SSd[3]+anovaThetaReal$ANOVA$SSn[3]),
    anovaThetaReal$ANOVA$SSn[4] / (anovaThetaReal$ANOVA$SSd[4]+anovaThetaReal$ANOVA$SSn[4])
  ); print(anovaThetaReal)
capture.output(print(anovaThetaReal), file = paste0(pathname, "/supplement/05s_theta_timeFactor_real_anovaFreq.doc"))

# bayesian ANOVA on theta at Fz in classical conditioning group
set.seed(rngSeed); anovaBFThetaReal <- generalTestBF(
  formula = theta ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataThetaLong[dataThetaLong$usGroup == "real",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
)

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFThetaReal@bayesFactor$bf <- log(exp(anovaBFThetaReal@bayesFactor$bf) / 
                                        exp(anovaBFThetaReal@bayesFactor$bf[length(anovaBFThetaReal@bayesFactor$bf)]))
anovaBFThetaReal@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFThetaReal)
capture.output(print(anovaBFThetaReal), file = paste0(pathname, "/supplement/05s_theta_timeFactor_real_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFThetaReal)
capture.output(bf_inclusion(anovaBFThetaReal), file = paste0(pathname, "/supplement/05s_theta_timeFactor_real_BFinclusion.doc"))

# quick graph of CS Type x time ANOVA for Theta at Fz in classical conditioning group
plotThetaReal <- ezPlot(
  data = dataThetaLong[dataThetaLong$usGroup == "real",],
  dv = theta,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotThetaReal 
ggsave(plot = plotThetaReal, filename = paste0(pathname, "/supplement/05s_theta_timeFactor_real_plot.jpg"),
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on Theta at Fz in classical conditioning group
### 1stBl
# CS+av vs CS+neu
thetaRealAvNeu1stBl_t <- t.test(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                               y = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvNeu1stBl_d <- cohens_d(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvNeu1stBl_BF <- ttestBF(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaRealAvMin1stBl_t <- t.test(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                               y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvMin1stBl_d <- cohens_d(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvMin1stBl_BF <- ttestBF(x = dataTheta$Av_1stBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaRealNeuMin1stBl_t <- t.test(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                                y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
thetaRealNeuMin1stBl_d <- cohens_d(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                                  paired = TRUE)
thetaRealNeuMin1stBl_BF <- ttestBF(x = dataTheta$Neu_1stBl[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_1stBl[dataTheta$usGroup == "real"],
                                  nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
thetaRealAvNeu2ndBl_t <- t.test(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                               y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvNeu2ndBl_d <- cohens_d(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvNeu2ndBl_BF <- ttestBF(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
thetaRealAvMin2ndBl_t <- t.test(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                               y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
thetaRealAvMin2ndBl_d <- cohens_d(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                                 paired = TRUE)
thetaRealAvMin2ndBl_BF <- ttestBF(x = dataTheta$Av_2ndBl[dataTheta$usGroup == "real"],
                                 y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
thetaRealNeuMin2ndBl_t <- t.test(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                                y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
thetaRealNeuMin2ndBl_d <- cohens_d(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                                  paired = TRUE)
thetaRealNeuMin2ndBl_BF <- ttestBF(x = dataTheta$Neu_2ndBl[dataTheta$usGroup == "real"],
                                  y = dataTheta$Min_2ndBl[dataTheta$usGroup == "real"],
                                  nullInterval = NULL, paired = TRUE) # two-sided

tableThetaReal <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(thetaRealAvNeu1stBl_t$statistic, thetaRealAvMin1stBl_t$statistic, thetaRealNeuMin1stBl_t$statistic,
        thetaRealAvNeu2ndBl_t$statistic, thetaRealAvMin2ndBl_t$statistic, thetaRealNeuMin2ndBl_t$statistic),
  df = c(thetaRealAvNeu1stBl_t$parameter, thetaRealAvMin1stBl_t$parameter, thetaRealNeuMin1stBl_t$parameter,
         thetaRealAvNeu2ndBl_t$parameter, thetaRealAvMin2ndBl_t$parameter, thetaRealNeuMin2ndBl_t$parameter), 
  p = c(thetaRealAvNeu1stBl_t$p.value, thetaRealAvMin1stBl_t$p.value, thetaRealNeuMin1stBl_t$p.value,
        thetaRealAvNeu2ndBl_t$p.value, thetaRealAvMin2ndBl_t$p.value, thetaRealNeuMin2ndBl_t$p.value),
  d = c(thetaRealAvNeu1stBl_d$Cohens_d, thetaRealAvMin1stBl_d$Cohens_d, thetaRealNeuMin1stBl_d$Cohens_d,
        thetaRealAvNeu2ndBl_d$Cohens_d, thetaRealAvMin2ndBl_d$Cohens_d, thetaRealNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(thetaRealAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(thetaRealAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(thetaRealNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(thetaRealAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(thetaRealAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(thetaRealNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableThetaReal, file = paste0(pathname, "/supplement/05s_theta_timeFactor_real_tTable.doc"))



###################################
### Across groups - theta at Fz ###
###################################

# descriptive statistics for Theta at Fz across conditioning groups
describe(dataTheta)

# frequentist CS x Time ANOVA on Theta at Fz across conditioning groups
anovaTheta <- ezANOVA(
  data = dataThetaLong,
  dv = theta,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaTheta$ANOVA$pEtaSq <- c(
  anovaTheta$ANOVA$SSn[1] / (anovaTheta$ANOVA$SSd[1]+anovaTheta$ANOVA$SSn[1]),
  anovaTheta$ANOVA$SSn[2] / (anovaTheta$ANOVA$SSd[2]+anovaTheta$ANOVA$SSn[2]),
  anovaTheta$ANOVA$SSn[3] / (anovaTheta$ANOVA$SSd[3]+anovaTheta$ANOVA$SSn[3]),
  anovaTheta$ANOVA$SSn[4] / (anovaTheta$ANOVA$SSd[4]+anovaTheta$ANOVA$SSn[4]),
  anovaTheta$ANOVA$SSn[5] / (anovaTheta$ANOVA$SSd[5]+anovaTheta$ANOVA$SSn[5]),
  anovaTheta$ANOVA$SSn[6] / (anovaTheta$ANOVA$SSd[6]+anovaTheta$ANOVA$SSn[6]),
  anovaTheta$ANOVA$SSn[7] / (anovaTheta$ANOVA$SSd[7]+anovaTheta$ANOVA$SSn[7]),
  anovaTheta$ANOVA$SSn[8] / (anovaTheta$ANOVA$SSd[8]+anovaTheta$ANOVA$SSn[8])
); print(anovaTheta)
capture.output(print(anovaTheta), file = paste0(pathname, "/supplement/05s_theta_timeFactor_acrossGroups_anovaFreq.doc"))

# bayesian ANOVA on Theta at Fz across conditioning groups
set.seed(rngSeed); anovaBFTheta <- generalTestBF(
  formula = theta ~ usGroup*CS*time + partInd + partInd:CS + partInd:time,
  data = dataThetaLong,
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 10000 # only 10,000 iterations because it has to compute 128 models
)

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFTheta@bayesFactor$bf <- log(exp(anovaBFTheta@bayesFactor$bf) / 
                                        exp(anovaBFTheta@bayesFactor$bf[length(anovaBFTheta@bayesFactor$bf)]))
anovaBFTheta@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFTheta)
capture.output(print(anovaBFTheta), file = paste0(pathname, "/supplement/05s_theta_timeFactor_acrossGroups_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFTheta)
capture.output(bf_inclusion(anovaBFTheta), file = paste0(pathname, "/supplement/05s_theta_timeFactor_acrossGroups_BFinclusion.doc"))

# quick graph of group x CS ANOVA on Theta at Fz
plotTheta <- ezPlot(
  data = dataThetaLong,
  dv = theta,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,
  col = usGroup
); plotTheta
ggsave(plot = plotTheta, filename = paste0(pathname, "/supplement/05s_theta_timeFactor_acrossGroups_plot.jpg"),
       width = 20, height = 10, units = "cm")
