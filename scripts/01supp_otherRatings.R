# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Mar 2024
# --- content: Main analyses on unpleasantness, arousal, anger, & disgust ratings (ANOVAs, pairwise comparisons, plotting)

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



########################
### data preparation ###
########################

# load rating data from text file
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")



# create data frames in wide & long format for unpleasantness ratings
dataUnpleas <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resolve issue in which Bayes ANOVA crashes using alphanumeric codes as participant ID
  partCode = factor(importRatings$partCode),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_Pre = importRatings$unpleas_csplus_av_2,
  Av_Mid = importRatings$unpleas_csplus_av_3,
  Av_Post = importRatings$unpleas_csplus_av_4,
  Neu_Pre = importRatings$unpleas_csplus_neu_2,
  Neu_Mid = importRatings$unpleas_csplus_neu_3,
  Neu_Post = importRatings$unpleas_csplus_neu_4,
  Min_Pre = importRatings$unpleas_csminus_2,
  Min_Mid = importRatings$unpleas_csminus_3,
  Min_Post = importRatings$unpleas_csminus_4
)  
dataUnpleasLong <- gather(data = dataUnpleas, key = "cond", value = "unpleasantness",
                          Av_Pre:Min_Post)
dataUnpleasLong <- separate(data = dataUnpleasLong, col = cond,
                            into = c("CS","time"), sep = "_")
dataUnpleasLong$CS <- factor(dataUnpleasLong$CS, levels = c("Av", "Neu", "Min"))
dataUnpleasLong$time <- factor(dataUnpleasLong$time, levels = c("Pre", "Mid", "Post"))

# create data frames in wide & long format for arousal ratings
dataArousal <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resolve issue in which Bayes ANOVA crashes using alphanumeric codes as participant ID
  partCode = factor(importRatings$partCode),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_Pre = importRatings$arous_csplus_av_2,
  Av_Mid = importRatings$arous_csplus_av_3,
  Av_Post = importRatings$arous_csplus_av_4,
  Neu_Pre = importRatings$arous_csplus_neu_2,
  Neu_Mid = importRatings$arous_csplus_neu_3,
  Neu_Post = importRatings$arous_csplus_neu_4,
  Min_Pre = importRatings$arous_csminus_2,
  Min_Mid = importRatings$arous_csminus_3,
  Min_Post = importRatings$arous_csminus_4
)  
dataArousalLong <- gather(data = dataArousal, key = "cond", value = "arousal",
                          Av_Pre:Min_Post)
dataArousalLong <- separate(data = dataArousalLong, col = cond,
                            into = c("CS","time"), sep = "_")
dataArousalLong$CS <- factor(dataArousalLong$CS, levels = c("Av", "Neu", "Min"))
dataArousalLong$time <- factor(dataArousalLong$time, levels = c("Pre", "Mid", "Post"))

# create data frames in wide & long format for anger ratings
dataAnger <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resolve issue in which Bayes ANOVA crashes using alphanumeric codes as participant ID
  partCode = factor(importRatings$partCode),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_Pre = importRatings$ang_csplus_av_2,
  Av_Mid = importRatings$ang_csplus_av_3,
  Av_Post = importRatings$ang_csplus_av_4,
  Neu_Pre = importRatings$ang_csplus_neu_2,
  Neu_Mid = importRatings$ang_csplus_neu_3,
  Neu_Post = importRatings$ang_csplus_neu_4,
  Min_Pre = importRatings$ang_csminus_2,
  Min_Mid = importRatings$ang_csminus_3,
  Min_Post = importRatings$ang_csminus_4
)  
dataAngerLong <- gather(data = dataAnger, key = "cond", value = "anger",
                       Av_Pre:Min_Post)
dataAngerLong <- separate(data = dataAngerLong, col = cond, into = c("CS","time"),
                         sep = "_")
dataAngerLong$CS <- factor(dataAngerLong$CS, levels = c("Av", "Neu", "Min"))
dataAngerLong$time <- factor(dataAngerLong$time, levels = c("Pre", "Mid", "Post"))

# create data frames in wide & long format for disgust ratings
dataDisgust <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resolve issue in which Bayes ANOVA crashes using alphanumeric codes as participant ID
  partCode = factor(importRatings$partCode),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
  Av_Pre = importRatings$dis_csplus_av_2,
  Av_Mid = importRatings$dis_csplus_av_3,
  Av_Post = importRatings$dis_csplus_av_4,
  Neu_Pre = importRatings$dis_csplus_neu_2,
  Neu_Mid = importRatings$dis_csplus_neu_3,
  Neu_Post = importRatings$dis_csplus_neu_4,
  Min_Pre = importRatings$dis_csminus_2,
  Min_Mid = importRatings$dis_csminus_3,
  Min_Post = importRatings$dis_csminus_4
)  
dataDisgustLong <- gather(data = dataDisgust, key = "cond", value = "disgust",
                        Av_Pre:Min_Post)
dataDisgustLong <- separate(data = dataDisgustLong, col = cond, into = c("CS","time"),
                          sep = "_")
dataDisgustLong$CS <- factor(dataDisgustLong$CS, levels = c("Av", "Neu", "Min"))
dataDisgustLong$time <- factor(dataDisgustLong$time, levels = c("Pre", "Mid", "Post"))



##############################################################################
### Imagery-based conditioning - unpleasantness ratings - primary analyses ###
##############################################################################

# descriptive statistics for unpleasantness ratings in imagery-based conditioning group
describe(dataUnpleas[dataUnpleas$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = unpleasantness rating
anovaUnpleasIma <- ezANOVA(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "ima" & dataUnpleasLong$time == "Post",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaUnpleasIma$ANOVA$pEtaSq <-
    c(anovaUnpleasIma$ANOVA$SSn[1] / (anovaUnpleasIma$ANOVA$SSd[1]+anovaUnpleasIma$ANOVA$SSn[1]),
      anovaUnpleasIma$ANOVA$SSn[2] / (anovaUnpleasIma$ANOVA$SSd[2]+anovaUnpleasIma$ANOVA$SSn[2])
); print(anovaUnpleasIma)

# bayesian ANOVA on unpleasantness ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFUnpleasIma <- anovaBF(
  formula = unpleasantness ~ CS + partInd,
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "ima" & dataUnpleasLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFUnpleasIma)

# frequentist & bayesian t-tests on unpleasantness ratings in imagery-based conditioning group
# CS+av vs CS+neu
unpleasImaAvNeu_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                       y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                       alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvNeu_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                              paired = TRUE)
unpleasImaAvNeu_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasImaAvMin_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                            y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvMin_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                              paired = TRUE)
unpleasImaAvMin_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasImaNeuMin_t <- t.test(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                             y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                             alternative = "two.sided", paired = TRUE) # two-sided
unpleasImaNeuMin_d <- cohens_d(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               paired = TRUE)
unpleasImaNeuMin_BF <- ttestBF(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided



#######################################################################
### Imagery-based conditioning - arousal ratings - primary analyses ###
#######################################################################

# descriptive statistics for arousal ratings in imagery-based conditioning group
describe(dataArousal[dataArousal$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = arousal rating
anovaArousalIma <- ezANOVA(
  data = dataArousalLong[dataArousalLong$usGroup == "ima" & dataArousalLong$time == "Post",],
  dv = arousal,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaArousalIma$ANOVA$pEtaSq <-
  c(anovaArousalIma$ANOVA$SSn[1] / (anovaArousalIma$ANOVA$SSd[1]+anovaArousalIma$ANOVA$SSn[1]),
    anovaArousalIma$ANOVA$SSn[2] / (anovaArousalIma$ANOVA$SSd[2]+anovaArousalIma$ANOVA$SSn[2])
  ); print(anovaArousalIma)

# bayesian ANOVA on arousal ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFArousalIma <- anovaBF(
  formula = arousal ~ CS + partInd,
  data = dataArousalLong[dataArousalLong$usGroup == "ima" & dataArousalLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFArousalIma)

# frequentist & bayesian t-tests on arousal ratings in imagery-based conditioning group
# CS+av vs CS+neu
arousalImaAvNeu_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                            y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
arousalImaAvNeu_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                              y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                              paired = TRUE)
arousalImaAvNeu_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                              y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalImaAvMin_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                            y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
arousalImaAvMin_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                              y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                              paired = TRUE)
arousalImaAvMin_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                              y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y

# CS+neu vs CS-
arousalImaNeuMin_t <- t.test(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                             y = dataArousal$Min_Post[dataArousal$usGroup == "ima"], 
                             alternative = "two.sided", paired = TRUE) # two-sided
arousalImaNeuMin_d <- cohens_d(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                               y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                               paired = TRUE)
arousalImaNeuMin_BF <- ttestBF(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                               y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided



#####################################################################
### Imagery-based conditioning - anger ratings - primary analyses ###
#####################################################################

# descriptive statistics for anger ratings in imagery-based conditioning group
describe(dataAnger[dataAnger$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = anger rating
anovaAngerIma <- ezANOVA(
  data = dataAngerLong[dataAngerLong$usGroup == "ima" & dataAngerLong$time == "Post",],
  dv = anger,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaAngerIma$ANOVA$pEtaSq <-
  c(anovaAngerIma$ANOVA$SSn[1] / (anovaAngerIma$ANOVA$SSd[1]+anovaAngerIma$ANOVA$SSn[1]),
    anovaAngerIma$ANOVA$SSn[2] / (anovaAngerIma$ANOVA$SSd[2]+anovaAngerIma$ANOVA$SSn[2])
  ); print(anovaAngerIma)

# bayesian ANOVA on anger ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFAngerIma <- anovaBF(
  formula = anger ~ CS + partInd,
  data = dataAngerLong[dataAngerLong$usGroup == "ima" & dataAngerLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFAngerIma)

# frequentist & bayesian t-tests on anger ratings in imagery-based conditioning group
# CS+av vs CS+neu
angerImaAvNeu_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                            y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
angerImaAvNeu_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                              y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                              paired = TRUE)
angerImaAvNeu_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                              y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
angerImaAvMin_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                            y = dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
angerImaAvMin_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                              y = dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                              paired = TRUE)
angerImaAvMin_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "ima"],
                              y = dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y

# CS+neu vs CS-
angerImaNeuMin_t <- t.test(x = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                             y = dataAnger$Min_Post[dataAnger$usGroup == "ima"], 
                             alternative = "two.sided", paired = TRUE) # two-sided
angerImaNeuMin_d <- cohens_d(x = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                               y = dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                               paired = TRUE)
angerImaNeuMin_BF <- ttestBF(x = dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                               y = dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided



#######################################################################
### Imagery-based conditioning - disgust ratings - primary analyses ###
#######################################################################

# descriptive statistics for disgust ratings in imagery-based conditioning group
describe(dataDisgust[dataDisgust$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = disgust rating
anovaDisgustIma <- ezANOVA(
  data = dataDisgustLong[dataDisgustLong$usGroup == "ima" & dataDisgustLong$time == "Post",],
  dv = disgust,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaDisgustIma$ANOVA$pEtaSq <-
  c(anovaDisgustIma$ANOVA$SSn[1] / (anovaDisgustIma$ANOVA$SSd[1]+anovaDisgustIma$ANOVA$SSn[1]),
    anovaDisgustIma$ANOVA$SSn[2] / (anovaDisgustIma$ANOVA$SSd[2]+anovaDisgustIma$ANOVA$SSn[2])
  ); print(anovaDisgustIma)

# bayesian ANOVA on disgust ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFDisgustIma <- anovaBF(
  formula = disgust ~ CS + partInd,
  data = dataDisgustLong[dataDisgustLong$usGroup == "ima" & dataDisgustLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFDisgustIma)

# frequentist & bayesian t-tests on disgust ratings in imagery-based conditioning group
# CS+av vs CS+neu
disgustImaAvNeu_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                          y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                          alternative = "greater", paired = TRUE) # one-sided
disgustImaAvNeu_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                            y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                            paired = TRUE)
disgustImaAvNeu_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                            y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
disgustImaAvMin_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                          y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                          alternative = "greater", paired = TRUE) # one-sided
disgustImaAvMin_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                            y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                            paired = TRUE)
disgustImaAvMin_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"],
                            y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y

# CS+neu vs CS-
disgustImaNeuMin_t <- t.test(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                           y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"], 
                           alternative = "two.sided", paired = TRUE) # two-sided
disgustImaNeuMin_d <- cohens_d(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                             y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                             paired = TRUE)
disgustImaNeuMin_BF <- ttestBF(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                             y = dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                             nullIntervall = NULL, paired = TRUE) # two-sided



##########################################################################
### Classical conditioning - unpleasantness ratings - primary analyses ###
##########################################################################

# descriptive statistics for unpleasantness ratings in classical conditioning group
describe(dataUnpleas[dataUnpleas$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = unpleasantness rating
anovaUnpleasReal <- ezANOVA(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "real" & dataUnpleasLong$time == "Post",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaUnpleasReal$ANOVA$pEtaSq <-
  c(anovaUnpleasReal$ANOVA$SSn[1] / (anovaUnpleasReal$ANOVA$SSd[1]+anovaUnpleasReal$ANOVA$SSn[1]),
    anovaUnpleasReal$ANOVA$SSn[2] / (anovaUnpleasReal$ANOVA$SSd[2]+anovaUnpleasReal$ANOVA$SSn[2])
  ); print(anovaUnpleasReal)

# bayesian ANOVA on unpleasantness ratings in classical conditioning group
set.seed(rngSeed); anovaBFUnpleasReal <- anovaBF(
  formula = unpleasantness ~ CS + partInd,
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "real" & dataUnpleasLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFUnpleasReal)

# frequentist & bayesian t-tests on unpleasantness ratings in classical conditioning group
# CS+av vs CS+neu
unpleasRealAvNeu_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                             y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvNeu_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                               paired = TRUE)
unpleasRealAvNeu_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasRealAvMin_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                             y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvMin_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                               paired = TRUE)
unpleasRealAvMin_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasRealNeuMin_t <- t.test(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                              y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                              alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealNeuMin_d <- cohens_d(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                paired = TRUE)
unpleasRealNeuMin_BF <- ttestBF(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                nullIntervall = NULL, paired = TRUE) # two-sided



###################################################################
### Classical conditioning - arousal ratings - primary analyses ###
###################################################################

# descriptive statistics for arousal ratings in classical conditioning group
describe(dataArousal[dataArousal$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = arousal rating
anovaArousalReal <- ezANOVA(
  data = dataArousalLong[dataArousalLong$usGroup == "real" & dataArousalLong$time == "Post",],
  dv = arousal,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaArousalReal$ANOVA$pEtaSq <-
  c(anovaArousalReal$ANOVA$SSn[1] / (anovaArousalReal$ANOVA$SSd[1]+anovaArousalReal$ANOVA$SSn[1]),
    anovaArousalReal$ANOVA$SSn[2] / (anovaArousalReal$ANOVA$SSd[2]+anovaArousalReal$ANOVA$SSn[2])
  ); print(anovaArousalReal)

# bayesian ANOVA on arousal ratings in classical conditioning group
set.seed(rngSeed); anovaBFArousalReal <- anovaBF(
  formula = arousal ~ CS + partInd,
  data = dataArousalLong[dataArousalLong$usGroup == "real" & dataArousalLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFArousalReal)

# frequentist & bayesian t-tests on arousal ratings in classical conditioning group
# CS+av vs CS+neu
arousalRealAvNeu_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                             y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
arousalRealAvNeu_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                               paired = TRUE)
arousalRealAvNeu_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalRealAvMin_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                              y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
arousalRealAvMin_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                               paired = TRUE)
arousalRealAvMin_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
arousalRealNeuMin_t <- t.test(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                              y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                              alternative = "two.sided", paired = TRUE) # two-sided
arousalRealNeuMin_d <- cohens_d(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Min_Post[dataArousal$usGroup == "real"], 
                                paired = TRUE)
arousalRealNeuMin_BF <- ttestBF(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                nullIntervall = NULL, paired = TRUE) # two-sided



#################################################################
### Classical conditioning - anger ratings - primary analyses ###
#################################################################

# descriptive statistics for anger ratings in imagery-based conditioning group
describe(dataAnger[dataAnger$usGroup == "real",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = anger rating
anovaAngerReal <- ezANOVA(
  data = dataAngerLong[dataAngerLong$usGroup == "real" & dataAngerLong$time == "Post",],
  dv = anger,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaAngerReal$ANOVA$pEtaSq <-
  c(anovaAngerReal$ANOVA$SSn[1] / (anovaAngerReal$ANOVA$SSd[1]+anovaAngerReal$ANOVA$SSn[1]),
    anovaAngerReal$ANOVA$SSn[2] / (anovaAngerReal$ANOVA$SSd[2]+anovaAngerReal$ANOVA$SSn[2])
  ); print(anovaAngerReal)

# bayesian ANOVA on anger ratings in classical conditioning group
set.seed(rngSeed); anovaBFAngerReal <- anovaBF(
  formula = anger ~ CS + partInd,
  data = dataAngerLong[dataAngerLong$usGroup == "real" & dataAngerLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFAngerReal)

# frequentist & bayesian t-tests on anger ratings in classical conditioning group
# CS+av vs CS+neu
angerRealAvNeu_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                          y = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                          alternative = "greater", paired = TRUE) # one-sided
angerRealAvNeu_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                            paired = TRUE)
angerRealAvNeu_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
angerRealAvMin_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                          y = dataAnger$Min_Post[dataAnger$usGroup == "real"],
                          alternative = "greater", paired = TRUE) # one-sided
angerRealAvMin_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Min_Post[dataAnger$usGroup == "real"],
                            paired = TRUE)
angerRealAvMin_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Min_Post[dataAnger$usGroup == "real"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y

# CS+neu vs CS-
angerRealNeuMin_t <- t.test(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                           y = dataAnger$Min_Post[dataAnger$usGroup == "real"], 
                           alternative = "two.sided", paired = TRUE) # two-sided
angerRealNeuMin_d <- cohens_d(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                             y = dataAnger$Min_Post[dataAnger$usGroup == "real"],
                             paired = TRUE)
angerRealNeuMin_BF <- ttestBF(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                             y = dataAnger$Min_Post[dataAnger$usGroup == "real"],
                             nullIntervall = NULL, paired = TRUE) # two-sided



###################################################################
### Classical conditioning - disgust ratings - primary analyses ###
###################################################################

# descriptive statistics for disgust ratings in classical conditioning group
describe(dataDisgust[dataDisgust$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = disgust rating
anovaDisgustReal <- ezANOVA(
  data = dataDisgustLong[dataDisgustLong$usGroup == "real" & dataDisgustLong$time == "Post",],
  dv = disgust,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaDisgustReal$ANOVA$pEtaSq <-
  c(anovaDisgustReal$ANOVA$SSn[1] / (anovaDisgustReal$ANOVA$SSd[1]+anovaDisgustReal$ANOVA$SSn[1]),
    anovaDisgustReal$ANOVA$SSn[2] / (anovaDisgustReal$ANOVA$SSd[2]+anovaDisgustReal$ANOVA$SSn[2])
  ); print(anovaDisgustReal)

# bayesian ANOVA on disgust ratings in classical conditioning group
set.seed(rngSeed); anovaBFDisgustReal <- anovaBF(
  formula = disgust ~ CS + partInd,
  data = dataDisgustLong[dataDisgustLong$usGroup == "real" & dataDisgustLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFDisgustReal)

# frequentist & bayesian t-tests on disgust ratings in classical conditioning group
# CS+av vs CS+neu
disgustRealAvNeu_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                            alternative = "greater", paired = TRUE) # one-sided
disgustRealAvNeu_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                              y = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                              paired = TRUE)
disgustRealAvNeu_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                              y = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
disgustRealAvMin_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                            alternative = "greater", paired = TRUE) # one-sided
disgustRealAvMin_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                              y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                              paired = TRUE)
disgustRealAvMin_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"],
                              y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y

# CS+neu vs CS-
disgustRealNeuMin_t <- t.test(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                             y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"], 
                             alternative = "two.sided", paired = TRUE) # two-sided
disgustRealNeuMin_d <- cohens_d(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                               y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                               paired = TRUE)
disgustRealNeuMin_BF <- ttestBF(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                               y = dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                               nullIntervall = NULL, paired = TRUE) # two-sided



#################################################################
### Across groups - unpleasantness ratings - primary analyses ###
#################################################################

# descriptive statistics  for unpleasantess ratings across conditioning groups
describe(dataUnpleas)

# frequentist ANOVA on unpleasantness ratings across conditioning groups
anovaUnpleas <- ezANOVA(
  data = dataUnpleasLong[dataUnpleasLong$time == "Post",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaUnpleas$ANOVA$pEtaSq <- c(anovaUnpleas$ANOVA$SSn[1] /
                                    (anovaUnpleas$ANOVA$SSd[1]+anovaUnpleas$ANOVA$SSn[1]),
                                  anovaUnpleas$ANOVA$SSn[2] /
                                    (anovaUnpleas$ANOVA$SSd[2]+anovaUnpleas$ANOVA$SSn[2]),
                                  anovaUnpleas$ANOVA$SSn[3] /
                                    (anovaUnpleas$ANOVA$SSd[3]+anovaUnpleas$ANOVA$SSn[3]),
                                  anovaUnpleas$ANOVA$SSn[4] /
                                    (anovaUnpleas$ANOVA$SSd[4]+anovaUnpleas$ANOVA$SSn[4])
); print(anovaUnpleas)

# bayesian ANOVA on unpleasantness ratings across conditioning groups
set.seed(rngSeed); anovaBFUnpleas <- anovaBF(
  formula = unpleasantness ~ usGroup*CS + partInd,
  data = dataUnpleasLong[dataUnpleasLong$time == "Post",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFUnpleas)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFUnpleas)

# quick graph of group x CS ANOVA on unpleasantness ratings
ezPlot(
  data = dataUnpleasLong[dataUnpleasLong$time == "Post",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on unpleasantness ratings (difference scores) across groups
# delta [CS+av - CS+neu]
unpleasBothAvNeu_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                 dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                             y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                 dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
unpleasBothAvNeu_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                   dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                   dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               paired = FALSE)
unpleasBothAvNeu_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                   dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                   dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
unpleasBothAvMin_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                 dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                             y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                 dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
unpleasBothAvMin_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                   dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                   dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               paired = FALSE)
unpleasBothAvMin_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"] -
                                   dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"] -
                                   dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
unpleasBothNeuMin_t <- t.test(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"] -
                                  dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                              y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"] -
                                  dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
unpleasBothNeuMin_d <- cohens_d(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"] -
                                    dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"] -
                                    dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                                paired = FALSE)
unpleasBothNeuMin_BF <- ttestBF(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"] -
                                    dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"] -
                                    dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided



##########################################################
### Across groups - arousal ratings - primary analyses ###
##########################################################

# descriptive statistics for arousal ratings across conditioning groups
describe(dataArousal)

# frequentist ANOVA on arousal ratings across conditioning groups
anovaArousal <- ezANOVA(
  data = dataArousalLong[dataArousalLong$time == "Post",],
  dv = arousal,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaArousal$ANOVA$pEtaSq <- c(anovaArousal$ANOVA$SSn[1] /
                                    (anovaArousal$ANOVA$SSd[1]+anovaArousal$ANOVA$SSn[1]),
                                  anovaArousal$ANOVA$SSn[2] /
                                    (anovaArousal$ANOVA$SSd[2]+anovaArousal$ANOVA$SSn[2]),
                                  anovaArousal$ANOVA$SSn[3] /
                                    (anovaArousal$ANOVA$SSd[3]+anovaArousal$ANOVA$SSn[3]),
                                  anovaArousal$ANOVA$SSn[4] /
                                    (anovaArousal$ANOVA$SSd[4]+anovaArousal$ANOVA$SSn[4])
); print(anovaArousal)

# bayesian ANOVA on arousal ratings across conditioning groups
set.seed(rngSeed); anovaBFArousal <- anovaBF(
  formula = arousal ~ usGroup*CS + partInd,
  data = dataArousalLong[dataArousalLong$time == "Post",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFArousal)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFArousal)

# quick graph of group x CS ANOVA on arousal ratings
ezPlot(
  data = dataArousalLong[dataArousalLong$time == "Post",],
  dv = arousal,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on arousal ratings (difference scores) across groups
# delta [CS+av - CS+neu]
arousalBothAvNeu_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                 dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                             y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                 dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
arousalBothAvNeu_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                   dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                   dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                               paired = FALSE)
arousalBothAvNeu_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                   dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                   dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
arousalBothAvMin_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                 dataArousal$Min_Post[dataArousal$usGroup == "real"],
                             y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                 dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
arousalBothAvMin_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                   dataArousal$Min_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                   dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                               paired = FALSE)
arousalBothAvMin_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"] -
                                   dataArousal$Min_Post[dataArousal$usGroup == "real"],
                               y = dataArousal$Av_Post[dataArousal$usGroup == "ima"] -
                                   dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
arousalBothNeuMin_t <- t.test(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"] -
                                  dataArousal$Min_Post[dataArousal$usGroup == "real"],
                              y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"] -
                                  dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
arousalBothNeuMin_d <- cohens_d(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"] -
                                    dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"] -
                                    dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                paired = FALSE)
arousalBothNeuMin_BF <- ttestBF(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"] -
                                    dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"] -
                                    dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided



#######################################################
### Across groups - anger ratings - primary analyses ###
#######################################################

# descriptive statistics  for anger ratings across conditioning groups
describe(dataAnger)

# frequentist ANOVA on anger ratings across conditioning groups
anovaAnger <- ezANOVA(
  data = dataAngerLong[dataAngerLong$time == "Post",],
  dv = anger,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaAnger$ANOVA$pEtaSq <- c(anovaAnger$ANOVA$SSn[1] /
                                 (anovaAnger$ANOVA$SSd[1]+anovaAnger$ANOVA$SSn[1]),
                               anovaAnger$ANOVA$SSn[2] /
                                 (anovaAnger$ANOVA$SSd[2]+anovaAnger$ANOVA$SSn[2]),
                               anovaAnger$ANOVA$SSn[3] /
                                 (anovaAnger$ANOVA$SSd[3]+anovaAnger$ANOVA$SSn[3]),
                               anovaAnger$ANOVA$SSn[4] /
                                 (anovaAnger$ANOVA$SSd[4]+anovaAnger$ANOVA$SSn[4])
); print(anovaAnger)

# bayesian ANOVA on anger ratings across conditioning groups
set.seed(rngSeed); anovaBFAnger <- anovaBF(
  formula = anger ~ usGroup*CS + partInd,
  data = dataAngerLong[dataAngerLong$time == "Post",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFAnger)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFAnger)

# quick graph of group x CS ANOVA on anger ratings
ezPlot(
  data = dataAngerLong[dataAngerLong$time == "Post",],
  dv = anger,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on anger ratings (difference scores) across groups
# delta [CS+av - CS+neu]
angerBothAvNeu_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                            dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                          y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                            dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
angerBothAvNeu_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                              dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                              dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                            paired = FALSE)
angerBothAvNeu_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                              dataAnger$Neu_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                              dataAnger$Neu_Post[dataAnger$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
angerBothAvMin_t <- t.test(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                            dataAnger$Min_Post[dataAnger$usGroup == "real"],
                          y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                            dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
angerBothAvMin_d <- cohens_d(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                              dataAnger$Min_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                              dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                            paired = FALSE)
angerBothAvMin_BF <- ttestBF(x = dataAnger$Av_Post[dataAnger$usGroup == "real"] -
                              dataAnger$Min_Post[dataAnger$usGroup == "real"],
                            y = dataAnger$Av_Post[dataAnger$usGroup == "ima"] -
                              dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
angerBothNeuMin_t <- t.test(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"] -
                             dataAnger$Min_Post[dataAnger$usGroup == "real"],
                           y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"] - 
                             dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                           alternative = "two.sided", paired = FALSE) # two-sided
angerBothNeuMin_d <- cohens_d(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"] -
                               dataAnger$Min_Post[dataAnger$usGroup == "real"],
                             y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"] -
                               dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                             paired = FALSE)
angerBothNeuMin_BF <- ttestBF(x = dataAnger$Neu_Post[dataAnger$usGroup == "real"] -
                               dataAnger$Min_Post[dataAnger$usGroup == "real"],
                             y = dataAnger$Neu_Post[dataAnger$usGroup == "ima"] - 
                               dataAnger$Min_Post[dataAnger$usGroup == "ima"],
                             nullInterval = NULL, paired = FALSE) # two-sided



##########################################################
### Across groups - disgust ratings - primary analyses ###
##########################################################

# descriptive statistics  for disgust ratings across conditioning groups
describe(dataDisgust)

# frequentist ANOVA on disgust ratings across conditioning groups
anovaDisgust <- ezANOVA(
  data = dataDisgustLong[dataDisgustLong$time == "Post",],
  dv = disgust,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaDisgust$ANOVA$pEtaSq <- c(anovaDisgust$ANOVA$SSn[1] /
                                 (anovaDisgust$ANOVA$SSd[1]+anovaDisgust$ANOVA$SSn[1]),
                               anovaDisgust$ANOVA$SSn[2] /
                                 (anovaDisgust$ANOVA$SSd[2]+anovaDisgust$ANOVA$SSn[2]),
                               anovaDisgust$ANOVA$SSn[3] /
                                 (anovaDisgust$ANOVA$SSd[3]+anovaDisgust$ANOVA$SSn[3]),
                               anovaDisgust$ANOVA$SSn[4] /
                                 (anovaDisgust$ANOVA$SSd[4]+anovaDisgust$ANOVA$SSn[4])
); print(anovaDisgust)

# bayesian ANOVA on disgust ratings across conditioning groups
set.seed(rngSeed); anovaBFDisgust <- anovaBF(
  formula = disgust ~ usGroup*CS + partInd,
  data = dataDisgustLong[dataDisgustLong$time == "Post",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFDisgust)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFDisgust)

# quick graph of group x CS ANOVA on disgust ratings
ezPlot(
  data = dataDisgustLong[dataDisgustLong$time == "Post",],
  dv = disgust,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on disgust ratings (difference scores) across groups
# delta [CS+av - CS+neu]
disgustBothAvNeu_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                            dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                          y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                            dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
disgustBothAvNeu_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                              dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                              dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                            paired = FALSE)
disgustBothAvNeu_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                              dataDisgust$Neu_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                              dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
disgustBothAvMin_t <- t.test(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                            dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                          y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                            dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
disgustBothAvMin_d <- cohens_d(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                              dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                              dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                            paired = FALSE)
disgustBothAvMin_BF <- ttestBF(x = dataDisgust$Av_Post[dataDisgust$usGroup == "real"] -
                              dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                            y = dataDisgust$Av_Post[dataDisgust$usGroup == "ima"] -
                              dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
disgustBothNeuMin_t <- t.test(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"] -
                             dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                           y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"] - 
                             dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                           alternative = "two.sided", paired = FALSE) # two-sided
disgustBothNeuMin_d <- cohens_d(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"] -
                               dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                             y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"] -
                               dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                             paired = FALSE)
disgustBothNeuMin_BF <- ttestBF(x = dataDisgust$Neu_Post[dataDisgust$usGroup == "real"] -
                               dataDisgust$Min_Post[dataDisgust$usGroup == "real"],
                             y = dataDisgust$Neu_Post[dataDisgust$usGroup == "ima"] - 
                               dataDisgust$Min_Post[dataDisgust$usGroup == "ima"],
                             nullInterval = NULL, paired = FALSE) # two-sided


#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-",
                     "groups: delta CS+av / CS+neu", "groups: delta CS+av / CS-", "groups: delta CSneu / CS-"), 4),
  t = c(unpleasImaAvNeu_t$statistic, unpleasImaAvMin_t$statistic, unpleasImaNeuMin_t$statistic,
        unpleasRealAvNeu_t$statistic, unpleasRealAvMin_t$statistic, unpleasRealNeuMin_t$statistic,
        unpleasBothAvNeu_t$statistic, unpleasBothAvMin_t$statistic, unpleasBothNeuMin_t$statistic,
        arousalImaAvNeu_t$statistic, arousalImaAvMin_t$statistic, arousalImaNeuMin_t$statistic,
        arousalRealAvNeu_t$statistic, arousalRealAvMin_t$statistic, arousalRealNeuMin_t$statistic,
        arousalBothAvNeu_t$statistic, arousalBothAvMin_t$statistic, arousalBothNeuMin_t$statistic,
        angerImaAvNeu_t$statistic, angerImaAvMin_t$statistic, angerImaNeuMin_t$statistic,
        angerRealAvNeu_t$statistic, angerRealAvMin_t$statistic, angerRealNeuMin_t$statistic,
        angerBothAvNeu_t$statistic, angerBothAvMin_t$statistic, angerBothNeuMin_t$statistic,
        disgustImaAvNeu_t$statistic, disgustImaAvMin_t$statistic, disgustImaNeuMin_t$statistic,
        disgustRealAvNeu_t$statistic, disgustRealAvMin_t$statistic, disgustRealNeuMin_t$statistic,
        disgustBothAvNeu_t$statistic, disgustBothAvMin_t$statistic, disgustBothNeuMin_t$statistic), 
  df = c(unpleasImaAvNeu_t$parameter, unpleasImaAvMin_t$parameter, unpleasImaNeuMin_t$parameter,
         unpleasRealAvNeu_t$parameter, unpleasRealAvMin_t$parameter, unpleasRealNeuMin_t$parameter,
         unpleasBothAvNeu_t$parameter, unpleasBothAvMin_t$parameter, unpleasBothNeuMin_t$parameter,
         arousalImaAvNeu_t$parameter, arousalImaAvMin_t$parameter, arousalImaNeuMin_t$parameter,
         arousalRealAvNeu_t$parameter, arousalRealAvMin_t$parameter, arousalRealNeuMin_t$parameter,
         arousalBothAvNeu_t$parameter, arousalBothAvMin_t$parameter, arousalBothNeuMin_t$parameter,
         angerImaAvNeu_t$parameter, angerImaAvMin_t$parameter, angerImaNeuMin_t$parameter,
         angerRealAvNeu_t$parameter, angerRealAvMin_t$parameter, angerRealNeuMin_t$parameter,
         angerBothAvNeu_t$parameter, angerBothAvMin_t$parameter, angerBothNeuMin_t$parameter,
         disgustImaAvNeu_t$parameter, disgustImaAvMin_t$parameter, disgustImaNeuMin_t$parameter,
         disgustRealAvNeu_t$parameter, disgustRealAvMin_t$parameter, disgustRealNeuMin_t$parameter,
         disgustBothAvNeu_t$parameter, disgustBothAvMin_t$parameter, disgustBothNeuMin_t$parameter), 
  p = c(unpleasImaAvNeu_t$p.value, unpleasImaAvMin_t$p.value, unpleasImaNeuMin_t$p.value,
        unpleasRealAvNeu_t$p.value, unpleasRealAvMin_t$p.value, unpleasRealNeuMin_t$p.value,
        unpleasBothAvNeu_t$p.value*3, unpleasBothAvMin_t$p.value*3, unpleasBothNeuMin_t$p.value*3, # Bonferroni
        arousalImaAvNeu_t$p.value, arousalImaAvMin_t$p.value, arousalImaNeuMin_t$p.value,
        arousalRealAvNeu_t$p.value, arousalRealAvMin_t$p.value, arousalRealNeuMin_t$p.value,
        arousalBothAvNeu_t$p.value*3, arousalBothAvMin_t$p.value*3, arousalBothNeuMin_t$p.value*3,  # Bonferroni
        angerImaAvNeu_t$p.value, angerImaAvMin_t$p.value, angerImaNeuMin_t$p.value,
        angerRealAvNeu_t$p.value, angerRealAvMin_t$p.value, angerRealNeuMin_t$p.value,
        angerBothAvNeu_t$p.value*3, angerBothAvMin_t$p.value*3, angerBothNeuMin_t$p.value*3,  # Bonferroni
        disgustImaAvNeu_t$p.value, disgustImaAvMin_t$p.value, disgustImaNeuMin_t$p.value,
        disgustRealAvNeu_t$p.value, disgustRealAvMin_t$p.value, disgustRealNeuMin_t$p.value,
        disgustBothAvNeu_t$p.value*3, disgustBothAvMin_t$p.value*3, disgustBothNeuMin_t$p.value*3),  # Bonferroni
  d = c(unpleasImaAvNeu_d$Cohens_d, unpleasImaAvMin_d$Cohens_d, unpleasImaNeuMin_d$Cohens_d,
        unpleasRealAvNeu_d$Cohens_d, unpleasRealAvMin_d$Cohens_d, unpleasRealNeuMin_d$Cohens_d,
        unpleasBothAvNeu_d$Cohens_d, unpleasBothAvMin_d$Cohens_d, unpleasBothNeuMin_d$Cohens_d,
        arousalImaAvNeu_d$Cohens_d, arousalImaAvMin_d$Cohens_d, arousalImaNeuMin_d$Cohens_d,
        arousalRealAvNeu_d$Cohens_d, arousalRealAvMin_d$Cohens_d, arousalRealNeuMin_d$Cohens_d,
        arousalBothAvNeu_d$Cohens_d, arousalBothAvMin_d$Cohens_d, arousalBothNeuMin_d$Cohens_d,
        angerImaAvNeu_d$Cohens_d, angerImaAvMin_d$Cohens_d, angerImaNeuMin_d$Cohens_d,
        angerRealAvNeu_d$Cohens_d, angerRealAvMin_d$Cohens_d, angerRealNeuMin_d$Cohens_d,
        angerBothAvNeu_d$Cohens_d, angerBothAvMin_d$Cohens_d, angerBothNeuMin_d$Cohens_d,
        disgustImaAvNeu_d$Cohens_d, disgustImaAvMin_d$Cohens_d, disgustImaNeuMin_d$Cohens_d,
        disgustRealAvNeu_d$Cohens_d, disgustRealAvMin_d$Cohens_d, disgustRealNeuMin_d$Cohens_d,
        disgustBothAvNeu_d$Cohens_d, disgustBothAvMin_d$Cohens_d, disgustBothNeuMin_d$Cohens_d), 
  BF = c(exp(unpleasImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(unpleasImaAvMin_BF@bayesFactor[["bf"]][1]), exp(unpleasImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(unpleasRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(unpleasRealAvMin_BF@bayesFactor[["bf"]][1]), exp(unpleasRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(unpleasBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(unpleasBothAvMin_BF@bayesFactor[["bf"]][1]), exp(unpleasBothNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(arousalImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(arousalImaAvMin_BF@bayesFactor[["bf"]][1]), exp(arousalImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(arousalRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(arousalRealAvMin_BF@bayesFactor[["bf"]][1]), exp(arousalRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(arousalBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(arousalBothAvMin_BF@bayesFactor[["bf"]][1]), exp(arousalBothNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(angerImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(angerImaAvMin_BF@bayesFactor[["bf"]][1]), exp(angerImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(angerRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(angerRealAvMin_BF@bayesFactor[["bf"]][1]), exp(angerRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(angerBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(angerBothAvMin_BF@bayesFactor[["bf"]][1]), exp(angerBothNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(disgustImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(disgustImaAvMin_BF@bayesFactor[["bf"]][1]), exp(disgustImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(disgustRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(disgustRealAvMin_BF@bayesFactor[["bf"]][1]), exp(disgustRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(disgustBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(disgustBothAvMin_BF@bayesFactor[["bf"]][1]), exp(disgustBothNeuMin_BF@bayesFactor[["bf"]][1]))
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

tableUnpleas <- flextable(tableData[1:9,])
tableUnpleas <- add_header_lines(tableUnpleas, top = TRUE, values = "unpleasantness")
tableUnpleas <- align(tableUnpleas, align = "center")
tableArousal <- flextable(tableData[10:18,])
tableArousal <- add_header_lines(tableArousal, top = TRUE, values = "arousal")
tableAnger <- flextable(tableData[19:27,])
tableAnger <- add_header_lines(tableAnger, top = TRUE, values = "anger")
tableDisgust <- flextable(tableData[28:36,])
tableDisgust <- add_header_lines(tableDisgust, top = TRUE, values = "disgust")

save_as_docx(tableUnpleas, path = paste0(pathname, "/supplement/t_table_Unpleasant.docx"))
save_as_docx(tableArousal, path = paste0(pathname, "/supplement/t_table_Arousal.docx"))
save_as_docx(tableAnger, path = paste0(pathname, "/supplement/t_table_Anger.docx"))
save_as_docx(tableDisgust, path = paste0(pathname, "/supplement/t_table_Disgust.docx"))



#############################
### Plots for all ratings ###
#############################
# remove between-subject variance for plotting standard errors based on
# within-subject variance
dataUnpleasLongPost <- dataUnpleas[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]
dataArousalLongPost <- dataArousal[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]
dataAngerLongPost <- dataAnger[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]
dataDisgustLongPost <- dataDisgust[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]

# remove each participant's average from each single value
dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "ima",6:8] <- as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "ima",3:5])) + mean(as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "ima",3:5]))
dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "real",6:8] <- as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "real",3:5])) + mean(as.matrix(dataUnpleasLongPost[dataUnpleasLongPost$usGroup == "real",3:5]))
names(dataUnpleasLongPost) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")

dataArousalLongPost[dataArousalLongPost$usGroup == "ima",6:8] <- as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "ima",3:5])) + mean(as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "ima",3:5]))
dataArousalLongPost[dataArousalLongPost$usGroup == "real",6:8] <- as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "real",3:5])) + mean(as.matrix(dataArousalLongPost[dataArousalLongPost$usGroup == "real",3:5]))
names(dataArousalLongPost) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")

dataAngerLongPost[dataAngerLongPost$usGroup == "ima",6:8] <- as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "ima",3:5])) + mean(as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "ima",3:5]))
dataAngerLongPost[dataAngerLongPost$usGroup == "real",6:8] <- as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "real",3:5])) + mean(as.matrix(dataAngerLongPost[dataAngerLongPost$usGroup == "real",3:5]))
names(dataAngerLongPost) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")

dataDisgustLongPost[dataDisgustLongPost$usGroup == "ima",6:8] <- as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "ima",3:5])) + mean(as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "ima",3:5]))
dataDisgustLongPost[dataDisgustLongPost$usGroup == "real",6:8] <- as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "real",3:5])) + mean(as.matrix(dataDisgustLongPost[dataDisgustLongPost$usGroup == "real",3:5]))
names(dataDisgustLongPost) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")

# into long format
dataUnpleasLongPost <- pivot_longer(data = dataUnpleasLongPost, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "fear")
dataUnpleasLongPost <- pivot_wider(data = dataUnpleasLongPost, names_from = "variance", values_from = "fear")
dataUnpleasLongPost$CS <- factor(dataUnpleasLongPost$CS, levels = c("Av","Neu","Min"))
levels(dataUnpleasLongPost$usGroup) <- c("Imagery-based","Classical")

dataArousalLongPost <- pivot_longer(data = dataArousalLongPost, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "fear")
dataArousalLongPost <- pivot_wider(data = dataArousalLongPost, names_from = "variance", values_from = "fear")
dataArousalLongPost$CS <- factor(dataArousalLongPost$CS, levels = c("Av","Neu","Min"))
levels(dataArousalLongPost$usGroup) <- c("Imagery-based","Classical")

dataAngerLongPost <- pivot_longer(data = dataAngerLongPost, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "fear")
dataAngerLongPost <- pivot_wider(data = dataAngerLongPost, names_from = "variance", values_from = "fear")
dataAngerLongPost$CS <- factor(dataAngerLongPost$CS, levels = c("Av","Neu","Min"))
levels(dataAngerLongPost$usGroup) <- c("Imagery-based","Classical")

dataDisgustLongPost <- pivot_longer(data = dataDisgustLongPost, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "fear")
dataDisgustLongPost <- pivot_wider(data = dataDisgustLongPost, names_from = "variance", values_from = "fear")
dataDisgustLongPost$CS <- factor(dataDisgustLongPost$CS, levels = c("Av","Neu","Min"))
levels(dataDisgustLongPost$usGroup) <- c("Imagery-based","Classical")


# some general settings
plotFS <- 8
showSig <- TRUE
csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")

# graphs of group x CS effects on unpleasantness ratings
graphUnpleas <- ggplot(data = dataUnpleasLongPost, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.1, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.25, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Unpleasantness rating (1-5)", limits = c(0.4,6), breaks = 1:5, oob = rescale_none, expand = c(0,0)) +
  geom_vline(xintercept = 0.41) +
  geom_rect(xmin = 0.4, xmax = 2.6, ymin = 0.35, ymax = 0.8, fill = "white", color = "white") +
  geom_hline(yintercept = 0.8) +
  geom_beeswarm(aes(color = CS), dodge.width = 0.8, cex = 0.6, size = .1, color = "gray70") +
  geom_violin(alpha = .2, color = NA, bw = .5, position = position_dodge(0.8), width = 0.75) +
  geom_text(x = 0.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(aes(label = usGroup, y = 5.9), colour = "black", size = (plotFS-2)/.pt, fontface = "bold") +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y  = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

if (showSig == TRUE){
  graphUnpleas <- graphUnpleas +
    geom_segment(aes(x = 0.74, y = 5.1, xend = 1.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 0.87, y =5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 0.74, y = 5.4, xend = 1.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.1, xend = 2.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.87, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.4, xend = 2.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 2.0, y = 5.45), size = plotFS/4, color = "gray20")
}
graphUnpleas

# graphs of group x CS effects on arousal ratings
graphArousal <- ggplot(data = dataArousalLongPost, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.1, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.25, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Arousal rating (1-5)", limits = c(0.4,6), breaks = 1:5, oob = rescale_none, expand = c(0,0)) +
  geom_vline(xintercept = 0.41) +
  geom_rect(xmin = 0.4, xmax = 2.6, ymin = 0.35, ymax = 0.8, fill = "white", color = "white") +
  geom_hline(yintercept = 0.8) +
  geom_beeswarm(aes(color = CS), dodge.width = 0.8, cex = 0.6, size = .1, color = "gray70") +
  geom_violin(alpha = .2, color = NA, bw = .5, position = position_dodge(0.8), width = 0.75) +
  geom_text(x = 0.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(aes(label = usGroup, y = 5.9), colour = "black", size = (plotFS-2)/.pt, fontface = "bold") +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y  = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

if (showSig == TRUE){
  graphArousal <- graphArousal +
    geom_segment(aes(x = 1.0, y = 5.1, xend = 1.26, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "**", x = 1.13, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 0.74, y = 5.4, xend = 1.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.1, xend = 1.99, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.87, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.4, xend = 2.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 2.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 2.01, y = 5.1, xend = 2.26, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "**", x = 2.13, y = 5.15), size = plotFS/4, color = "gray20")
}
graphArousal

# graphs of group x CS effects on anger ratings
graphAnger <- ggplot(data = dataAngerLongPost, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.1, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.25, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Anger rating (1-5)", limits = c(0.4,6), breaks = 1:5, oob = rescale_none, expand = c(0,0)) +
  geom_vline(xintercept = 0.41) +
  geom_rect(xmin = 0.4, xmax = 2.6, ymin = 0.35, ymax = 0.8, fill = "white", color = "white") +
  geom_hline(yintercept = 0.8) +
  geom_beeswarm(aes(color = CS), dodge.width = 0.8, cex = 0.6, size = .1, color = "gray70") +
  geom_violin(alpha = .2, color = NA, bw = .5, position = position_dodge(0.8), width = 0.75) +
  geom_text(x = 0.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(aes(label = usGroup, y = 5.9), colour = "black", size = (plotFS-2)/.pt, fontface = "bold") +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y  = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

if (showSig == TRUE){
  graphAnger <- graphAnger +
    geom_segment(aes(x = 0.74, y = 5.1, xend = 1.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 0.87, y =5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 0.74, y = 5.4, xend = 1.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.1, xend = 2.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.87, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.4, xend = 2.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 2.0, y = 5.45), size = plotFS/4, color = "gray20")
}
graphAnger

# graphs of group x CS effects on disgust ratings
graphDisgust <- ggplot(data = dataDisgustLongPost, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.1, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.25, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Disgust rating (1-5)", limits = c(0.4,6), breaks = 1:5, oob = rescale_none, expand = c(0,0)) +
  geom_vline(xintercept = 0.41) +
  geom_rect(xmin = 0.4, xmax = 2.6, ymin = 0.35, ymax = 0.8, fill = "white", color = "white") +
  geom_hline(yintercept = 0.8) +
  geom_beeswarm(aes(color = CS), dodge.width = 0.8, cex = 0.6, size = .1, color = "gray70") +
  geom_violin(alpha = .2, color = NA, bw = .5, position = position_dodge(0.8), width = 0.75) +
  geom_text(x = 0.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.74, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.26, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(aes(label = usGroup, y = 5.9), colour = "black", size = (plotFS-2)/.pt, fontface = "bold") +
  theme(legend.position = "none",
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y  = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

if (showSig == TRUE){
  graphDisgust <- graphDisgust +
    geom_segment(aes(x = 0.74, y = 5.1, xend = 1.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 0.87, y =5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 0.74, y = 5.4, xend = 1.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.1, xend = 2.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "**", x = 1.87, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.4, xend = 2.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "**", x = 2.0, y = 5.45), size = plotFS/4, color = "gray20")
}
graphDisgust



### combining graphs into one figure
# add margins to subplots
graphUnpleas <- graphUnpleas + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphArousal <- graphArousal + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphAnger <- graphAnger + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphDisgust <- graphDisgust + theme(plot.margin = unit(c(5,5,5,5), "mm"))

# merging the subplots
graphRatings <- ggarrange(graphUnpleas, graphArousal,
                          graphAnger, graphDisgust,
                          labels = c("A", "B", "C", "D"),
                          ncol = 2, nrow = 2)



# saving it
ggsave(filename = paste0(pathname, "/supplement/FigureSupp_plotOtherRatings.eps"),
       plot = graphRatings,
       width = 210,
       height = 180,
       units = "mm",
       dpi = 300
)

ggsave(filename = paste0(pathname, "/supplement/FigureSupp_plotOtherRatings.pdf"),
       plot = graphRatings,
       device = "pdf",
       width = 210,
       height = 180,
       units = "mm",
       dpi = 300
)
