# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Jan 2023
# --- content: supplementary rating analyses - ANOVA with factor Time

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

########################
### data preparation ###
########################

# load rating data from text file
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# create data frames in wide & long format for fear ratings
dataFear <- data.frame(
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
  Av_Pre = importRatings$anx_csplus_av_2,
  Av_Mid = importRatings$anx_csplus_av_3,
  Av_Post = importRatings$anx_csplus_av_4,
  Neu_Pre = importRatings$anx_csplus_neu_2,
  Neu_Mid = importRatings$anx_csplus_neu_3,
  Neu_Post = importRatings$anx_csplus_neu_4,
  Min_Pre = importRatings$anx_csminus_2,
  Min_Mid = importRatings$anx_csminus_3,
  Min_Post = importRatings$anx_csminus_4
)  
dataFearLong <- gather(data = dataFear, key = "cond", value = "fear",
                       Av_Pre:Min_Post)
dataFearLong <- separate(data = dataFearLong, col = cond, into = c("CS","time"),
                         sep = "_")
dataFearLong$CS <- factor(dataFearLong$CS, levels = c("Av", "Neu", "Min"))
dataFearLong$time <- factor(dataFearLong$time, levels = c("Pre", "Mid", "Post"))

# create data frames in wide & long format for unpleasantness ratings
dataUnpleas <- data.frame(
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
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
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
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

# create data frames in wide & long format for arousal ratings
dataAnger <- data.frame(
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
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
dataAngerLong <- separate(data = dataAngerLong, col = cond,
                            into = c("CS","time"), sep = "_")
dataAngerLong$CS <- factor(dataAngerLong$CS, levels = c("Av", "Neu", "Min"))
dataAngerLong$time <- factor(dataAngerLong$time, levels = c("Pre", "Mid", "Post"))

# create data frames in wide & long format for fear ratings
dataDisgust <- data.frame(
  partInd = factor(1:dim(importRatings)[1]),
  usGroup = factor(importRatings$group, labels = c("ima", "real")),
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



###################################################################################
### Imagery-based conditioning - fear ratings - supplemental analyses over time ###
###################################################################################

# descriptive statistics for fear ratings in imagery-based conditioning group
describe(dataFear[dataFear$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = fear rating
anovaFearIma <- ezANOVA(
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  dv = fear,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaFearIma$ANOVA$pEtaSq <-
  c(anovaFearIma$ANOVA$SSn[1] / (anovaFearIma$ANOVA$SSd[1]+anovaFearIma$ANOVA$SSn[1]),
    anovaFearIma$ANOVA$SSn[2] / (anovaFearIma$ANOVA$SSd[2]+anovaFearIma$ANOVA$SSn[2]),
    anovaFearIma$ANOVA$SSn[3] / (anovaFearIma$ANOVA$SSd[3]+anovaFearIma$ANOVA$SSn[3]),
    anovaFearIma$ANOVA$SSn[4] / (anovaFearIma$ANOVA$SSd[4]+anovaFearIma$ANOVA$SSn[4])
  ); print(anovaFearIma)
capture.output(print(anovaFearIma), file = "supplement/01s1_timeFactor_fear_ima_anovaFreq.doc")

# bayesian ANOVA on fear ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFFearIma <- anovaBF(
  formula = fear ~ CS*time + partInd,
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFearIma)
capture.output(print(anovaBFFearIma), file = "supplement/01s1_timeFactor_fear_ima_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for fear ratings in imagery-based conditioning group
plotFearIma <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearIma 
ggsave(plot = plotFearIma, filename = "supplement/01s1_timeFactor_fear_ima_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on fear ratings in imagery-based conditioning group
### Pre
# CS+av vs CS+neu
fearImaAvNeuPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                            y = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                            alternative = "two.sided", paired = TRUE) # two-sided
fearImaAvNeuPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                              y = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                              paired = TRUE)
fearImaAvNeuPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                              y = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
fearImaAvMinPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                            y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                            alternative = "two.sided", paired = TRUE) # two-sided
fearImaAvMinPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                              y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                              paired = TRUE)
fearImaAvMinPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "ima"],
                              y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
fearImaNeuMinPre_t <- t.test(x = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                             y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = TRUE) # two-sided
fearImaNeuMinPre_d <- cohens_d(x = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                               paired = TRUE)
fearImaNeuMinPre_BF <- ttestBF(x = dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Pre[dataFear$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided
### Mid
# CS+av vs CS+neu
fearImaAvNeuMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                            y = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
fearImaAvNeuMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                              y = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                              paired = TRUE)
fearImaAvNeuMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                              y = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearImaAvMinMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                            y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
fearImaAvMinMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                              y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                              paired = TRUE)
fearImaAvMinMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "ima"],
                              y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearImaNeuMinMid_t <- t.test(x = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                             y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = TRUE) # two-sided
fearImaNeuMinMid_d <- cohens_d(x = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                               paired = TRUE)
fearImaNeuMinMid_BF <- ttestBF(x = dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Mid[dataFear$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
fearImaAvNeuPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                             y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
fearImaAvNeuPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                               y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                               paired = TRUE)
fearImaAvNeuPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                               y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearImaAvMinPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                             y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
fearImaAvMinPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                               paired = TRUE)
fearImaAvMinPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                               y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearImaNeuMinPost_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                              y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
fearImaNeuMinPost_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                                y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                                paired = TRUE)
fearImaNeuMinPost_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                                y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                                nullIntervall = NULL, paired = TRUE) # two-sided

tableFearIma <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(fearImaAvNeuPre_t$statistic, fearImaAvMinPre_t$statistic, fearImaNeuMinPre_t$statistic,
        fearImaAvNeuMid_t$statistic, fearImaAvMinMid_t$statistic, fearImaNeuMinMid_t$statistic,
        fearImaAvNeuPost_t$statistic, fearImaAvMinPost_t$statistic, fearImaNeuMinPost_t$statistic),
  df = c(fearImaAvNeuPre_t$parameter, fearImaAvMinPre_t$parameter, fearImaNeuMinPre_t$parameter,
         fearImaAvNeuMid_t$parameter, fearImaAvMinMid_t$parameter, fearImaNeuMinMid_t$parameter,
         fearImaAvNeuPost_t$parameter, fearImaAvMinPost_t$parameter, fearImaNeuMinPost_t$parameter), 
  p = c(fearImaAvNeuPre_t$p.value, fearImaAvMinPre_t$p.value, fearImaNeuMinPre_t$p.value,
        fearImaAvNeuMid_t$p.value, fearImaAvMinMid_t$p.value, fearImaNeuMinMid_t$p.value,
        fearImaAvNeuPost_t$p.value, fearImaAvMinPost_t$p.value, fearImaNeuMinPost_t$p.value),
  d = c(fearImaAvNeuPre_d$Cohens_d, fearImaAvMinPre_d$Cohens_d, fearImaNeuMinPre_d$Cohens_d,
        fearImaAvNeuMid_d$Cohens_d, fearImaAvMinMid_d$Cohens_d, fearImaNeuMinMid_d$Cohens_d,
        fearImaAvNeuPost_d$Cohens_d, fearImaAvMinPost_d$Cohens_d, fearImaNeuMinPost_d$Cohens_d),
  BF = c(exp(fearImaAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(fearImaAvMinPre_BF@bayesFactor[["bf"]][1]), exp(fearImaNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(fearImaAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(fearImaAvMinMid_BF@bayesFactor[["bf"]][1]), exp(fearImaNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(fearImaAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(fearImaAvMinPost_BF@bayesFactor[["bf"]][1]), exp(fearImaNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableFearIma, file = "supplement/01s1_timeFactor_fear_ima_tTable.doc")



#############################################################################################
### Imagery-based conditioning - unpleasantness ratings - supplemental analyses over time ###
#############################################################################################

# descriptive statistics for unpleasantness ratings in imagery-based conditioning group
describe(dataUnpleas[dataUnpleas$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = unpleasantness rating
anovaUnpleasIma <- ezANOVA(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "ima",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaUnpleasIma$ANOVA$pEtaSq <-
    c(anovaUnpleasIma$ANOVA$SSn[1] / (anovaUnpleasIma$ANOVA$SSd[1]+anovaUnpleasIma$ANOVA$SSn[1]),
      anovaUnpleasIma$ANOVA$SSn[2] / (anovaUnpleasIma$ANOVA$SSd[2]+anovaUnpleasIma$ANOVA$SSn[2]),
      anovaUnpleasIma$ANOVA$SSn[3] / (anovaUnpleasIma$ANOVA$SSd[3]+anovaUnpleasIma$ANOVA$SSn[3]),
      anovaUnpleasIma$ANOVA$SSn[4] / (anovaUnpleasIma$ANOVA$SSd[4]+anovaUnpleasIma$ANOVA$SSn[4])
); print(anovaUnpleasIma)
capture.output(print(anovaUnpleasIma), file = "supplement/01s1_timeFactor_unpleas_ima_anovaFreq.doc")

# bayesian ANOVA on unpleasantness ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFUnpleasIma <- anovaBF(
  formula = unpleasantness ~ CS*time + partInd,
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "ima",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFUnpleasIma)
capture.output(print(anovaBFUnpleasIma), file = "supplement/01s1_timeFactor_unpleas_ima_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for unpleasantness ratings in imagery-based conditioning group
plotUnpleasIma <- ezPlot(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "ima",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotUnpleasIma 
ggsave(plot = plotUnpleasIma, filename = "supplement/01s1_timeFactor_unpleas_ima_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on unpleasantness ratings in imagery-based conditioning group
### Pre
# CS+av vs CS+neu
unpleasImaAvNeuPre_t <- t.test(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                alternative = "two.sided", paired = TRUE) # two-sided
unpleasImaAvNeuPre_d <- cohens_d(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                  paired = TRUE)
unpleasImaAvNeuPre_BF <- ttestBF(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                  nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
unpleasImaAvMinPre_t <- t.test(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                alternative = "two.sided", paired = TRUE) # twp-sided
unpleasImaAvMinPre_d <- cohens_d(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                  paired = TRUE)
unpleasImaAvMinPre_BF <- ttestBF(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                  nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
unpleasImaNeuMinPre_t <- t.test(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                 y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                 alternative = "two.sided", paired = TRUE) # two-sided
unpleasImaNeuMinPre_d <- cohens_d(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                   y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                   paired = TRUE)
unpleasImaNeuMinPre_BF <- ttestBF(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "ima"],
                                   y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "ima"],
                                   nullIntervall = NULL, paired = TRUE) # two-sided

### Mid
# CS+av vs CS+neu
unpleasImaAvNeuMid_t <- t.test(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvNeuMid_d <- cohens_d(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                  paired = TRUE)
unpleasImaAvNeuMid_BF <- ttestBF(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasImaAvMinMid_t <- t.test(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvMinMid_d <- cohens_d(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                  paired = TRUE)
unpleasImaAvMinMid_BF <- ttestBF(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "ima"],
                                  y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasImaNeuMinMid_t <- t.test(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                 y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                 alternative = "two.sided", paired = TRUE) # two-sided
unpleasImaNeuMinMid_d <- cohens_d(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                   y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                   paired = TRUE)
unpleasImaNeuMinMid_BF <- ttestBF(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "ima"],
                                   y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "ima"],
                                   nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
unpleasImaAvNeuPost_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                       y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                       alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvNeuPost_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                              paired = TRUE)
unpleasImaAvNeuPost_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasImaAvMinPost_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                            y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                            alternative = "greater", paired = TRUE) # one-sided
unpleasImaAvMinPost_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                              paired = TRUE)
unpleasImaAvMinPost_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "ima"],
                              y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasImaNeuMinPost_t <- t.test(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                             y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                             alternative = "two.sided", paired = TRUE) # two-sided
unpleasImaNeuMinPost_d <- cohens_d(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               paired = TRUE)
unpleasImaNeuMinPost_BF <- ttestBF(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "ima"],
                               y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "ima"],
                               nullIntervall = NULL, paired = TRUE) # two-sided

tableUnpleasIma <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(unpleasImaAvNeuPre_t$statistic, unpleasImaAvMinPre_t$statistic, unpleasImaNeuMinPre_t$statistic,
        unpleasImaAvNeuMid_t$statistic, unpleasImaAvMinMid_t$statistic, unpleasImaNeuMinMid_t$statistic,
        unpleasImaAvNeuPost_t$statistic, unpleasImaAvMinPost_t$statistic, unpleasImaNeuMinPost_t$statistic),
  df = c(unpleasImaAvNeuPre_t$parameter, unpleasImaAvMinPre_t$parameter, unpleasImaNeuMinPre_t$parameter,
         unpleasImaAvNeuMid_t$parameter, unpleasImaAvMinMid_t$parameter, unpleasImaNeuMinMid_t$parameter,
         unpleasImaAvNeuPost_t$parameter, unpleasImaAvMinPost_t$parameter, unpleasImaNeuMinPost_t$parameter), 
  p = c(unpleasImaAvNeuPre_t$p.value, unpleasImaAvMinPre_t$p.value, unpleasImaNeuMinPre_t$p.value,
        unpleasImaAvNeuMid_t$p.value, unpleasImaAvMinMid_t$p.value, unpleasImaNeuMinMid_t$p.value,
        unpleasImaAvNeuPost_t$p.value, unpleasImaAvMinPost_t$p.value, unpleasImaNeuMinPost_t$p.value),
  d = c(unpleasImaAvNeuPre_d$Cohens_d, unpleasImaAvMinPre_d$Cohens_d, unpleasImaNeuMinPre_d$Cohens_d,
        unpleasImaAvNeuMid_d$Cohens_d, unpleasImaAvMinMid_d$Cohens_d, unpleasImaNeuMinMid_d$Cohens_d,
        unpleasImaAvNeuPost_d$Cohens_d, unpleasImaAvMinPost_d$Cohens_d, unpleasImaNeuMinPost_d$Cohens_d),
  BF = c(exp(unpleasImaAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(unpleasImaAvMinPre_BF@bayesFactor[["bf"]][1]), exp(unpleasImaNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(unpleasImaAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(unpleasImaAvMinMid_BF@bayesFactor[["bf"]][1]), exp(unpleasImaNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(unpleasImaAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(unpleasImaAvMinPost_BF@bayesFactor[["bf"]][1]), exp(unpleasImaNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableUnpleasIma, file = "supplement/01s1_timeFactor_unpleas_ima_tTable.doc")



######################################################################################
### Imagery-based conditioning - arousal ratings - supplemental analyses over time ###
######################################################################################

# descriptive statistics for arousal ratings in imagery-based conditioning group
describe(dataArousal[dataArousal$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = arousal rating
anovaArousalIma <- ezANOVA(
  data = dataArousalLong[dataArousalLong$usGroup == "ima",],
  dv = arousal,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaArousalIma$ANOVA$pEtaSq <-
  c(anovaArousalIma$ANOVA$SSn[1] / (anovaArousalIma$ANOVA$SSd[1]+anovaArousalIma$ANOVA$SSn[1]),
    anovaArousalIma$ANOVA$SSn[2] / (anovaArousalIma$ANOVA$SSd[2]+anovaArousalIma$ANOVA$SSn[2]),
    anovaArousalIma$ANOVA$SSn[3] / (anovaArousalIma$ANOVA$SSd[3]+anovaArousalIma$ANOVA$SSn[3]),
    anovaArousalIma$ANOVA$SSn[4] / (anovaArousalIma$ANOVA$SSd[4]+anovaArousalIma$ANOVA$SSn[4])
  ); print(anovaArousalIma)
capture.output(print(anovaArousalIma), file = "supplement/01s1_timeFactor_arousal_ima_anovaFreq.doc")

# bayesian ANOVA on arousal ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFArousalIma <- anovaBF(
  formula = arousal ~ CS*time + partInd,
  data = dataArousalLong[dataArousalLong$usGroup == "ima",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFArousalIma)
capture.output(print(anovaBFArousalIma), file = "supplement/01s1_timeFactor_arousal_ima_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for arousal ratings in imagery-based conditioning group
plotArousalIma <- ezPlot(
  data = dataArousalLong[dataArousalLong$usGroup == "ima",],
  dv = arousal,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotArousalIma 
ggsave(plot = plotArousalIma, filename = "supplement/01s1_timeFactor_arousal_ima_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on arousal ratings in imagery-based conditioning group
### Pre
# CS+av vs CS+neu
arousalImaAvNeuPre_t <- t.test(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                               y = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                               alternative = "two.sided", paired = TRUE) # two-sided
arousalImaAvNeuPre_d <- cohens_d(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                                 paired = TRUE)
arousalImaAvNeuPre_BF <- ttestBF(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
arousalImaAvMinPre_t <- t.test(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                               y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                               alternative = "two.sided", paired = TRUE) # two-sided
arousalImaAvMinPre_d <- cohens_d(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                                 paired = TRUE)
arousalImaAvMinPre_BF <- ttestBF(x = dataArousal$Av_Pre[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
arousalImaNeuMinPre_t <- t.test(x = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                                y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                                alternative = "two.sided", paired = TRUE) # two-sided
arousalImaNeuMinPre_d <- cohens_d(x = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                                  paired = TRUE)
arousalImaNeuMinPre_BF <- ttestBF(x = dataArousal$Neu_Pre[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Pre[dataArousal$usGroup == "ima"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Mid
# CS+av vs CS+neu
arousalImaAvNeuMid_t <- t.test(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                               y = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                               alternative = "greater", paired = TRUE) # one-sided
arousalImaAvNeuMid_d <- cohens_d(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                                 paired = TRUE)
arousalImaAvNeuMid_BF <- ttestBF(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalImaAvMinMid_t <- t.test(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                               y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                               alternative = "greater", paired = TRUE) # one-sided
arousalImaAvMinMid_d <- cohens_d(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                                 paired = TRUE)
arousalImaAvMinMid_BF <- ttestBF(x = dataArousal$Av_Mid[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
arousalImaNeuMinMid_t <- t.test(x = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                                y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                                alternative = "two.sided", paired = TRUE) # two-sided
arousalImaNeuMinMid_d <- cohens_d(x = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                                  paired = TRUE)
arousalImaNeuMinMid_BF <- ttestBF(x = dataArousal$Neu_Mid[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Mid[dataArousal$usGroup == "ima"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
arousalImaAvNeuPost_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                alternative = "greater", paired = TRUE) # one-sided
arousalImaAvNeuPost_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                  paired = TRUE)
arousalImaAvNeuPost_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalImaAvMinPost_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                alternative = "greater", paired = TRUE) # one-sided
arousalImaAvMinPost_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                  paired = TRUE)
arousalImaAvMinPost_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "ima"],
                                  y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
arousalImaNeuMinPost_t <- t.test(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                 y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                 alternative = "two.sided", paired = TRUE) # two-sided
arousalImaNeuMinPost_d <- cohens_d(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                   y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                   paired = TRUE)
arousalImaNeuMinPost_BF <- ttestBF(x = dataArousal$Neu_Post[dataArousal$usGroup == "ima"],
                                   y = dataArousal$Min_Post[dataArousal$usGroup == "ima"],
                                   nullIntervall = NULL, paired = TRUE) # two-sided

tableArousalIma <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(arousalImaAvNeuPre_t$statistic, arousalImaAvMinPre_t$statistic, arousalImaNeuMinPre_t$statistic,
        arousalImaAvNeuMid_t$statistic, arousalImaAvMinMid_t$statistic, arousalImaNeuMinMid_t$statistic,
        arousalImaAvNeuPost_t$statistic, arousalImaAvMinPost_t$statistic, arousalImaNeuMinPost_t$statistic),
  df = c(arousalImaAvNeuPre_t$parameter, arousalImaAvMinPre_t$parameter, arousalImaNeuMinPre_t$parameter,
         arousalImaAvNeuMid_t$parameter, arousalImaAvMinMid_t$parameter, arousalImaNeuMinMid_t$parameter,
         arousalImaAvNeuPost_t$parameter, arousalImaAvMinPost_t$parameter, arousalImaNeuMinPost_t$parameter), 
  p = c(arousalImaAvNeuPre_t$p.value, arousalImaAvMinPre_t$p.value, arousalImaNeuMinPre_t$p.value,
        arousalImaAvNeuMid_t$p.value, arousalImaAvMinMid_t$p.value, arousalImaNeuMinMid_t$p.value,
        arousalImaAvNeuPost_t$p.value, arousalImaAvMinPost_t$p.value, arousalImaNeuMinPost_t$p.value),
  d = c(arousalImaAvNeuPre_d$Cohens_d, arousalImaAvMinPre_d$Cohens_d, arousalImaNeuMinPre_d$Cohens_d,
        arousalImaAvNeuMid_d$Cohens_d, arousalImaAvMinMid_d$Cohens_d, arousalImaNeuMinMid_d$Cohens_d,
        arousalImaAvNeuPost_d$Cohens_d, arousalImaAvMinPost_d$Cohens_d, arousalImaNeuMinPost_d$Cohens_d),
  BF = c(exp(arousalImaAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(arousalImaAvMinPre_BF@bayesFactor[["bf"]][1]), exp(arousalImaNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(arousalImaAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(arousalImaAvMinMid_BF@bayesFactor[["bf"]][1]), exp(arousalImaNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(arousalImaAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(arousalImaAvMinPost_BF@bayesFactor[["bf"]][1]), exp(arousalImaNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableArousalIma, file = "supplement/01s1_timeFactor_arousal_ima_tTable.doc")






############################################################################
### Classical conditioning - unpleasantness ratings - supplemental analyses over time ###
############################################################################

# descriptive statistics for unpleasantness ratings in classical conditioning group
describe(dataUnpleas[dataUnpleas$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = unpleasantness rating
anovaUnpleasReal <- ezANOVA(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "real",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaUnpleasReal$ANOVA$pEtaSq <-
  c(anovaUnpleasReal$ANOVA$SSn[1] / (anovaUnpleasReal$ANOVA$SSd[1]+anovaUnpleasReal$ANOVA$SSn[1]),
    anovaUnpleasReal$ANOVA$SSn[2] / (anovaUnpleasReal$ANOVA$SSd[2]+anovaUnpleasReal$ANOVA$SSn[2]),
    anovaUnpleasReal$ANOVA$SSn[3] / (anovaUnpleasReal$ANOVA$SSd[3]+anovaUnpleasReal$ANOVA$SSn[3]),
    anovaUnpleasReal$ANOVA$SSn[4] / (anovaUnpleasReal$ANOVA$SSd[4]+anovaUnpleasReal$ANOVA$SSn[4])
  ); print(anovaUnpleasReal)
capture.output(print(anovaUnpleasReal), file = "Supplement/01s_unpleas_real_anovaFreq.doc")

# bayesian ANOVA on unpleasantness ratings in classical conditioning group
set.seed(rngSeed); anovaBFUnpleasReal <- anovaBF(
  formula = unpleasantness ~ CS*time + partInd,
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "real",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFUnpleasReal)
capture.output(print(anovaBFUnpleasReal), file = "Supplement/01s_unpleas_real_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for unpleasantness ratings in classical conditioning group
plotUnpleasReal <- ezPlot(
  data = dataUnpleasLong[dataUnpleasLong$usGroup == "real",],
  dv = unpleasantness,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotUnpleasReal 
ggsave(plot = plotUnpleasReal, filename = "Supplement/01s_unpleas_real_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on unpleasantness ratings in classical conditioning group
### Pre
# CS+av vs CS+neu
unpleasRealAvNeuPre_t <- t.test(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealAvNeuPre_d <- cohens_d(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                                 paired = TRUE)
unpleasRealAvNeuPre_BF <- ttestBF(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
unpleasRealAvMinPre_t <- t.test(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealAvMinPre_d <- cohens_d(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                                 paired = TRUE)
unpleasRealAvMinPre_BF <- ttestBF(x = dataUnpleas$Av_Pre[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
unpleasRealNeuMinPre_t <- t.test(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealNeuMinPre_d <- cohens_d(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                                  paired = TRUE)
unpleasRealNeuMinPre_BF <- ttestBF(x = dataUnpleas$Neu_Pre[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Pre[dataUnpleas$usGroup == "real"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Mid
# CS+av vs CS+neu
unpleasRealAvNeuMid_t <- t.test(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvNeuMid_d <- cohens_d(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                                 paired = TRUE)
unpleasRealAvNeuMid_BF <- ttestBF(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasRealAvMinMid_t <- t.test(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                               y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvMinMid_d <- cohens_d(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                                 paired = TRUE)
unpleasRealAvMinMid_BF <- ttestBF(x = dataUnpleas$Av_Mid[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasRealNeuMinMid_t <- t.test(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealNeuMinMid_d <- cohens_d(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                                  paired = TRUE)
unpleasRealNeuMinMid_BF <- ttestBF(x = dataUnpleas$Neu_Mid[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Mid[dataUnpleas$usGroup == "real"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
unpleasRealAvNeuPost_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvNeuPost_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                  paired = TRUE)
unpleasRealAvNeuPost_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
unpleasRealAvMinPost_t <- t.test(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                alternative = "greater", paired = TRUE) # one-sided
unpleasRealAvMinPost_d <- cohens_d(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                  paired = TRUE)
unpleasRealAvMinPost_BF <- ttestBF(x = dataUnpleas$Av_Post[dataUnpleas$usGroup == "real"],
                                  y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
unpleasRealNeuMinPost_t <- t.test(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                 y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                 alternative = "two.sided", paired = TRUE) # two-sided
unpleasRealNeuMinPost_d <- cohens_d(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                   y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                   paired = TRUE)
unpleasRealNeuMinPost_BF <- ttestBF(x = dataUnpleas$Neu_Post[dataUnpleas$usGroup == "real"],
                                   y = dataUnpleas$Min_Post[dataUnpleas$usGroup == "real"],
                                   nullIntervall = NULL, paired = TRUE) # two-sided

tableUnpleasReal <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(unpleasRealAvNeuPre_t$statistic, unpleasRealAvMinPre_t$statistic, unpleasRealNeuMinPre_t$statistic,
        unpleasRealAvNeuMid_t$statistic, unpleasRealAvMinMid_t$statistic, unpleasRealNeuMinMid_t$statistic,
        unpleasRealAvNeuPost_t$statistic, unpleasRealAvMinPost_t$statistic, unpleasRealNeuMinPost_t$statistic),
  df = c(unpleasRealAvNeuPre_t$parameter, unpleasRealAvMinPre_t$parameter, unpleasRealNeuMinPre_t$parameter,
         unpleasRealAvNeuMid_t$parameter, unpleasRealAvMinMid_t$parameter, unpleasRealNeuMinMid_t$parameter,
         unpleasRealAvNeuPost_t$parameter, unpleasRealAvMinPost_t$parameter, unpleasRealNeuMinPost_t$parameter), 
  p = c(unpleasRealAvNeuPre_t$p.value, unpleasRealAvMinPre_t$p.value, unpleasRealNeuMinPre_t$p.value,
        unpleasRealAvNeuMid_t$p.value, unpleasRealAvMinMid_t$p.value, unpleasRealNeuMinMid_t$p.value,
        unpleasRealAvNeuPost_t$p.value, unpleasRealAvMinPost_t$p.value, unpleasRealNeuMinPost_t$p.value),
  d = c(unpleasRealAvNeuPre_d$Cohens_d, unpleasRealAvMinPre_d$Cohens_d, unpleasRealNeuMinPre_d$Cohens_d,
        unpleasRealAvNeuMid_d$Cohens_d, unpleasRealAvMinMid_d$Cohens_d, unpleasRealNeuMinMid_d$Cohens_d,
        unpleasRealAvNeuPost_d$Cohens_d, unpleasRealAvMinPost_d$Cohens_d, unpleasRealNeuMinPost_d$Cohens_d),
  BF = c(exp(unpleasRealAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(unpleasRealAvMinPre_BF@bayesFactor[["bf"]][1]), exp(unpleasRealNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(unpleasRealAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(unpleasRealAvMinMid_BF@bayesFactor[["bf"]][1]), exp(unpleasRealNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(unpleasRealAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(unpleasRealAvMinPost_BF@bayesFactor[["bf"]][1]), exp(unpleasRealNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableUnpleasReal, file = "Supplement/01s_unpleas_real_tTable.doc")

#####################################################################
### Classical conditioning - arousal ratings - supplemental analyses over time ###
#####################################################################

# descriptive statistics for arousal ratings in classical conditioning group
describe(dataArousal[dataArousal$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = arousal rating
anovaArousalReal <- ezANOVA(
  data = dataArousalLong[dataArousalLong$usGroup == "real",],
  dv = arousal,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaArousalReal$ANOVA$pEtaSq <-
  c(anovaArousalReal$ANOVA$SSn[1] / (anovaArousalReal$ANOVA$SSd[1]+anovaArousalReal$ANOVA$SSn[1]),
    anovaArousalReal$ANOVA$SSn[2] / (anovaArousalReal$ANOVA$SSd[2]+anovaArousalReal$ANOVA$SSn[2]),
    anovaArousalReal$ANOVA$SSn[3] / (anovaArousalReal$ANOVA$SSd[3]+anovaArousalReal$ANOVA$SSn[3]),
    anovaArousalReal$ANOVA$SSn[4] / (anovaArousalReal$ANOVA$SSd[4]+anovaArousalReal$ANOVA$SSn[4])
  ); print(anovaArousalReal)
capture.output(print(anovaArousalReal), file = "Supplement/01s_arousal_real_anovaFreq.doc")

# bayesian ANOVA on arousal ratings in classical conditioning group
set.seed(rngSeed); anovaBFArousalReal <- anovaBF(
  formula = arousal ~ CS*time + partInd,
  data = dataArousalLong[dataArousalLong$usGroup == "real",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFArousalReal)
capture.output(print(anovaBFArousalReal), file = "Supplement/01s_arousal_real_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for arousal ratings in classical conditioning group
plotArousalReal <- ezPlot(
  data = dataArousalLong[dataArousalLong$usGroup == "real",],
  dv = arousal,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotArousalReal 
ggsave(plot = plotArousalReal, filename = "Supplement/01s_arousal_real_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on arousal ratings in classical conditioning group
### Pre
# CS+av vs CS+neu
arousalRealAvNeuPre_t <- t.test(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                               y = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
arousalRealAvNeuPre_d <- cohens_d(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                                 y = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                                 paired = TRUE)
arousalRealAvNeuPre_BF <- ttestBF(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                                 y = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
arousalRealAvMinPre_t <- t.test(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                               y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
arousalRealAvMinPre_d <- cohens_d(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                                 y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                                 paired = TRUE)
arousalRealAvMinPre_BF <- ttestBF(x = dataArousal$Av_Pre[dataArousal$usGroup == "real"],
                                 y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
arousalRealNeuMinPre_t <- t.test(x = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                                y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
arousalRealNeuMinPre_d <- cohens_d(x = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                                  paired = TRUE)
arousalRealNeuMinPre_BF <- ttestBF(x = dataArousal$Neu_Pre[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Pre[dataArousal$usGroup == "real"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Mid
# CS+av vs CS+neu
arousalRealAvNeuMid_t <- t.test(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                               y = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
arousalRealAvNeuMid_d <- cohens_d(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                                 y = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                                 paired = TRUE)
arousalRealAvNeuMid_BF <- ttestBF(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                                 y = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalRealAvMinMid_t <- t.test(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                               y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                               alternative = "greater", paired = TRUE) # one-sided
arousalRealAvMinMid_d <- cohens_d(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                                 y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                                 paired = TRUE)
arousalRealAvMinMid_BF <- ttestBF(x = dataArousal$Av_Mid[dataArousal$usGroup == "real"],
                                 y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                                 nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
arousalRealNeuMinMid_t <- t.test(x = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                                y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                                alternative = "two.sided", paired = TRUE) # two-sided
arousalRealNeuMinMid_d <- cohens_d(x = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                                  paired = TRUE)
arousalRealNeuMinMid_BF <- ttestBF(x = dataArousal$Neu_Mid[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Mid[dataArousal$usGroup == "real"],
                                  nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
arousalRealAvNeuPost_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                alternative = "greater", paired = TRUE) # one-sided
arousalRealAvNeuPost_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                  y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                  paired = TRUE)
arousalRealAvNeuPost_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                  y = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
arousalRealAvMinPost_t <- t.test(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                alternative = "greater", paired = TRUE) # one-sided
arousalRealAvMinPost_d <- cohens_d(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                  paired = TRUE)
arousalRealAvMinPost_BF <- ttestBF(x = dataArousal$Av_Post[dataArousal$usGroup == "real"],
                                  y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                  nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
arousalRealNeuMinPost_t <- t.test(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                 y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                 alternative = "two.sided", paired = TRUE) # two-sided
arousalRealNeuMinPost_d <- cohens_d(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                   y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                   paired = TRUE)
arousalRealNeuMinPost_BF <- ttestBF(x = dataArousal$Neu_Post[dataArousal$usGroup == "real"],
                                   y = dataArousal$Min_Post[dataArousal$usGroup == "real"],
                                   nullIntervall = NULL, paired = TRUE) # two-sided

tableArousalReal <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(arousalRealAvNeuPre_t$statistic, arousalRealAvMinPre_t$statistic, arousalRealNeuMinPre_t$statistic,
        arousalRealAvNeuMid_t$statistic, arousalRealAvMinMid_t$statistic, arousalRealNeuMinMid_t$statistic,
        arousalRealAvNeuPost_t$statistic, arousalRealAvMinPost_t$statistic, arousalRealNeuMinPost_t$statistic),
  df = c(arousalRealAvNeuPre_t$parameter, arousalRealAvMinPre_t$parameter, arousalRealNeuMinPre_t$parameter,
         arousalRealAvNeuMid_t$parameter, arousalRealAvMinMid_t$parameter, arousalRealNeuMinMid_t$parameter,
         arousalRealAvNeuPost_t$parameter, arousalRealAvMinPost_t$parameter, arousalRealNeuMinPost_t$parameter), 
  p = c(arousalRealAvNeuPre_t$p.value, arousalRealAvMinPre_t$p.value, arousalRealNeuMinPre_t$p.value,
        arousalRealAvNeuMid_t$p.value, arousalRealAvMinMid_t$p.value, arousalRealNeuMinMid_t$p.value,
        arousalRealAvNeuPost_t$p.value, arousalRealAvMinPost_t$p.value, arousalRealNeuMinPost_t$p.value),
  d = c(arousalRealAvNeuPre_d$Cohens_d, arousalRealAvMinPre_d$Cohens_d, arousalRealNeuMinPre_d$Cohens_d,
        arousalRealAvNeuMid_d$Cohens_d, arousalRealAvMinMid_d$Cohens_d, arousalRealNeuMinMid_d$Cohens_d,
        arousalRealAvNeuPost_d$Cohens_d, arousalRealAvMinPost_d$Cohens_d, arousalRealNeuMinPost_d$Cohens_d),
  BF = c(exp(arousalRealAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(arousalRealAvMinPre_BF@bayesFactor[["bf"]][1]), exp(arousalRealNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(arousalRealAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(arousalRealAvMinMid_BF@bayesFactor[["bf"]][1]), exp(arousalRealNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(arousalRealAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(arousalRealAvMinPost_BF@bayesFactor[["bf"]][1]), exp(arousalRealNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableArousalReal, file = "Supplement/01s_arousal_real_tTable.doc")



##################################################################
### Classical conditioning - fear ratings - supplemental analyses over time ###
##################################################################

# descriptive statistics for fear ratings in classical conditioning group
describe(dataFear[dataFear$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = fear rating
anovaFearReal <- ezANOVA(
  data = dataFearLong[dataFearLong$usGroup == "real",],
  dv = fear,
  wid = partInd,
  within = .(CS, time),
  type = 3,
  detailed = TRUE
); anovaFearReal$ANOVA$pEtaSq <-
  c(anovaFearReal$ANOVA$SSn[1] / (anovaFearReal$ANOVA$SSd[1]+anovaFearReal$ANOVA$SSn[1]),
    anovaFearReal$ANOVA$SSn[2] / (anovaFearReal$ANOVA$SSd[2]+anovaFearReal$ANOVA$SSn[2]),
    anovaFearReal$ANOVA$SSn[3] / (anovaFearReal$ANOVA$SSd[3]+anovaFearReal$ANOVA$SSn[3]),
    anovaFearReal$ANOVA$SSn[4] / (anovaFearReal$ANOVA$SSd[4]+anovaFearReal$ANOVA$SSn[4])
  ); print(anovaFearReal)
capture.output(print(anovaFearReal), file = "Supplement/01s_fear_real_anovaFreq.doc")

# bayesian ANOVA on fear ratings in classical conditioning group
set.seed(rngSeed); anovaBFFearReal <- anovaBF(
  formula = fear ~ CS*time + partInd,
  data = dataFearLong[dataFearLong$usGroup == "real",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFearReal)
capture.output(print(anovaBFFearReal), file = "Supplement/01s_fear_real_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for fear ratings in classical conditioning group
plotFearReal <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "real",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearReal 
ggsave(plot = plotFearReal, filename = "Supplement/01s_fear_real_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on fear ratings in classical conditioning group
### Pre
# CS+av vs CS+neu
fearRealAvNeuPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                            y = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                            alternative = "two.sided", paired = TRUE) # two-sided
fearRealAvNeuPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                              paired = TRUE)
fearRealAvNeuPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
fearRealAvMinPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                            y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                            alternative = "two.sided", paired = TRUE) # two-sided
fearRealAvMinPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                              y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                              paired = TRUE)
fearRealAvMinPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "real"],
                              y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
fearRealNeuMinPre_t <- t.test(x = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                             y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                             alternative = "two.sided", paired = TRUE) # two-sided
fearRealNeuMinPre_d <- cohens_d(x = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                               paired = TRUE)
fearRealNeuMinPre_BF <- ttestBF(x = dataFear$Neu_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Min_Pre[dataFear$usGroup == "real"],
                               nullIntervall = NULL, paired = TRUE) # two-sided
### Mid
# CS+av vs CS+neu
fearRealAvNeuMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                            y = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                            alternative = "greater", paired = TRUE) # one-sided
fearRealAvNeuMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                              paired = TRUE)
fearRealAvNeuMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearRealAvMinMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                            y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                            alternative = "greater", paired = TRUE) # one-sided
fearRealAvMinMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                              y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                              paired = TRUE)
fearRealAvMinMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "real"],
                              y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearRealNeuMinMid_t <- t.test(x = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                             y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                             alternative = "two.sided", paired = TRUE) # two-sided
fearRealNeuMinMid_d <- cohens_d(x = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                               paired = TRUE)
fearRealNeuMinMid_BF <- ttestBF(x = dataFear$Neu_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Min_Mid[dataFear$usGroup == "real"],
                               nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
fearRealAvNeuPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                             y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
fearRealAvNeuPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                               y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                               paired = TRUE)
fearRealAvNeuPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                               y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearRealAvMinPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                             y = dataFear$Min_Post[dataFear$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
fearRealAvMinPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                               y = dataFear$Min_Post[dataFear$usGroup == "real"],
                               paired = TRUE)
fearRealAvMinPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                               y = dataFear$Min_Post[dataFear$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearRealNeuMinPost_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                              y = dataFear$Min_Post[dataFear$usGroup == "real"],
                              alternative = "two.sided", paired = TRUE) # two-sided
fearRealNeuMinPost_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                                y = dataFear$Min_Post[dataFear$usGroup == "real"],
                                paired = TRUE)
fearRealNeuMinPost_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                                y = dataFear$Min_Post[dataFear$usGroup == "real"],
                                nullIntervall = NULL, paired = TRUE) # two-sided

tableFearReal <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(fearRealAvNeuPre_t$statistic, fearRealAvMinPre_t$statistic, fearRealNeuMinPre_t$statistic,
        fearRealAvNeuMid_t$statistic, fearRealAvMinMid_t$statistic, fearRealNeuMinMid_t$statistic,
        fearRealAvNeuPost_t$statistic, fearRealAvMinPost_t$statistic, fearRealNeuMinPost_t$statistic),
  df = c(fearRealAvNeuPre_t$parameter, fearRealAvMinPre_t$parameter, fearRealNeuMinPre_t$parameter,
         fearRealAvNeuMid_t$parameter, fearRealAvMinMid_t$parameter, fearRealNeuMinMid_t$parameter,
         fearRealAvNeuPost_t$parameter, fearRealAvMinPost_t$parameter, fearRealNeuMinPost_t$parameter), 
  p = c(fearRealAvNeuPre_t$p.value, fearRealAvMinPre_t$p.value, fearRealNeuMinPre_t$p.value,
        fearRealAvNeuMid_t$p.value, fearRealAvMinMid_t$p.value, fearRealNeuMinMid_t$p.value,
        fearRealAvNeuPost_t$p.value, fearRealAvMinPost_t$p.value, fearRealNeuMinPost_t$p.value),
  d = c(fearRealAvNeuPre_d$Cohens_d, fearRealAvMinPre_d$Cohens_d, fearRealNeuMinPre_d$Cohens_d,
        fearRealAvNeuMid_d$Cohens_d, fearRealAvMinMid_d$Cohens_d, fearRealNeuMinMid_d$Cohens_d,
        fearRealAvNeuPost_d$Cohens_d, fearRealAvMinPost_d$Cohens_d, fearRealNeuMinPost_d$Cohens_d),
  BF = c(exp(fearRealAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(fearRealAvMinPre_BF@bayesFactor[["bf"]][1]), exp(fearRealNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(fearRealAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(fearRealAvMinMid_BF@bayesFactor[["bf"]][1]), exp(fearRealNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(fearRealAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(fearRealAvMinPost_BF@bayesFactor[["bf"]][1]), exp(fearRealNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableFearReal, file = "Supplement/01s_fear_real_tTable.doc")



###################################################################
### Across groups - unpleasantness ratings - supplemental analyses over time ###
###################################################################

# descriptive statistics  for unpleasantess ratings across conditioning groups
describe(dataUnpleas)

# frequentist ANOVA on unpleasantness ratings across conditioning groups
anovaUnpleas <- ezANOVA(
  data = dataUnpleasLong,
  dv = unpleasantness,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaUnpleas$ANOVA$pEtaSq <- 
  c(anovaUnpleas$ANOVA$SSn[1] / (anovaUnpleas$ANOVA$SSd[1]+anovaUnpleas$ANOVA$SSn[1]),
    anovaUnpleas$ANOVA$SSn[2] / (anovaUnpleas$ANOVA$SSd[2]+anovaUnpleas$ANOVA$SSn[2]),
    anovaUnpleas$ANOVA$SSn[3] / (anovaUnpleas$ANOVA$SSd[3]+anovaUnpleas$ANOVA$SSn[3]),
    anovaUnpleas$ANOVA$SSn[4] / (anovaUnpleas$ANOVA$SSd[4]+anovaUnpleas$ANOVA$SSn[4]),
    anovaUnpleas$ANOVA$SSn[5] / (anovaUnpleas$ANOVA$SSd[5]+anovaUnpleas$ANOVA$SSn[5]),
    anovaUnpleas$ANOVA$SSn[6] / (anovaUnpleas$ANOVA$SSd[6]+anovaUnpleas$ANOVA$SSn[6]),
    anovaUnpleas$ANOVA$SSn[7] / (anovaUnpleas$ANOVA$SSd[7]+anovaUnpleas$ANOVA$SSn[7]),
    anovaUnpleas$ANOVA$SSn[8] / (anovaUnpleas$ANOVA$SSd[8]+anovaUnpleas$ANOVA$SSn[8])
); print(anovaUnpleas)
capture.output(print(anovaUnpleas), file = "Supplement/01s_unpleas_both_anovaFreq.doc")

# bayesian ANOVA on unpleasantness ratings across conditioning groups
set.seed(rngSeed); anovaBFUnpleas <- anovaBF(
  formula = unpleasantness ~ usGroup*CS*time + partInd,
  data = dataUnpleasLong,
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFUnpleas)
capture.output(print(anovaBFUnpleas), file = "Supplement/01s_unpleas_both_anovaBayes.doc")

# quick & dirty graph of group x CS ANOVA on valence ratings
plotUnpleas <- ezPlot(
  data = dataUnpleasLong,
  dv = unpleasantness,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,  
  col = usGroup
)  
ggsave(plot = plotUnpleas, filename = "Supplement/01s_unpleas_both_plot.jpg",
       width = 20, height = 10, units = "cm")



############################################################
### Across groups - arousal ratings - supplemental analyses over time ###
############################################################

# descriptive statistics for arousal ratings across conditioning groups
describe(dataArousal)

# frequentist ANOVA on arousal ratings across conditioning groups
anovaArousal <- ezANOVA(
  data = dataArousalLong,
  dv = arousal,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaArousal$ANOVA$pEtaSq <- 
  c(anovaArousal$ANOVA$SSn[1] / (anovaArousal$ANOVA$SSd[1]+anovaArousal$ANOVA$SSn[1]),
    anovaArousal$ANOVA$SSn[2] / (anovaArousal$ANOVA$SSd[2]+anovaArousal$ANOVA$SSn[2]),
    anovaArousal$ANOVA$SSn[3] / (anovaArousal$ANOVA$SSd[3]+anovaArousal$ANOVA$SSn[3]),
    anovaArousal$ANOVA$SSn[4] / (anovaArousal$ANOVA$SSd[4]+anovaArousal$ANOVA$SSn[4]),
    anovaArousal$ANOVA$SSn[5] / (anovaArousal$ANOVA$SSd[5]+anovaArousal$ANOVA$SSn[5]),
    anovaArousal$ANOVA$SSn[6] / (anovaArousal$ANOVA$SSd[6]+anovaArousal$ANOVA$SSn[6]),
    anovaArousal$ANOVA$SSn[7] / (anovaArousal$ANOVA$SSd[7]+anovaArousal$ANOVA$SSn[7]),
    anovaArousal$ANOVA$SSn[8] / (anovaArousal$ANOVA$SSd[8]+anovaArousal$ANOVA$SSn[8])
  ); print(anovaArousal)
capture.output(print(anovaArousal), file = "Supplement/01s_arousal_both_anovaFreq.doc")

# bayesian ANOVA on arousal ratings across conditioning groups
set.seed(rngSeed); anovaBFArousal <- anovaBF(
  formula = arousal ~ usGroup*CS*time + partInd,
  data = dataArousalLong,
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFArousal)
capture.output(print(anovaBFArousal), file = "Supplement/01s_arousal_both_anovaBayes.doc")

# quick & dirty graph of group x CS ANOVA on valence ratings
plotArousal <- ezPlot(
  data = dataArousalLong,
  dv = arousal,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,  
  col = usGroup
)  
ggsave(plot = plotArousal, filename = "Supplement/01s_arousal_both_plot.jpg",
       width = 20, height = 10, units = "cm")



#########################################################
### Across groups - fear ratings - supplemental analyses over time ###
#########################################################

# descriptive statistics for fear ratings across conditioning groups
describe(dataFear)

# frequentist ANOVA on fear ratings across conditioning groups
anovaFear <- ezANOVA(
  data = dataFearLong,
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaFear$ANOVA$pEtaSq <- 
  c(anovaFear$ANOVA$SSn[1] / (anovaFear$ANOVA$SSd[1]+anovaFear$ANOVA$SSn[1]),
    anovaFear$ANOVA$SSn[2] / (anovaFear$ANOVA$SSd[2]+anovaFear$ANOVA$SSn[2]),
    anovaFear$ANOVA$SSn[3] / (anovaFear$ANOVA$SSd[3]+anovaFear$ANOVA$SSn[3]),
    anovaFear$ANOVA$SSn[4] / (anovaFear$ANOVA$SSd[4]+anovaFear$ANOVA$SSn[4]),
    anovaFear$ANOVA$SSn[5] / (anovaFear$ANOVA$SSd[5]+anovaFear$ANOVA$SSn[5]),
    anovaFear$ANOVA$SSn[6] / (anovaFear$ANOVA$SSd[6]+anovaFear$ANOVA$SSn[6]),
    anovaFear$ANOVA$SSn[7] / (anovaFear$ANOVA$SSd[7]+anovaFear$ANOVA$SSn[7]),
    anovaFear$ANOVA$SSn[8] / (anovaFear$ANOVA$SSd[8]+anovaFear$ANOVA$SSn[8])
  ); print(anovaFear)
capture.output(print(anovaFear), file = "Supplement/01s_fear_both_anovaFreq.doc")

# bayesian ANOVA on fear ratings across conditioning groups
set.seed(rngSeed); anovaBFFear <- anovaBF(
  formula = fear ~ usGroup*CS*time + partInd,
  data = dataFearLong,
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFear)
capture.output(print(anovaBFFear), file = "Supplement/01s_fear_both_anovaBayes.doc")

# quick & dirty graph of group x CS ANOVA on valence ratings
plotFear <- ezPlot(
  data = dataFearLong,
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,  
  col = usGroup
)  
ggsave(plot = plotFear, filename = "Supplement/01s_fear_both_plot.jpg",
       width = 20, height = 10, units = "cm")