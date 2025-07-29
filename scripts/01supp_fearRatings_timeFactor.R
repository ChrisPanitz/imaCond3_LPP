# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: Supplementary analyses on fear ratings: adding factor "Time"

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




##############################################################################
### Across groups - fear ratings - supplementary analyses with time factor ###
##############################################################################

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
capture.output(print(anovaFear), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_acrossGroups_anovaFreq.doc"))

# bayesian ANOVA on fear ratings across conditioning groups
set.seed(rngSeed); anovaBFFear <- generalTestBF(
  formula = fear ~ usGroup*CS*time + partInd + partInd:CS + partInd:time,
  data = dataFearLong,
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 10000 # only 10,000 iterations because it has to compute 128 models
) 

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFFear@bayesFactor$bf <- log(exp(anovaBFFear@bayesFactor$bf) / 
                                    exp(anovaBFFear@bayesFactor$bf[length(anovaBFFear@bayesFactor$bf)]))
anovaBFFear@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFFear)
capture.output(print(anovaBFFear), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_acrossGroups_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFFear)
capture.output(bf_inclusion(anovaBFFear), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_acrossGroups_BFinclusion.doc"))

# quick graphs of group x CS ANOVA on fear ratings
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
ggsave(plot = plotFear, filename = paste0(pathname, "/supplement/01s_fearRating_timeFactor_acrossGroups_plot.jpg"),
       width = 20, height = 10, units = "cm")

# frequentist & bayesian t-tests on fear ratings (difference scores) across groups
# Pre
# delta [CS+av - CS+neu]
fearBetweenAvNeuPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                               dataFear$Neu_Pre[dataFear$usGroup == "real"],
                             y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                               dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvNeuPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvNeuPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Pre[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
fearBetweenAvMinPre_t <- t.test(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                               dataFear$Min_Pre[dataFear$usGroup == "real"],
                             y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                               dataFear$Min_Pre[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvMinPre_d <- cohens_d(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                                 dataFear$Min_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Pre[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvMinPre_BF <- ttestBF(x = dataFear$Av_Pre[dataFear$usGroup == "real"] -
                                 dataFear$Min_Pre[dataFear$usGroup == "real"],
                               y = dataFear$Av_Pre[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Pre[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
fearBetweenNeuMinPre_t <- t.test(x = dataFear$Neu_Pre[dataFear$usGroup == "real"] -
                                dataFear$Min_Pre[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Pre[dataFear$usGroup == "ima"] - 
                                dataFear$Min_Pre[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenNeuMinPre_d <- cohens_d(x = dataFear$Neu_Pre[dataFear$usGroup == "real"] -
                                  dataFear$Min_Pre[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Pre[dataFear$usGroup == "ima"] -
                                  dataFear$Min_Pre[dataFear$usGroup == "ima"],
                                paired = FALSE)
fearBetweenNeuMinPre_BF <- ttestBF(x = dataFear$Neu_Pre[dataFear$usGroup == "real"] -
                                  dataFear$Min_Pre[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Pre[dataFear$usGroup == "ima"] - 
                                  dataFear$Min_Pre[dataFear$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
### Mid
# delta [CS+av - CS+neu]
fearBetweenAvNeuMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                               dataFear$Neu_Mid[dataFear$usGroup == "real"],
                             y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                               dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvNeuMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvNeuMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Mid[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
fearBetweenAvMinMid_t <- t.test(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                               dataFear$Min_Mid[dataFear$usGroup == "real"],
                             y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                               dataFear$Min_Mid[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvMinMid_d <- cohens_d(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                                 dataFear$Min_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Mid[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvMinMid_BF <- ttestBF(x = dataFear$Av_Mid[dataFear$usGroup == "real"] -
                                 dataFear$Min_Mid[dataFear$usGroup == "real"],
                               y = dataFear$Av_Mid[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Mid[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
fearBetweenNeuMinMid_t <- t.test(x = dataFear$Neu_Mid[dataFear$usGroup == "real"] -
                                dataFear$Min_Mid[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Mid[dataFear$usGroup == "ima"] - 
                                dataFear$Min_Mid[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenNeuMinMid_d <- cohens_d(x = dataFear$Neu_Mid[dataFear$usGroup == "real"] -
                                  dataFear$Min_Mid[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Mid[dataFear$usGroup == "ima"] -
                                  dataFear$Min_Mid[dataFear$usGroup == "ima"],
                                paired = FALSE)
fearBetweenNeuMinMid_BF <- ttestBF(x = dataFear$Neu_Mid[dataFear$usGroup == "real"] -
                                  dataFear$Min_Mid[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Mid[dataFear$usGroup == "ima"] - 
                                  dataFear$Min_Mid[dataFear$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
### Post
# delta [CS+av - CS+neu]
fearBetweenAvNeuPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                dataFear$Neu_Post[dataFear$usGroup == "real"],
                              y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                dataFear$Neu_Post[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvNeuPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                  dataFear$Neu_Post[dataFear$usGroup == "real"],
                                y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                  dataFear$Neu_Post[dataFear$usGroup == "ima"],
                                paired = FALSE)
fearBetweenAvNeuPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                  dataFear$Neu_Post[dataFear$usGroup == "real"],
                                y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                  dataFear$Neu_Post[dataFear$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
fearBetweenAvMinPost_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                dataFear$Min_Post[dataFear$usGroup == "real"],
                              y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                dataFear$Min_Post[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvMinPost_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                  dataFear$Min_Post[dataFear$usGroup == "real"],
                                y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                  dataFear$Min_Post[dataFear$usGroup == "ima"],
                                paired = FALSE)
fearBetweenAvMinPost_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                  dataFear$Min_Post[dataFear$usGroup == "real"],
                                y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                  dataFear$Min_Post[dataFear$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
fearBetweenNeuMinPost_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                 dataFear$Min_Post[dataFear$usGroup == "real"],
                               y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                                 dataFear$Min_Post[dataFear$usGroup == "ima"],
                               alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenNeuMinPost_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                   dataFear$Min_Post[dataFear$usGroup == "real"],
                                 y = dataFear$Neu_Post[dataFear$usGroup == "ima"] -
                                   dataFear$Min_Post[dataFear$usGroup == "ima"],
                                 paired = FALSE)
fearBetweenNeuMinPost_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                   dataFear$Min_Post[dataFear$usGroup == "real"],
                                 y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                                   dataFear$Min_Post[dataFear$usGroup == "ima"],
                                 nullInterval = NULL, paired = FALSE) # two-sided

tableFearBetween <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(fearBetweenAvNeuPre_t$statistic, fearBetweenAvMinPre_t$statistic, fearBetweenNeuMinPre_t$statistic,
        fearBetweenAvNeuMid_t$statistic, fearBetweenAvMinMid_t$statistic, fearBetweenNeuMinMid_t$statistic,
        fearBetweenAvNeuPost_t$statistic, fearBetweenAvMinPost_t$statistic, fearBetweenNeuMinPost_t$statistic),
  df = c(fearBetweenAvNeuPre_t$parameter, fearBetweenAvMinPre_t$parameter, fearBetweenNeuMinPre_t$parameter,
         fearBetweenAvNeuMid_t$parameter, fearBetweenAvMinMid_t$parameter, fearBetweenNeuMinMid_t$parameter,
         fearBetweenAvNeuPost_t$parameter, fearBetweenAvMinPost_t$parameter, fearBetweenNeuMinPost_t$parameter), 
  p = c(fearBetweenAvNeuPre_t$p.value*3, fearBetweenAvMinPre_t$p.value*3, fearBetweenNeuMinPre_t$p.value*3,
        fearBetweenAvNeuMid_t$p.value*3, fearBetweenAvMinMid_t$p.value*3, fearBetweenNeuMinMid_t$p.value*3,
        fearBetweenAvNeuPost_t$p.value*3, fearBetweenAvMinPost_t$p.value*3, fearBetweenNeuMinPost_t$p.value*3),
  d = c(fearBetweenAvNeuPre_d$Cohens_d, fearBetweenAvMinPre_d$Cohens_d, fearBetweenNeuMinPre_d$Cohens_d,
        fearBetweenAvNeuMid_d$Cohens_d, fearBetweenAvMinMid_d$Cohens_d, fearBetweenNeuMinMid_d$Cohens_d,
        fearBetweenAvNeuPost_d$Cohens_d, fearBetweenAvMinPost_d$Cohens_d, fearBetweenNeuMinPost_d$Cohens_d),
  BF = c(exp(fearBetweenAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(fearBetweenAvMinPre_BF@bayesFactor[["bf"]][1]), exp(fearBetweenNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(fearBetweenAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(fearBetweenAvMinMid_BF@bayesFactor[["bf"]][1]), exp(fearBetweenNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(fearBetweenAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(fearBetweenAvMinPost_BF@bayesFactor[["bf"]][1]), exp(fearBetweenNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = rep("two.sided",9)
)
tableFearBetween$p[tableFearBetween$p > 1] <- 1
capture.output(tableFearBetween, file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_betweenGroups_tTable.doc"))

# frequentist & bayesian t-tests on fear ratings across conditioning groups
### Pre
# CS+av vs CS+neu
fearAcrossAvNeuPre_t <- t.test(x = dataFear$Av_Pre,
                            y = dataFear$Neu_Pre,
                            alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossAvNeuPre_d <- cohens_d(x = dataFear$Av_Pre,
                              y = dataFear$Neu_Pre,
                              paired = TRUE)
fearAcrossAvNeuPre_BF <- ttestBF(x = dataFear$Av_Pre,
                              y = dataFear$Neu_Pre,
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+av vs CS-
fearAcrossAvMinPre_t <- t.test(x = dataFear$Av_Pre,
                            y = dataFear$Min_Pre,
                            alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossAvMinPre_d <- cohens_d(x = dataFear$Av_Pre,
                              y = dataFear$Min_Pre,
                              paired = TRUE)
fearAcrossAvMinPre_BF <- ttestBF(x = dataFear$Av_Pre,
                              y = dataFear$Min_Pre,
                              nullInterval = NULL, paired = TRUE) # two-sided
# CS+neu vs CS-
fearAcrossNeuMinPre_t <- t.test(x = dataFear$Neu_Pre,
                             y = dataFear$Min_Pre,
                             alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossNeuMinPre_d <- cohens_d(x = dataFear$Neu_Pre,
                               y = dataFear$Min_Pre,
                               paired = TRUE)
fearAcrossNeuMinPre_BF <- ttestBF(x = dataFear$Neu_Pre,
                               y = dataFear$Min_Pre,
                               nullIntervall = NULL, paired = TRUE) # two-sided
### Mid
# CS+av vs CS+neu
fearAcrossAvNeuMid_t <- t.test(x = dataFear$Av_Mid,
                            y = dataFear$Neu_Mid,
                            alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvNeuMid_d <- cohens_d(x = dataFear$Av_Mid,
                              y = dataFear$Neu_Mid,
                              paired = TRUE)
fearAcrossAvNeuMid_BF <- ttestBF(x = dataFear$Av_Mid,
                              y = dataFear$Neu_Mid,
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearAcrossAvMinMid_t <- t.test(x = dataFear$Av_Mid,
                            y = dataFear$Min_Mid,
                            alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvMinMid_d <- cohens_d(x = dataFear$Av_Mid,
                              y = dataFear$Min_Mid,
                              paired = TRUE)
fearAcrossAvMinMid_BF <- ttestBF(x = dataFear$Av_Mid,
                              y = dataFear$Min_Mid,
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearAcrossNeuMinMid_t <- t.test(x = dataFear$Neu_Mid,
                             y = dataFear$Min_Mid,
                             alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossNeuMinMid_d <- cohens_d(x = dataFear$Neu_Mid,
                               y = dataFear$Min_Mid,
                               paired = TRUE)
fearAcrossNeuMinMid_BF <- ttestBF(x = dataFear$Neu_Mid,
                               y = dataFear$Min_Mid,
                               nullIntervall = NULL, paired = TRUE) # two-sided

### Post
# CS+av vs CS+neu
fearAcrossAvNeuPost_t <- t.test(x = dataFear$Av_Post,
                             y = dataFear$Neu_Post,
                             alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvNeuPost_d <- cohens_d(x = dataFear$Av_Post,
                               y = dataFear$Neu_Post,
                               paired = TRUE)
fearAcrossAvNeuPost_BF <- ttestBF(x = dataFear$Av_Post,
                               y = dataFear$Neu_Post,
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearAcrossAvMinPost_t <- t.test(x = dataFear$Av_Post,
                             y = dataFear$Min_Post,
                             alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvMinPost_d <- cohens_d(x = dataFear$Av_Post,
                               y = dataFear$Min_Post,
                               paired = TRUE)
fearAcrossAvMinPost_BF <- ttestBF(x = dataFear$Av_Post,
                               y = dataFear$Min_Post,
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearAcrossNeuMinPost_t <- t.test(x = dataFear$Neu_Post,
                              y = dataFear$Min_Post,
                              alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossNeuMinPost_d <- cohens_d(x = dataFear$Neu_Post,
                                y = dataFear$Min_Post,
                                paired = TRUE)
fearAcrossNeuMinPost_BF <- ttestBF(x = dataFear$Neu_Post,
                                y = dataFear$Min_Post,
                                nullIntervall = NULL, paired = TRUE) # two-sided

tableFearAcross <- data.frame(
  time = c(rep("Pre",3), rep("Mid",3), rep("Post",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 3),
  t = c(fearAcrossAvNeuPre_t$statistic, fearAcrossAvMinPre_t$statistic, fearAcrossNeuMinPre_t$statistic,
        fearAcrossAvNeuMid_t$statistic, fearAcrossAvMinMid_t$statistic, fearAcrossNeuMinMid_t$statistic,
        fearAcrossAvNeuPost_t$statistic, fearAcrossAvMinPost_t$statistic, fearAcrossNeuMinPost_t$statistic),
  df = c(fearAcrossAvNeuPre_t$parameter, fearAcrossAvMinPre_t$parameter, fearAcrossNeuMinPre_t$parameter,
         fearAcrossAvNeuMid_t$parameter, fearAcrossAvMinMid_t$parameter, fearAcrossNeuMinMid_t$parameter,
         fearAcrossAvNeuPost_t$parameter, fearAcrossAvMinPost_t$parameter, fearAcrossNeuMinPost_t$parameter), 
  p = c(fearAcrossAvNeuPre_t$p.value, fearAcrossAvMinPre_t$p.value, fearAcrossNeuMinPre_t$p.value,
        fearAcrossAvNeuMid_t$p.value, fearAcrossAvMinMid_t$p.value, fearAcrossNeuMinMid_t$p.value,
        fearAcrossAvNeuPost_t$p.value, fearAcrossAvMinPost_t$p.value, fearAcrossNeuMinPost_t$p.value),
  d = c(fearAcrossAvNeuPre_d$Cohens_d, fearAcrossAvMinPre_d$Cohens_d, fearAcrossNeuMinPre_d$Cohens_d,
        fearAcrossAvNeuMid_d$Cohens_d, fearAcrossAvMinMid_d$Cohens_d, fearAcrossNeuMinMid_d$Cohens_d,
        fearAcrossAvNeuPost_d$Cohens_d, fearAcrossAvMinPost_d$Cohens_d, fearAcrossNeuMinPost_d$Cohens_d),
  BF = c(exp(fearAcrossAvNeuPre_BF@bayesFactor[["bf"]][1]), exp(fearAcrossAvMinPre_BF@bayesFactor[["bf"]][1]), exp(fearAcrossNeuMinPre_BF@bayesFactor[["bf"]][1]),
         exp(fearAcrossAvNeuMid_BF@bayesFactor[["bf"]][1]), exp(fearAcrossAvMinMid_BF@bayesFactor[["bf"]][1]), exp(fearAcrossNeuMinMid_BF@bayesFactor[["bf"]][1]),
         exp(fearAcrossAvNeuPost_BF@bayesFactor[["bf"]][1]), exp(fearAcrossAvMinPost_BF@bayesFactor[["bf"]][1]), exp(fearAcrossNeuMinPost_BF@bayesFactor[["bf"]][1])),
  testDir = c("two.sided","two.sided","two.sided", rep(c("one.sided","one.sided","two.sided"),2))
)
capture.output(tableFearAcross, file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_acrossGroups_tTable.doc"))



###########################################################################################
### Imagery-based conditioning - fear ratings - supplementary analyses with time factor ###
###########################################################################################

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
capture.output(print(anovaFearIma), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_ima_anovaFreq.doc"))

# bayesian ANOVA on fear ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFFearIma <- generalTestBF(
  formula = fear ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
) 

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFFearIma@bayesFactor$bf <- log(exp(anovaBFFearIma@bayesFactor$bf) / 
                                       exp(anovaBFFearIma@bayesFactor$bf[length(anovaBFFearIma@bayesFactor$bf)]))
anovaBFFearIma@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFFearIma)
capture.output(print(anovaBFFearIma), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_ima_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFFearIma)
capture.output(bf_inclusion(anovaBFFearIma), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_ima_BFinclusion.doc"))

# quick graph of CS Type x time ANOVA for fear ratings in imagery-based conditioning group
plotFearIma <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearIma 
ggsave(plot = plotFearIma, filename = paste0(pathname, "/supplement/01s_fearRating_timeFactor_ima_plot.jpg"),
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
capture.output(tableFearIma, file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_ima_tTable.doc"))



#######################################################################################
### Classical conditioning - fear ratings - supplementary analyses with time factor ###
#######################################################################################

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
capture.output(print(anovaFearReal), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_real_anovaFreq.doc"))

# bayesian ANOVA on fear ratings in classical conditioning group
set.seed(rngSeed); anovaBFFearReal <- generalTestBF(
  formula = fear ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataFearLong[dataFearLong$usGroup == "real",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
) 

# compute Bayes factor relative to null model including random slopes instead
# of intercept-only null model
anovaBFFearReal@bayesFactor$bf <- log(exp(anovaBFFearReal@bayesFactor$bf) / 
                                      exp(anovaBFFearReal@bayesFactor$bf[length(anovaBFFearReal@bayesFactor$bf)]))
anovaBFFearReal@denominator@longName <- "Intercept and random slopes only"

# show and save results
print(anovaBFFearReal)
capture.output(print(anovaBFFearReal), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_real_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFFearReal)
capture.output(bf_inclusion(anovaBFFearReal), file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_real_BFinclusion.doc"))

# quick graph of CS Type x time ANOVA for fear ratings in classical conditioning group
plotFearReal <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "real",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearReal 
ggsave(plot = plotFearReal, filename = paste0(pathname, "/supplement/01s_fearRating_timeFactor_real_plot.jpg"),
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
capture.output(tableFearReal, file = paste0(pathname, "/supplement/01s_fearRating_timeFactor_real_tTable.doc"))

