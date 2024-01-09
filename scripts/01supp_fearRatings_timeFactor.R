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
capture.output(print(anovaFearIma), file = "supplement/01s_fearRating_timeFactor_ima_anovaFreq.doc")

# bayesian ANOVA on fear ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFFearIma <- anovaBF(
  formula = fear ~ CS*time + partInd,
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFearIma)
capture.output(print(anovaBFFearIma), file = "supplement/01s_fearRating_timeFactor_ima_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for fear ratings in imagery-based conditioning group
plotFearIma <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "ima",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearIma 
ggsave(plot = plotFearIma, filename = "supplement/01s_fearRating_timeFactor_ima_plot.jpg",
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
capture.output(tableFearIma, file = "supplement/01s_fearRating_timeFactor_ima_tTable.doc")



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
capture.output(print(anovaFearReal), file = "supplement/01s_fearRating_timeFactor_real_anovaFreq.doc")

# bayesian ANOVA on fear ratings in classical conditioning group
set.seed(rngSeed); anovaBFFearReal <- anovaBF(
  formula = fear ~ CS*time + partInd,
  data = dataFearLong[dataFearLong$usGroup == "real",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFearReal)
capture.output(print(anovaBFFearReal), file = "supplement/01s_fearRating_timeFactor_real_anovaBayes.doc")

# quick & dirty graph of CS Type x time ANOVA for fear ratings in classical conditioning group
plotFearReal <- ezPlot(
  data = dataFearLong[dataFearLong$usGroup == "real",],
  dv = fear,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotFearReal 
ggsave(plot = plotFearReal, filename = "supplement/01s_fearRating_timeFactor_real_plot.jpg",
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
capture.output(tableFearReal, file = "supplement/01s_fearRating_timeFactor_real_tTable.doc")



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
capture.output(print(anovaFear), file = "supplement/01s_fearRating_timeFactor_acrossGroups_anovaFreq.doc")

# bayesian ANOVA on fear ratings across conditioning groups
set.seed(rngSeed); anovaBFFear <- anovaBF(
  formula = fear ~ usGroup*CS*time + partInd,
  data = dataFearLong,
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFear)
capture.output(print(anovaBFFear), file = "supplement/01s_fearRating_timeFactor_acrossGroups_anovaBayes.doc")

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
ggsave(plot = plotFear, filename = "Supplement/01s_fearRating_timeFactor_acrossGroups_plot.jpg",
       width = 20, height = 10, units = "cm")