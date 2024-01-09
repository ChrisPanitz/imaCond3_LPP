# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Feb 2022
# --- content: IBI analyses of imagery-based conditioning data in Panitz & Mueller (2022)

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

########################
### data preparation ###
########################

# Rating data contains group membership
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname,"/experimentData/imaCond3_demographicsAndRatings.txt"), sep = ",")

# load preprocessed IBI data (see readme file for explanation of content & structure)
importIBI <- read.csv(paste0(pathname,"/experimentData/imaCond3_ibimatrix_cs.txt"), sep = "\t")
importIBI <- merge(importRatings[,c("partCode","group")], importIBI, by = "partCode")
importIBI <- importIBI[order(importIBI$partCode),]
importIBI$partInd <- factor(1:48)

# create IBI data frame; IBI values are means from 2 to 5 s post-CS
dataIBI <- data.frame(
  partInd = factor(importIBI$partInd),
  usGroup = factor(importIBI$group, labels = c("ima", "real")),
  #Av_allTr = rowMeans(subset(importIBI, select = ibi_Akq_51_5:ibi_Akq_51_10)),
  Av_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_51_5:ibi_Akq1_51_10)),
  Av_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_51_5:ibi_Akq2_51_10)),
  #Neu_allTr = rowMeans(subset(importIBI, select = ibi_Akq_52_5:ibi_Akq_52_10)),
  Neu_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_52_5:ibi_Akq1_52_10)),
  Neu_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_52_5:ibi_Akq2_52_10)),
  #Min_allTr = rowMeans(subset(importIBI, select = ibi_Akq_53_5:ibi_Akq_53_10)),
  Min_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_53_5:ibi_Akq1_53_10)),
  Min_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_53_5:ibi_Akq2_53_10))
)  

# transform into long format
dataIBILong <- gather(data = dataIBI, key = "cond", value = "IBI", Av_1stBl:Min_2ndBl)
dataIBILong <- separate(data = dataIBILong, col = cond, into = c("CS","time"), sep = "_")
dataIBILong$CS <- factor(dataIBILong$CS, levels = c("Av","Neu","Min"))
dataIBILong$time <- factor(dataIBILong$time, levels = c("1stBl","2ndBl"))



############################################################################
### Imagery-based conditioning - supplementary analyses with time factor ###
############################################################################

# descriptive statistics for IBI in imagery-based conditioning group
describe(dataIBI[dataIBI$usGroup == "ima",])

# frequentist CS x Time ANOVA on IBI in imagery-based conditioning group, including p. eta^2
# IV = CS, Time; DV = IBI
anovaIBIIma <- ezANOVA(
  data = dataIBILong[dataIBILong$usGroup == "ima",],
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaIBIIma$ANOVA$pEtaSq <- 
  c(anovaIBIIma$ANOVA$SSn[1] / (anovaIBIIma$ANOVA$SSd[1]+anovaIBIIma$ANOVA$SSn[1]),
    anovaIBIIma$ANOVA$SSn[2] / (anovaIBIIma$ANOVA$SSd[2]+anovaIBIIma$ANOVA$SSn[2]),
    anovaIBIIma$ANOVA$SSn[3] / (anovaIBIIma$ANOVA$SSd[3]+anovaIBIIma$ANOVA$SSn[3]),
    anovaIBIIma$ANOVA$SSn[4] / (anovaIBIIma$ANOVA$SSd[4]+anovaIBIIma$ANOVA$SSn[4])
  ); print(anovaIBIIma)
capture.output(print(anovaIBIIma), file = "supplement/03s_ibi_timeFactor_ima_anovaFreq.doc")

# bayesian CS x Time ANOVA on IBI in imagery-based conditioning group
set.seed(rngSeed); anovaBFIBIIma <- generalTestBF(
  formula = IBI ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataIBILong[dataIBILong$usGroup == "ima",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
); print(anovaBFIBIIma)
capture.output(print(anovaBFIBIIma), file = paste0(pathname, "/supplement/03s_ibi_timeFactor_ima_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFIBIIma)
capture.output(bf_inclusion(anovaBFIBIIma), file = paste0(pathname, "/supplement/03s_ibi_timeFactor_ima_BFinclusion.doc"))

# quick graph of CS Type x Time ANOVA for IBI in imagery-based conditioning group
plotIBIIma <- ezPlot(
  data = dataIBILong[dataIBILong$usGroup == "ima",],
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotIBIIma 
ggsave(plot = plotIBIIma, filename = "supplement/03s_ibi_timeFactor_ima_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on IBI in imagery-based conditioning group
### 1stBl
# CS+av vs CS+neu
ibiImaAvNeu1stBl_t <- t.test(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                        y = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
ibiImaAvNeu1stBl_d <- cohens_d(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                          y = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                          paired = TRUE)
ibiImaAvNeu1stBl_BF <- ttestBF(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                          y = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiImaAvMin1stBl_t <- t.test(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                        y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
ibiImaAvMin1stBl_d <- cohens_d(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                          y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                          paired = TRUE)
ibiImaAvMin1stBl_BF <- ttestBF(x = dataIBI$Av_1stBl[dataIBI$usGroup == "ima"],
                          y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiImaNeuMin1stBl_t <- t.test(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                         y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
ibiImaNeuMin1stBl_d <- cohens_d(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                           y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                           paired = TRUE)
ibiImaNeuMin1stBl_BF <- ttestBF(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "ima"],
                           y = dataIBI$Min_1stBl[dataIBI$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
ibiImaAvNeu2ndBl_t <- t.test(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                             y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiImaAvNeu2ndBl_d <- cohens_d(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                               y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                               paired = TRUE)
ibiImaAvNeu2ndBl_BF <- ttestBF(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                               y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiImaAvMin2ndBl_t <- t.test(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                             y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiImaAvMin2ndBl_d <- cohens_d(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                               y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                               paired = TRUE)
ibiImaAvMin2ndBl_BF <- ttestBF(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "ima"],
                               y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiImaNeuMin2ndBl_t <- t.test(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                              y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                              alternative = "two.sided", paired = TRUE) # two-sided
ibiImaNeuMin2ndBl_d <- cohens_d(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                                y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                                paired = TRUE)
ibiImaNeuMin2ndBl_BF <- ttestBF(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "ima"],
                                y = dataIBI$Min_2ndBl[dataIBI$usGroup == "ima"],
                                nullInterval = NULL, paired = TRUE) # two-sided

tableIBIIma <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(ibiImaAvNeu1stBl_t$statistic, ibiImaAvMin1stBl_t$statistic, ibiImaNeuMin1stBl_t$statistic,
        ibiImaAvNeu2ndBl_t$statistic, ibiImaAvMin2ndBl_t$statistic, ibiImaNeuMin2ndBl_t$statistic),
  df = c(ibiImaAvNeu1stBl_t$parameter, ibiImaAvMin1stBl_t$parameter, ibiImaNeuMin1stBl_t$parameter,
         ibiImaAvNeu2ndBl_t$parameter, ibiImaAvMin2ndBl_t$parameter, ibiImaNeuMin2ndBl_t$parameter), 
  p = c(ibiImaAvNeu1stBl_t$p.value, ibiImaAvMin1stBl_t$p.value, ibiImaNeuMin1stBl_t$p.value,
        ibiImaAvNeu2ndBl_t$p.value, ibiImaAvMin2ndBl_t$p.value, ibiImaNeuMin2ndBl_t$p.value),
  d = c(ibiImaAvNeu1stBl_d$Cohens_d, ibiImaAvMin1stBl_d$Cohens_d, ibiImaNeuMin1stBl_d$Cohens_d,
        ibiImaAvNeu2ndBl_d$Cohens_d, ibiImaAvMin2ndBl_d$Cohens_d, ibiImaNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(ibiImaAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(ibiImaAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(ibiImaNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(ibiImaAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(ibiImaAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(ibiImaNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableIBIIma, file = "supplement/03s_ibi_timeFactor_ima_tTable.doc")



########################################################################
### Classical conditioning - supplementary analyses with time factor ###
########################################################################

# descriptive statistics for IBI in classical conditioning group
describe(dataIBI[dataIBI$usGroup == "real",])

# frequentist CS x Time ANOVA on IBI in classical conditioning group, including p. eta^2
# IV = CS, Time; DV = IBI
anovaIBIReal <- ezANOVA(
  data = dataIBILong[dataIBILong$usGroup == "real",],
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaIBIReal$ANOVA$pEtaSq <- 
  c(anovaIBIReal$ANOVA$SSn[1] / (anovaIBIReal$ANOVA$SSd[1]+anovaIBIReal$ANOVA$SSn[1]),
    anovaIBIReal$ANOVA$SSn[2] / (anovaIBIReal$ANOVA$SSd[2]+anovaIBIReal$ANOVA$SSn[2]),
    anovaIBIReal$ANOVA$SSn[3] / (anovaIBIReal$ANOVA$SSd[3]+anovaIBIReal$ANOVA$SSn[3]),
    anovaIBIReal$ANOVA$SSn[4] / (anovaIBIReal$ANOVA$SSd[4]+anovaIBIReal$ANOVA$SSn[4])
  ); print(anovaIBIReal)
capture.output(print(anovaIBIReal), file = "supplement/03s_ibi_timeFactor_real_anovaFreq.doc")

# bayesian CS x Time ANOVA on IBI in classical conditioning group
set.seed(rngSeed); anovaBFIBIReal <- generalTestBF(
  formula = IBI ~ CS*time + partInd + partInd:CS + partInd:time,
  data = dataIBILong[dataIBILong$usGroup == "real",],
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 100000
); print(anovaBFIBIReal)
capture.output(print(anovaBFIBIReal), file = paste0(pathname, "/supplement/03s_ibi_timeFactor_real_anovaBayes.doc"))

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFIBIReal)
capture.output(bf_inclusion(anovaBFIBIReal), file = paste0(pathname, "/supplement/03s_ibi_timeFactor_real_BFinclusion.doc"))

# quick graph of CS Type x Time ANOVA for IBI in imagery-based conditioning group
plotIBIReal <- ezPlot(
  data = dataIBILong[dataIBILong$usGroup == "real",],
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotIBIReal
ggsave(plot = plotIBIReal, filename = "supplement/03s_ibi_timeFactor_real_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on IBI in classical conditioning group
### 1stBl
# CS+av vs CS+neu
ibiRealAvNeu1stBl_t <- t.test(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                             y = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiRealAvNeu1stBl_d <- cohens_d(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                               paired = TRUE)
ibiRealAvNeu1stBl_BF <- ttestBF(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiRealAvMin1stBl_t <- t.test(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                             y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiRealAvMin1stBl_d <- cohens_d(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                               paired = TRUE)
ibiRealAvMin1stBl_BF <- ttestBF(x = dataIBI$Av_1stBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiRealNeuMin1stBl_t <- t.test(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                              y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                              alternative = "two.sided", paired = TRUE) # two-sided
ibiRealNeuMin1stBl_d <- cohens_d(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                                y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                                paired = TRUE)
ibiRealNeuMin1stBl_BF <- ttestBF(x = dataIBI$Neu_1stBl[dataIBI$usGroup == "real"],
                                y = dataIBI$Min_1stBl[dataIBI$usGroup == "real"],
                                nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
ibiRealAvNeu2ndBl_t <- t.test(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                             y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiRealAvNeu2ndBl_d <- cohens_d(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                               paired = TRUE)
ibiRealAvNeu2ndBl_BF <- ttestBF(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiRealAvMin2ndBl_t <- t.test(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                             y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                             alternative = "greater", paired = TRUE) # one-sided
ibiRealAvMin2ndBl_d <- cohens_d(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                               paired = TRUE)
ibiRealAvMin2ndBl_BF <- ttestBF(x = dataIBI$Av_2ndBl[dataIBI$usGroup == "real"],
                               y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                               nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiRealNeuMin2ndBl_t <- t.test(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                              y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                              alternative = "two.sided", paired = TRUE) # two-sided
ibiRealNeuMin2ndBl_d <- cohens_d(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                                y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                                paired = TRUE)
ibiRealNeuMin2ndBl_BF <- ttestBF(x = dataIBI$Neu_2ndBl[dataIBI$usGroup == "real"],
                                y = dataIBI$Min_2ndBl[dataIBI$usGroup == "real"],
                                nullInterval = NULL, paired = TRUE) # two-sided

tableIBIReal <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(ibiRealAvNeu1stBl_t$statistic, ibiRealAvMin1stBl_t$statistic, ibiRealNeuMin1stBl_t$statistic,
        ibiRealAvNeu2ndBl_t$statistic, ibiRealAvMin2ndBl_t$statistic, ibiRealNeuMin2ndBl_t$statistic),
  df = c(ibiRealAvNeu1stBl_t$parameter, ibiRealAvMin1stBl_t$parameter, ibiRealNeuMin1stBl_t$parameter,
         ibiRealAvNeu2ndBl_t$parameter, ibiRealAvMin2ndBl_t$parameter, ibiRealNeuMin2ndBl_t$parameter), 
  p = c(ibiRealAvNeu1stBl_t$p.value, ibiRealAvMin1stBl_t$p.value, ibiRealNeuMin1stBl_t$p.value,
        ibiRealAvNeu2ndBl_t$p.value, ibiRealAvMin2ndBl_t$p.value, ibiRealNeuMin2ndBl_t$p.value),
  d = c(ibiRealAvNeu1stBl_d$Cohens_d, ibiRealAvMin1stBl_d$Cohens_d, ibiRealNeuMin1stBl_d$Cohens_d,
        ibiRealAvNeu2ndBl_d$Cohens_d, ibiRealAvMin2ndBl_d$Cohens_d, ibiRealNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(ibiRealAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(ibiRealAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(ibiRealNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(ibiRealAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(ibiRealAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(ibiRealNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableIBIReal, file = "Supplement/03s_ibi_timeFactor_real_tTable.doc")



###############################################################
### Across groups - supplementary analyses with time factor ###
###############################################################

# descriptive statistics for IBI ratings across conditioning groups
describe(dataIBI)

# frequentist Group x CS x Time ANOVA on IBI across conditioning groups
anovaIBI <- ezANOVA(
  data = dataIBILong,
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaIBI$ANOVA$pEtaSq <- c(
  anovaIBI$ANOVA$SSn[1] / (anovaIBI$ANOVA$SSd[1]+anovaIBI$ANOVA$SSn[1]),
  anovaIBI$ANOVA$SSn[2] / (anovaIBI$ANOVA$SSd[2]+anovaIBI$ANOVA$SSn[2]),
  anovaIBI$ANOVA$SSn[3] / (anovaIBI$ANOVA$SSd[3]+anovaIBI$ANOVA$SSn[3]),
  anovaIBI$ANOVA$SSn[4] / (anovaIBI$ANOVA$SSd[4]+anovaIBI$ANOVA$SSn[4]),
  anovaIBI$ANOVA$SSn[5] / (anovaIBI$ANOVA$SSd[5]+anovaIBI$ANOVA$SSn[5]),
  anovaIBI$ANOVA$SSn[6] / (anovaIBI$ANOVA$SSd[6]+anovaIBI$ANOVA$SSn[6]),
  anovaIBI$ANOVA$SSn[7] / (anovaIBI$ANOVA$SSd[7]+anovaIBI$ANOVA$SSn[7]),
  anovaIBI$ANOVA$SSn[8] / (anovaIBI$ANOVA$SSd[8]+anovaIBI$ANOVA$SSn[8])
); print(anovaIBI)
capture.output(print(anovaIBI), file = "Supplement/03s_ibi_acrossGroups_anovaFreq.doc")

# bayesian ANOVA on IBI across conditioning groups
# set.seed(rngSeed); anovaBFIBI <- anovaBF(
#   formula = IBI ~ usGroup*CS*time + partInd,
#   data = dataIBILong,
#   whichRandom = "partInd",
#   whichModels = "all",
#   iterations = 100000
# ); print(anovaBFIBI)
set.seed(rngSeed); anovaBFIBI <- generalTestBF(
  formula = IBI ~ usGroup*CS*time + partInd + partInd:CS + partInd:time,
  data = dataIBILong,
  whichRandom = c("partInd", "partInd:CS", "partInd:time"),
  neverExclude = c("partInd", "partInd:CS", "partInd:time"),
  whichModels = "all",
  iterations = 10000 # only 10,000 iterations because it has to compute 128 models
); print(anovaBFIBI)
capture.output(print(anovaBFIBI), file = "Supplement/03s_ibi_timeFactor_acrossGroups_anovaBayes.doc")

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFIBI)
capture.output(bf_inclusion(anovaBFIBI), file = paste0(pathname, "/supplement/03s_ibi_timeFactor_acrossGroups_BFinclusion.doc"))

# quick graph of US Group x CS Type x Time ANOVA for IBI across groups
plotIBI <- ezPlot(
  data = dataIBILong,
  dv = IBI,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,
  col = usGroup
) ; plotIBI
ggsave(plot = plotIBI, filename = "supplement/03s_ibi_timeFactor_both_plot.jpg",
       width = 20, height = 10, units = "cm")
