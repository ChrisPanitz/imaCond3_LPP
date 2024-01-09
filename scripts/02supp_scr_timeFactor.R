# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Feb 2022
# --- content: secondary SCR analyses of imagery-based conditioning data in Panitz & Mueller (2022)

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

# Rating data contains group membership
# (see imaCond3_demographicsAndRatings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# load SCR data from text file
# (see imaCond3_scr_readme.txt for more details)
importSCR <- read.csv(paste0(pathname, "/experimentData/imaCond3_scrmatrix_cs.txt"), sep = "")
importSCR <- merge(importRatings[,c("partCode","group")], importSCR, by = "partCode")
importSCR <- importSCR[order(importSCR$partCode),]

# create data frames in wide & long format for SCR
dataSCR <- data.frame(
  partInd = factor(1:dim(importSCR)[1]),
  usGroup = factor(importSCR$group, labels = c("ima", "real")),
  #Av_allTr = importSCR$scr_Akqall_51_norm,
  Av_1stBl = importSCR$scr_Akq1_51_norm,
  Av_2ndBl = importSCR$scr_Akq2_51_norm,
  #Neu_allTr = importSCR$scr_Akqall_52_norm,
  Neu_1stBl = importSCR$scr_Akq1_52_norm,
  Neu_2ndBl = importSCR$scr_Akq2_52_norm,
  #Min_allTr = importSCR$scr_Akqall_53_norm,
  Min_1stBl = importSCR$scr_Akq1_53_norm,
  Min_2ndBl = importSCR$scr_Akq2_53_norm
)  
dataSCRLong <- gather(data = dataSCR, key = "cond", value = "SCR",
                      Av_1stBl:Min_2ndBl)
dataSCRLong <- separate(data = dataSCRLong, col = cond,
                        into = c("CS","time"), sep = "_")
dataSCRLong$CS <- factor(dataSCRLong$CS, levels = c("Av","Neu","Min"))
dataSCRLong$time <- factor(dataSCRLong$time, levels = c("1stBl","2ndBl"))


#######################################################
### Imagery-based conditioning - secondary analyses ###
#######################################################

# descriptive statistics for SCR in imagery-based conditioning group
describe(dataSCR[dataSCR$usGroup == "ima",])

# frequentist CS x Time ANOVA on SCR in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = SCR 
anovaSCRIma <- ezANOVA(
  data = dataSCRLong[dataSCRLong$usGroup == "ima",],
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaSCRIma$ANOVA$pEtaSq <- 
  c(anovaSCRIma$ANOVA$SSn[1] / (anovaSCRIma$ANOVA$SSd[1]+anovaSCRIma$ANOVA$SSn[1]),
    anovaSCRIma$ANOVA$SSn[2] / (anovaSCRIma$ANOVA$SSd[2]+anovaSCRIma$ANOVA$SSn[2]),
    anovaSCRIma$ANOVA$SSn[3] / (anovaSCRIma$ANOVA$SSd[3]+anovaSCRIma$ANOVA$SSn[3]),
    anovaSCRIma$ANOVA$SSn[4] / (anovaSCRIma$ANOVA$SSd[4]+anovaSCRIma$ANOVA$SSn[4])
  ); print(anovaSCRIma)
capture.output(print(anovaSCRIma), file = "supplement/02s_scr_timeFactor_ima_anovaFreq.doc")

# bayesian CS x Time ANOVA on SCR in imagery-based conditioning group
set.seed(rngSeed); anovaBFSCRIma <- anovaBF(
  formula = SCR ~ CS*time + partInd,
  data = dataSCRLong[dataSCRLong$usGroup == "ima",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFSCRIma)
capture.output(print(anovaBFSCRIma), file = "supplement/02s_scr_timeFactor_ima_anovaBayes.doc")

# quick graph of CS Type x Time ANOVA for SCR in imagery-based conditioning group
plotSCRIma <- ezPlot(
  data = dataSCRLong[dataSCRLong$usGroup == "ima",],
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
) ; plotSCRIma 
ggsave(plot = plotSCRIma, filename = "supplement/02s_scr_timeFactor_ima_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on SCR in imagery-based conditioning group
### 1stBl
# CS+av vs CS+neu
scrImaAvNeu1stBl_t <- t.test(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                        y = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvNeu1stBl_d <- cohens_d(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvNeu1stBl_BF <- ttestBF(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrImaAvMin1stBl_t <- t.test(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                        y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvMin1stBl_d <- cohens_d(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvMin1stBl_BF <- ttestBF(x = dataSCR$Av_1stBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrImaNeuMin1stBl_t <- t.test(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                         y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
scrImaNeuMin1stBl_d <- cohens_d(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                           paired = TRUE)
scrImaNeuMin1stBl_BF <- ttestBF(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_1stBl[dataSCR$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
scrImaAvNeu2ndBl_t <- t.test(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                        y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvNeu2ndBl_d <- cohens_d(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvNeu2ndBl_BF <- ttestBF(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrImaAvMin2ndBl_t <- t.test(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                        y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvMin2ndBl_d <- cohens_d(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvMin2ndBl_BF <- ttestBF(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrImaNeuMin2ndBl_t <- t.test(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                         y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
scrImaNeuMin2ndBl_d <- cohens_d(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                           paired = TRUE)
scrImaNeuMin2ndBl_BF <- ttestBF(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_2ndBl[dataSCR$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided

tableSCRIma <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(scrImaAvNeu1stBl_t$statistic, scrImaAvMin1stBl_t$statistic, scrImaNeuMin1stBl_t$statistic,
        scrImaAvNeu2ndBl_t$statistic, scrImaAvMin2ndBl_t$statistic, scrImaNeuMin2ndBl_t$statistic),
  df = c(scrImaAvNeu1stBl_t$parameter, scrImaAvMin1stBl_t$parameter, scrImaNeuMin1stBl_t$parameter,
         scrImaAvNeu2ndBl_t$parameter, scrImaAvMin2ndBl_t$parameter, scrImaNeuMin2ndBl_t$parameter), 
  p = c(scrImaAvNeu1stBl_t$p.value, scrImaAvMin1stBl_t$p.value, scrImaNeuMin1stBl_t$p.value,
        scrImaAvNeu2ndBl_t$p.value, scrImaAvMin2ndBl_t$p.value, scrImaNeuMin2ndBl_t$p.value),
  d = c(scrImaAvNeu1stBl_d$Cohens_d, scrImaAvMin1stBl_d$Cohens_d, scrImaNeuMin1stBl_d$Cohens_d,
        scrImaAvNeu2ndBl_d$Cohens_d, scrImaAvMin2ndBl_d$Cohens_d, scrImaNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(scrImaAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(scrImaAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(scrImaNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(scrImaAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(scrImaAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(scrImaNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableSCRIma, file = "supplement/02s_scr_timeFactor_ima_tTable.doc")

###################################################
### Classical conditioning - secondary analyses ###
###################################################

# descriptive statistics for SCR in classical conditioning group
describe(dataSCR[dataSCR$usGroup == "real",])

# frequentist CS x Time ANOVA on SCR in classical conditioning group, including p. eta^2
# IV = CS; DV = SCR
anovaSCRReal <- ezANOVA(
  data = dataSCRLong[dataSCRLong$usGroup == "real",],
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  type = 3,
  detailed = TRUE
); anovaSCRReal$ANOVA$pEtaSq <- 
  c(anovaSCRReal$ANOVA$SSn[1] / (anovaSCRReal$ANOVA$SSd[1]+anovaSCRReal$ANOVA$SSn[1]),
    anovaSCRReal$ANOVA$SSn[2] / (anovaSCRReal$ANOVA$SSd[2]+anovaSCRReal$ANOVA$SSn[2]),
    anovaSCRReal$ANOVA$SSn[3] / (anovaSCRReal$ANOVA$SSd[3]+anovaSCRReal$ANOVA$SSn[3]),
    anovaSCRReal$ANOVA$SSn[4] / (anovaSCRReal$ANOVA$SSd[4]+anovaSCRReal$ANOVA$SSn[4])
  ); print(anovaSCRReal)
capture.output(print(anovaSCRReal), file = "supplement/02s_scr_timeFactor_real_anovaFreq.doc")

# bayesian CS x Time ANOVA on SCR in classical conditioning group
set.seed(rngSeed); anovaBFSCRReal <- anovaBF(
  formula = SCR ~ CS*time + partInd,
  data = dataSCRLong[dataSCRLong$usGroup == "real",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFSCRReal)
capture.output(print(anovaBFSCRReal), file = "supplement/02s_scr_timeFactor_real_anovaBayes.doc")

# quick graph of CS Type x Time ANOVA for SCR in classical conditioning group
plotSCRReal <- ezPlot(
  data = dataSCRLong[dataSCRLong$usGroup == "real",],
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  x = time,
  split = CS
); plotSCRReal
ggsave(plot = plotSCRReal, filename = "supplement/02s_scr_timeFactor_real_plot.jpg",
       width = 10, height = 10, units = "cm")

# frequentist & bayesian t-tests on SCR in classical conditioning group
### 1stBl
# CS+av vs CS+neu
scrRealAvNeu1stBl_t <- t.test(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                         y = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
scrRealAvNeu1stBl_d <- cohens_d(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                           y = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                           paired = TRUE)
scrRealAvNeu1stBl_BF <- ttestBF(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                           y = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrRealAvMin1stBl_t <- t.test(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                         y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
scrRealAvMin1stBl_d <- cohens_d(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                           y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                           paired = TRUE)
scrRealAvMin1stBl_BF <- ttestBF(x = dataSCR$Av_1stBl[dataSCR$usGroup == "real"],
                           y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrRealNeuMin1stBl_t <- t.test(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                          y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                          alternative = "two.sided", paired = TRUE) # two-sided
scrRealNeuMin1stBl_d <- cohens_d(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                            y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                            paired = TRUE)
scrRealNeuMin1stBl_BF <- ttestBF(x = dataSCR$Neu_1stBl[dataSCR$usGroup == "real"],
                            y = dataSCR$Min_1stBl[dataSCR$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # two-sided

### 2ndBl
# CS+av vs CS+neu
scrRealAvNeu2ndBl_t <- t.test(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                              y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
scrRealAvNeu2ndBl_d <- cohens_d(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                                y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                                paired = TRUE)
scrRealAvNeu2ndBl_BF <- ttestBF(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                                y = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrRealAvMin2ndBl_t <- t.test(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                              y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                              alternative = "greater", paired = TRUE) # one-sided
scrRealAvMin2ndBl_d <- cohens_d(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                                y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                                paired = TRUE)
scrRealAvMin2ndBl_BF <- ttestBF(x = dataSCR$Av_2ndBl[dataSCR$usGroup == "real"],
                                y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                                nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrRealNeuMin2ndBl_t <- t.test(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                               y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                               alternative = "two.sided", paired = TRUE) # two-sided
scrRealNeuMin2ndBl_d <- cohens_d(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                                 y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                                 paired = TRUE)
scrRealNeuMin2ndBl_BF <- ttestBF(x = dataSCR$Neu_2ndBl[dataSCR$usGroup == "real"],
                                 y = dataSCR$Min_2ndBl[dataSCR$usGroup == "real"],
                                 nullInterval = NULL, paired = TRUE) # two-sided

tableSCRReal <- data.frame(
  time = c(rep("1stBl",3), rep("2ndBl",3)),
  comparison = rep(c("CS+av vs CS+neu", "CS+av vs CS-", "CSneu vs CS-"), 2),
  t = c(scrRealAvNeu1stBl_t$statistic, scrRealAvMin1stBl_t$statistic, scrRealNeuMin1stBl_t$statistic,
        scrRealAvNeu2ndBl_t$statistic, scrRealAvMin2ndBl_t$statistic, scrRealNeuMin2ndBl_t$statistic),
  df = c(scrRealAvNeu1stBl_t$parameter, scrRealAvMin1stBl_t$parameter, scrRealNeuMin1stBl_t$parameter,
         scrRealAvNeu2ndBl_t$parameter, scrRealAvMin2ndBl_t$parameter, scrRealNeuMin2ndBl_t$parameter), 
  p = c(scrRealAvNeu1stBl_t$p.value, scrRealAvMin1stBl_t$p.value, scrRealNeuMin1stBl_t$p.value,
        scrRealAvNeu2ndBl_t$p.value, scrRealAvMin2ndBl_t$p.value, scrRealNeuMin2ndBl_t$p.value),
  d = c(scrRealAvNeu1stBl_d$Cohens_d, scrRealAvMin1stBl_d$Cohens_d, scrRealNeuMin1stBl_d$Cohens_d,
        scrRealAvNeu2ndBl_d$Cohens_d, scrRealAvMin2ndBl_d$Cohens_d, scrRealNeuMin2ndBl_d$Cohens_d),
  BF = c(exp(scrRealAvNeu1stBl_BF@bayesFactor[["bf"]][1]), exp(scrRealAvMin1stBl_BF@bayesFactor[["bf"]][1]), exp(scrRealNeuMin1stBl_BF@bayesFactor[["bf"]][1]),
         exp(scrRealAvNeu2ndBl_BF@bayesFactor[["bf"]][1]), exp(scrRealAvMin2ndBl_BF@bayesFactor[["bf"]][1]), exp(scrRealNeuMin2ndBl_BF@bayesFactor[["bf"]][1])),
  testDir = rep(c("one.sided","one.sided","two.sided"),2)
)
capture.output(tableSCRReal, file = "supplement/02s_scr_timeFactor_real_tTable.doc")



###############################################################
### Across groups - supplementary analyses with time factor ###
###############################################################

# descriptive statistics for SCR ratings across conditioning groups
describe(dataSCR)

# frequentist Group x CS x Time ANOVA on SCR 
anovaSCR <- ezANOVA(
  data = dataSCRLong,
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaSCR$ANOVA$pEtaSq <- c(
  anovaSCR$ANOVA$SSn[1] / (anovaSCR$ANOVA$SSd[1]+anovaSCR$ANOVA$SSn[1]),
  anovaSCR$ANOVA$SSn[2] / (anovaSCR$ANOVA$SSd[2]+anovaSCR$ANOVA$SSn[2]),
  anovaSCR$ANOVA$SSn[3] / (anovaSCR$ANOVA$SSd[3]+anovaSCR$ANOVA$SSn[3]),
  anovaSCR$ANOVA$SSn[4] / (anovaSCR$ANOVA$SSd[4]+anovaSCR$ANOVA$SSn[4]),
  anovaSCR$ANOVA$SSn[5] / (anovaSCR$ANOVA$SSd[5]+anovaSCR$ANOVA$SSn[5]),
  anovaSCR$ANOVA$SSn[6] / (anovaSCR$ANOVA$SSd[6]+anovaSCR$ANOVA$SSn[6]),
  anovaSCR$ANOVA$SSn[7] / (anovaSCR$ANOVA$SSd[7]+anovaSCR$ANOVA$SSn[7]),
  anovaSCR$ANOVA$SSn[8] / (anovaSCR$ANOVA$SSd[8]+anovaSCR$ANOVA$SSn[8])
); print(anovaSCR)
capture.output(print(anovaSCR), file = "supplement/02s_scr_timeFactor_both_anovaFreq.doc")

# bayesian Group x CS x Time ANOVA on SCR 
set.seed(rngSeed); anovaBFSCR <- anovaBF(
  formula = SCR ~ usGroup*CS*time + partInd,
  data = dataSCRLong,
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFSCR)
capture.output(print(anovaBFSCR), file = "supplement/02s_scr_timeFactor_both_anovaBayes.doc")

# quick graph of Group x CS x Time ANOVA on SCR
plotSCR <- ezPlot(
  data = dataSCRLong,
  dv = SCR,
  wid = partInd,
  within = .(CS,time),
  between = .(usGroup),
  x = time,
  split = CS,
  col = usGroup
); plotSCR
ggsave(plot = plotSCR, filename = "supplement/02s_scr_timeFactor_both_plot.jpg",
       width = 20, height = 10, units = "cm")