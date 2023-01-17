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
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# load preprocessed IBI data (see readme file for explanation of content & structure)
importIBI <- read.csv(paste0(pathname,"/experimentData/imaCond3_ibimatrix_cs.txt"), sep="")
importIBI <- merge(importRatings[,c("partCode","group")], importIBI, by = "partCode")
importIBI <- importIBI[order(importIBI$partCode),]

# create IBI data frame; IBI values are means from 2 to 5 s post-CS
dataIBI <- data.frame(
  partInd = factor(1:dim(importIBI)[1]),
  usGroup = factor(importIBI$group, labels = c("ima", "real")),
  Av_allTr = rowMeans(subset(importIBI, select = ibi_Akq_51_5:ibi_Akq_51_10)),
  Av_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_51_5:ibi_Akq1_51_10)),
  Av_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_51_5:ibi_Akq2_51_10)),
  Neu_allTr = rowMeans(subset(importIBI, select = ibi_Akq_52_5:ibi_Akq_52_10)),
  Neu_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_52_5:ibi_Akq1_52_10)),
  Neu_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_52_5:ibi_Akq2_52_10)),
  Min_allTr = rowMeans(subset(importIBI, select = ibi_Akq_53_5:ibi_Akq_53_10)),
  Min_1stBl = rowMeans(subset(importIBI, select = ibi_Akq1_53_5:ibi_Akq1_53_10)),
  Min_2ndBl = rowMeans(subset(importIBI, select = ibi_Akq2_53_5:ibi_Akq2_53_10))
)  

# transform into long format
dataIBILong <- gather(data = dataIBI, key = "cond", value = "IBI", Av_allTr:Min_2ndBl)
dataIBILong <- separate(data = dataIBILong, col = cond, into = c("CS","time"), sep = "_")
dataIBILong$CS <- factor(dataIBILong$CS)
dataIBILong$time <- factor(dataIBILong$time)



#####################################################
### Imagery-based conditioning - primary analyses ###
#####################################################

# descriptive statistics for IBI in imagery-based conditioning group
describe(dataIBI[dataIBI$usGroup == "ima",])

# frequentist ANOVA on IBI in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = IBI
anovaIBIIma <- ezANOVA(
  data = dataIBILong[dataIBILong$usGroup == "ima" & dataIBILong$time == "allTr",],
  dv = IBI,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaIBIIma$ANOVA$pEtaSq <- 
  c(anovaIBIIma$ANOVA$SSn[1] / (anovaIBIIma$ANOVA$SSd[1]+anovaIBIIma$ANOVA$SSn[1]),
    anovaIBIIma$ANOVA$SSn[2] / (anovaIBIIma$ANOVA$SSd[2]+anovaIBIIma$ANOVA$SSn[2])
  ); print(anovaIBIIma)

# bayesian ANOVA on IBI in imagery-based conditioning group
set.seed(rngSeed); anovaBFIBIIma <- anovaBF(
  formula = IBI ~ CS + partInd,
  data = dataIBILong[dataIBILong$usGroup == "ima" & dataIBILong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFIBIIma)

# frequentist & bayesian t-tests on IBI in imagery-based conditioning group
# CS+av vs CS+neu
ibiImaAvNeu_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                        y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
ibiImaAvNeu_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                          y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                          paired = TRUE)
ibiImaAvNeu_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                          y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiImaAvMin_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                        y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
ibiImaAvMin_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                          y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                          paired = TRUE)
ibiImaAvMin_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "ima"],
                          y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiImaNeuMin_t <- t.test(x = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                         y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
ibiImaNeuMin_d <- cohens_d(x = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                           y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                           paired = TRUE)
ibiImaNeuMin_BF <- ttestBF(x = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                           y = dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided



#################################################
### Classical conditioning - primary analyses ###
#################################################

# descriptive statistics for IBI in classical conditioning group
describe(dataIBI[dataIBI$usGroup == "real",])

# frequentist ANOVA on IBI in classical conditioning group, including p. eta^2
# IV = CS; DV = IBI
anovaIBIReal <- ezANOVA(
  data = dataIBILong[dataIBILong$usGroup == "real" & dataIBILong$time == "allTr",],
  dv = IBI,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaIBIReal$ANOVA$pEtaSq <- 
  c(anovaIBIReal$ANOVA$SSn[1] / (anovaIBIReal$ANOVA$SSd[1]+anovaIBIReal$ANOVA$SSn[1]),
    anovaIBIReal$ANOVA$SSn[2] / (anovaIBIReal$ANOVA$SSd[2]+anovaIBIReal$ANOVA$SSn[2])
  ); print(anovaIBIReal)

# bayesian ANOVA on IBI in classical conditioning group
set.seed(rngSeed); anovaBFIBIReal <- anovaBF(
  formula = IBI ~ CS + partInd,
  data = dataIBILong[dataIBILong$usGroup == "real" & dataIBILong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFIBIReal)

# frequentist & bayesian t-tests on IBI in classical conditioning group
# CS+av vs CS+neu
ibiRealAvNeu_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                         y = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
ibiRealAvNeu_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                           paired = TRUE)
ibiRealAvNeu_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
ibiRealAvMin_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                         y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
ibiRealAvMin_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                           paired = TRUE)
ibiRealAvMin_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
ibiRealNeuMin_t <- t.test(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                          y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                          alternative = "two.sided", paired = TRUE) # two-sided
ibiRealNeuMin_d <- cohens_d(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                            y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                            paired = TRUE)
ibiRealNeuMin_BF <- ttestBF(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                            y = dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # two-sided



########################################
### Across groups - primary analyses ###
########################################

# descriptive statistics for IBI ratings across conditioning groups
describe(dataIBI)

# frequentist ANOVA on IBI across conditioning groups
anovaIBI <- ezANOVA(
  data = dataIBILong[dataIBILong$time == "allTr",],
  dv = IBI,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaIBI$ANOVA$pEtaSq <- c(anovaIBI$ANOVA$SSn[1] /
                                (anovaIBI$ANOVA$SSd[1]+anovaIBI$ANOVA$SSn[1]),
                              anovaIBI$ANOVA$SSn[2] /
                                (anovaIBI$ANOVA$SSd[2]+anovaIBI$ANOVA$SSn[2]),
                              anovaIBI$ANOVA$SSn[3] /
                                (anovaIBI$ANOVA$SSd[3]+anovaIBI$ANOVA$SSn[3]),
                              anovaIBI$ANOVA$SSn[4] /
                                (anovaIBI$ANOVA$SSd[4]+anovaIBI$ANOVA$SSn[4])
); print(anovaIBI)

# bayesian ANOVA on IBI across conditioning groups
set.seed(rngSeed); anovaBFIBI <- anovaBF(
  formula = IBI ~ usGroup*CS + partInd,
  data = dataIBILong[dataIBILong$time == "allTr",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFIBI)

# inclusion factors for bayesian ANOVA effects
bf_nullModel <- 1
bf_usGroup <- exp(anovaBFIBI@bayesFactor$bf[1])
bf_cs <- exp(anovaBFIBI@bayesFactor$bf[2])
bf_interact <- exp(anovaBFIBI@bayesFactor$bf[3])
bf_usGroup_cs <- exp(anovaBFIBI@bayesFactor$bf[4])
bf_usGroup_interact <- exp(anovaBFIBI@bayesFactor$bf[5])
bf_cs_interact <- exp(anovaBFIBI@bayesFactor$bf[6])
bf_fullModel <- exp(anovaBFIBI@bayesFactor$bf[7])

# main effect US group: models [1] and [3] vs. null model and model [2]
bfIncGroupIBI <- (bf_usGroup + bf_usGroup_cs + bf_usGroup_interact + bf_fullModel) / 
  (bf_nullModel + bf_cs + bf_interact + bf_cs_interact); bfIncGroupIBI
# main effect CS type: models "main effect CS" & "main effects CS & group" vs.
#                      null model and "main effect group"
bfIncCsIBI <- (bf_cs + bf_usGroup_cs + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_interact + bf_usGroup_interact); bfIncCsIBI
# interaction: Full model vs. main-effects-only model
bfIncInteractIBI <- (bf_interact + bf_usGroup_interact + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_cs + bf_usGroup_cs); bfIncInteractIBI

# quick & dirty graph of group x CS ANOVA on IBI
ezPlot(
  data = dataIBILong[dataIBILong$time == "allTr",],
  dv = IBI,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on IBI (difference scores) across groups
# delta [CS+av - CS+neu]
ibiBothAvNeu_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                             dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                         y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                             dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
ibiBothAvNeu_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                               dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                               dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                               paired = FALSE)
ibiBothAvNeu_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                               dataIBI$Neu_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                               dataIBI$Neu_allTr[dataIBI$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
ibiBothAvMin_t <- t.test(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                             dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                         y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                             dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
ibiBothAvMin_d <- cohens_d(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                               dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                               dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                               paired = FALSE)
ibiBothAvMin_BF <- ttestBF(x = dataIBI$Av_allTr[dataIBI$usGroup == "real"] -
                               dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                           y = dataIBI$Av_allTr[dataIBI$usGroup == "ima"] -
                               dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
ibiBothNeuMin_t <- t.test(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"] -
                              dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                          y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"] -
                              dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
ibiBothNeuMin_d <- cohens_d(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"] -
                                dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                            y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"] -
                                dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                                paired = FALSE)
ibiBothNeuMin_BF <- ttestBF(x = dataIBI$Neu_allTr[dataIBI$usGroup == "real"] -
                                dataIBI$Min_allTr[dataIBI$usGroup == "real"],
                            y = dataIBI$Neu_allTr[dataIBI$usGroup == "ima"] -
                                dataIBI$Min_allTr[dataIBI$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided



# quick & dirty graph of group x CS ANOVA on valence ratings
ezPlot(
  data = dataIBILong,
  dv = IBI,
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
  t = c(ibiImaAvNeu_t$statistic, ibiImaAvMin_t$statistic, ibiImaNeuMin_t$statistic,
        ibiRealAvNeu_t$statistic, ibiRealAvMin_t$statistic, ibiRealNeuMin_t$statistic,
        ibiBothAvNeu_t$statistic, ibiBothAvMin_t$statistic, ibiBothNeuMin_t$statistic),
  df = c(ibiImaAvNeu_t$parameter, ibiImaAvMin_t$parameter, ibiImaNeuMin_t$parameter,
         ibiRealAvNeu_t$parameter, ibiRealAvMin_t$parameter, ibiRealNeuMin_t$parameter,
         ibiBothAvNeu_t$parameter, ibiBothAvMin_t$parameter, ibiBothNeuMin_t$parameter), 
  p = c(ibiImaAvNeu_t$p.value, ibiImaAvMin_t$p.value, ibiImaNeuMin_t$p.value,
        ibiRealAvNeu_t$p.value, ibiRealAvMin_t$p.value, ibiRealNeuMin_t$p.value,
        ibiBothAvNeu_t$p.value*3, ibiBothAvMin_t$p.value*3, ibiBothNeuMin_t$p.value*3),  # Bonferroni
  d = c(ibiImaAvNeu_d$Cohens_d, ibiImaAvMin_d$Cohens_d, ibiImaNeuMin_d$Cohens_d,
        ibiRealAvNeu_d$Cohens_d, ibiRealAvMin_d$Cohens_d, ibiRealNeuMin_d$Cohens_d,
        ibiBothAvNeu_d$Cohens_d, ibiBothAvMin_d$Cohens_d, ibiBothNeuMin_d$Cohens_d),
  BF = c(exp(ibiImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(ibiImaAvMin_BF@bayesFactor[["bf"]][1]), exp(ibiImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(ibiRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(ibiRealAvMin_BF@bayesFactor[["bf"]][1]), exp(ibiRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(ibiBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(ibiBothAvMin_BF@bayesFactor[["bf"]][1]), exp(ibiBothNeuMin_BF@bayesFactor[["bf"]][1]))
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

tableIBI <- flextable(tableData[1:9,])
tableIBI <- add_header_lines(tableIBI, top = TRUE, values = "IBI")
tableIBI <- align(tableIBI, align = "center")

save_as_docx(tableIBI, path = "Tables/tableIBI_raw.docx")



####################
### Plotting IBI ###
####################

csLabels = c("CS+av", "CS+neu", "CS-")
ibi4GA <- importIBI
ibi4GAima <- importIBI[importIBI$group == 0,]
ibi4GAreal <- importIBI[importIBI$group == 1,]

ibiGAacq <- data.frame(
              csType = factor(c(rep(1,16), rep(2,16), rep(3,16)), labels = csLabels),
              time = rep(seq(-.75,6.75,.5),3),
              IBI = c(colMeans(subset(ibi4GA, select = ibi_Akq_51_bl1:ibi_Akq_51_14)),
                      colMeans(subset(ibi4GA, select = ibi_Akq_52_bl1:ibi_Akq_52_14)),
                      colMeans(subset(ibi4GA, select = ibi_Akq_53_bl1:ibi_Akq_53_14))),
              row.names = NULL)
ibiGAacqIma <- data.frame(
              csType = factor(c(rep(1,16), rep(2,16), rep(3,16)), labels = csLabels),
              time = rep(seq(-.75,6.75,.5),3),
              IBI = c(colMeans(subset(ibi4GAima, select = ibi_Akq_51_bl1:ibi_Akq_51_14)),
                      colMeans(subset(ibi4GAima, select = ibi_Akq_52_bl1:ibi_Akq_52_14)),
                      colMeans(subset(ibi4GAima, select = ibi_Akq_53_bl1:ibi_Akq_53_14))),
              row.names = NULL)
ibiGAacqReal <- data.frame(
              csType = factor(c(rep(1,16), rep(2,16), rep(3,16)), labels = csLabels),
              time = rep(seq(-.75,6.75,.5),3),
              IBI = c(colMeans(subset(ibi4GAreal, select = ibi_Akq_51_bl1:ibi_Akq_51_14)),
                      colMeans(subset(ibi4GAreal, select = ibi_Akq_52_bl1:ibi_Akq_52_14)),
                      colMeans(subset(ibi4GAreal, select = ibi_Akq_53_bl1:ibi_Akq_53_14))),
              row.names = NULL)


dataIBIWithin <- dataIBI[,c("partInd","usGroup","Av_allTr","Neu_allTr","Min_allTr")]
# remove each participant's average from each single value
dataIBIWithin[,3:5] <- as.matrix(dataIBIWithin[,3:5]) -
  rowMeans(as.matrix(dataIBIWithin[,3:5])) 
# prepare data frame for bar plot with means from standard dataset and SE from
# dataset without betweem-subject variance
meanIBI <- data.frame(
  usGroup = factor(c(rep(1,3),rep(2,3)),
                   labels = c("Imagery-Based","Classical")),
  CS = factor(c(1,2,3,1,2,3),
              labels = c("CS+\nav","CS+\nneu","CS-\n")),
  mean = c(describe(dataIBI[dataIBI$usGroup == "ima", c(3,6,9)])$mean,
           describe(dataIBI[dataIBI$usGroup == "real", c(3,6,9)])$mean),
  se = c(describe(dataIBIWithin[dataIBIWithin$usGroup == "ima", 3:5])$se,
         describe(dataIBIWithin[dataIBIWithin$usGroup == "real", 3:5])$se)
)




lineSize = 1
yMin = -20
yMax = 42
plotFS <- 8
showSig <- TRUE

graphIBIima <- ggplot(data = ibiGAacqIma, aes(x = time, y = IBI, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = 2, xmax = 5, ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), size = lineSize) + 
  scale_x_continuous(breaks = seq(-1,7,1)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Imagery-Based Conditioning", x = "Time (s)", y = "IBI change (ms)", fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = FALSE) +
  theme(#legend.box = "horizontal",
        legend.position = "none",
        #legend.text = element_text(size = plotFS, color = "black"),
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

graphIBIreal <- ggplot(data = ibiGAacqReal, aes(x = time, y = IBI, colour = csType)) + 
  theme_classic() +
  geom_rect(xmin = 2, xmax = 5, ymin = yMin, ymax = yMax, fill = "gray90", colour = NA) +
  geom_line(aes(colour = csType), size = lineSize) + 
  scale_x_continuous(breaks = seq(-1,7,1)) +
  scale_colour_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  lims(y = c(yMin, yMax)) +
  labs(title = "Classical Conditioning", x = "Time (s)", y = "IBI change (ms)", fill = "", colour = "") +
  guides(colour = guide_legend(order = 1), fill = FALSE) +
  theme(#legend.box = "horizontal",
        legend.position = "none",
        #legend.text = element_text(size = plotFS, color = "black"),
        plot.title = element_text(size = plotFS, color = "black", face = "bold", hjust = .5),
        axis.title.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.text.x = element_text(margin = margin(t = 5), size = plotFS, color = "black"),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"))

graphIBImeans <- ggplot(data = meanIBI, aes(x = usGroup, y = mean, fill = factor(CS))) +
  theme_classic() +
  geom_col(aes(fill = CS), position = position_dodge(width = .9)) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = .1), position = position_dodge(width = .9)) +
  scale_x_discrete(aes(breaks = usGroup), name = "") +
  scale_y_continuous(name = "Mean IBI change (2-5 s)") +
  geom_hline(yintercept = 0) +
  #geom_text(aes(label = CS, y = 11.5), position = position_dodge(.9), colour = "black", size = plotFS/.pt) + 
  geom_text(aes(label = usGroup, y = 12.5), colour = "black", size = (plotFS/.pt)-.5, fontface = "bold") +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"))

# combining graphs into one figure
graphIBIima <- graphIBIima + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphIBIreal <- graphIBIreal + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphIBImeans <- graphIBImeans + theme(plot.margin = unit(c(5,5,5,5), "mm"))
graphIBI <- ggarrange(graphIBIima, graphIBIreal, graphIBImeans,
                          labels = c("A", "B", "C"),
                          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)
graphIBIupperRow <- ggarrange(graphIBIima, graphIBIreal,
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom"
)
graphIBI <- ggarrange(graphIBIupperRow, graphIBImeans,
                          labels = c("", "C"),
                          ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom"
)
graphIBI

# saving it
ggsave(filename = "Figures/Figure4_timeCourses_barPlot_IBITEST4.eps",
       plot = graphIBI,
       width = 200,
       height = 140,
       #height = 200,
       units = "mm",
       dpi = 300
)
