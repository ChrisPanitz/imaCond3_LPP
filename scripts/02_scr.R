# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Feb 2022
# --- content: SCR analyses of imagery-based conditioning data in Panitz & Mueller (2022)

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
#importSCR <- merge(importRatings[,c("vpcode","partInd","group")], importSCR, by = "vpcode")
importSCR <- merge(importRatings[,c("partCode","group")], importSCR, by = "partCode")
importSCR <- importSCR[order(importSCR$partCode),]

# create data frames in wide & long format for unpleasantness ratings
dataSCR <- data.frame(
  partInd = factor(1:dim(importSCR)[1]),
  usGroup = factor(importSCR$group, labels = c("ima", "real")),
  Av_allTr = importSCR$scr_Akqall_51_norm,
  Av_1stBl = importSCR$scr_Akq1_51_norm,
  Av_2ndBl = importSCR$scr_Akq2_51_norm,
  Neu_allTr = importSCR$scr_Akqall_52_norm,
  Neu_1stBl = importSCR$scr_Akq1_52_norm,
  Neu_2ndBl = importSCR$scr_Akq2_52_norm,
  Min_allTr = importSCR$scr_Akqall_53_norm,
  Min_1stBl = importSCR$scr_Akq1_53_norm,
  Min_2ndBl = importSCR$scr_Akq2_53_norm
)  
dataSCRLong <- gather(data = dataSCR, key = "cond", value = "SCR",
                      Av_allTr:Min_2ndBl)
dataSCRLong <- separate(data = dataSCRLong, col = cond,
                        into = c("CS","time"), sep = "_")
dataSCRLong$CS <- factor(dataSCRLong$CS)
dataSCRLong$time <- factor(dataSCRLong$time)


#####################################################
### Imagery-based conditioning - primary analyses ###
#####################################################

# descriptive statistics for SCR in imagery-based conditioning group
describe(dataSCR[dataSCR$usGroup == "ima",])

# frequentist ANOVA on SCR in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = SCR 
anovaSCRIma <- ezANOVA(
  data = dataSCRLong[dataSCRLong$usGroup == "ima" & dataSCRLong$time == "allTr",],
  dv = SCR,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaSCRIma$ANOVA$pEtaSq <- 
  c(anovaSCRIma$ANOVA$SSn[1] / (anovaSCRIma$ANOVA$SSd[1]+anovaSCRIma$ANOVA$SSn[1]),
    anovaSCRIma$ANOVA$SSn[2] / (anovaSCRIma$ANOVA$SSd[2]+anovaSCRIma$ANOVA$SSn[2])
  ); print(anovaSCRIma)

#  bayesian ANOVA on SCR in imagery-based conditioning group
set.seed(rngSeed); anovaBFSCRIma <- anovaBF(
  formula = SCR ~ CS + partInd,
  data = dataSCRLong[dataSCRLong$usGroup == "ima" & dataSCRLong$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFSCRIma)

# frequentist & bayesian t-tests on SCR in imagery-based conditioning group
# CS+av vs CS+neu
scrImaAvNeu_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                        y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvNeu_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvNeu_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                          y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrImaAvMin_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                        y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                        alternative = "greater", paired = TRUE) # one-sided
scrImaAvMin_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                          paired = TRUE)
scrImaAvMin_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "ima"],
                          y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrImaNeuMin_t <- t.test(x = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                         y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = TRUE) # two-sided
scrImaNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                           paired = TRUE)
scrImaNeuMin_BF <- ttestBF(x = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                           y = dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                           nullInterval = NULL, paired = TRUE) # two-sided



#################################################
### Classical conditioning - primary analyses ###
#################################################

# descriptive statistics for SCR in classical conditioning group
describe(dataSCR[dataSCR$usGroup == "real",])

# only using subjects without missing SCR data
dataSCR_noNA <- dataSCR[!is.na(rowMeans(subset(dataSCR, select = Av_allTr:Min_2ndBl))),]
dataSCRLong_noNA <- dataSCRLong[!is.na(dataSCRLong$SCR),]

# frequentist ANOVA on SCR in classical conditioning group, including p. eta^2
# IV = CS; DV = SCR
anovaSCRReal <- ezANOVA(
  #data = dataSCRLong_noNA[dataSCRLong_noNA$usGroup == "real" & dataSCRLong_noNA$time == "allTr",],
  data = dataSCRLong[dataSCRLong$usGroup == "real" & dataSCRLong$time == "allTr",],
  dv = SCR,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaSCRReal$ANOVA$pEtaSq <- 
  c(anovaSCRReal$ANOVA$SSn[1] / (anovaSCRReal$ANOVA$SSd[1]+anovaSCRReal$ANOVA$SSn[1]),
    anovaSCRReal$ANOVA$SSn[2] / (anovaSCRReal$ANOVA$SSd[2]+anovaSCRReal$ANOVA$SSn[2])
  ); print(anovaSCRReal)

# bayesian ANOVA on SCR in classical conditioning group
set.seed(rngSeed); anovaBFSCRReal <- anovaBF(
  formula = SCR ~ CS + partInd,
  data = dataSCRLong_noNA[dataSCRLong_noNA$usGroup == "real" & dataSCRLong_noNA$time == "allTr",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFSCRReal)

# frequentist & bayesian t-tests on SCR in classical conditioning group
# CS+av vs CS+neu
scrRealAvNeu_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
scrRealAvNeu_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                           paired = TRUE)
scrRealAvNeu_BF <- ttestBF(x = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "real"],
                           y = dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrRealAvMin_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
scrRealAvMin_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           paired = TRUE)
scrRealAvMin_BF <- ttestBF(x = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "real"],
                           y = dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrRealNeuMin_t <- t.test(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                          y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                          alternative = "two.sided", paired = TRUE) # two-sided
scrRealNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            paired = TRUE)
scrRealNeuMin_BF <- ttestBF(x = dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "real"],
                            y = dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # two-sided



########################################
### Across groups - primary analyses ###
########################################

# descriptive statistics for SCR ratings across conditioning groups
describe(dataSCR)

# frequentist ANOVA on SCR across conditioning groups
anovaSCR <- ezANOVA(
  #data = dataSCRLong_noNA[dataSCRLong_noNA$time == "allTr",],
  data = dataSCRLong[dataSCRLong$time == "allTr",],
  dv = SCR,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaSCR$ANOVA$pEtaSq <- c(anovaSCR$ANOVA$SSn[1] /
                                (anovaSCR$ANOVA$SSd[1]+anovaSCR$ANOVA$SSn[1]),
                              anovaSCR$ANOVA$SSn[2] /
                                (anovaSCR$ANOVA$SSd[2]+anovaSCR$ANOVA$SSn[2]),
                              anovaSCR$ANOVA$SSn[3] /
                                (anovaSCR$ANOVA$SSd[3]+anovaSCR$ANOVA$SSn[3]),
                              anovaSCR$ANOVA$SSn[4] /
                                (anovaSCR$ANOVA$SSd[4]+anovaSCR$ANOVA$SSn[4])
); print(anovaSCR)

# bayesian ANOVA on SCR across conditioning groups
set.seed(rngSeed); anovaBFSCR <- anovaBF(
  formula = SCR ~ usGroup*CS + partInd,
  #data = dataSCRLong_noNA[dataSCRLong_noNA$time == "allTr",],
  data = dataSCRLong[dataSCRLong$time == "allTr",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFSCR)

# inclusion factors for bayesian ANOVA effects
bf_nullModel <- 1
bf_usGroup <- exp(anovaBFSCR@bayesFactor$bf[1])
bf_cs <- exp(anovaBFSCR@bayesFactor$bf[2])
bf_interact <- exp(anovaBFSCR@bayesFactor$bf[3])
bf_usGroup_cs <- exp(anovaBFSCR@bayesFactor$bf[4])
bf_usGroup_interact <- exp(anovaBFSCR@bayesFactor$bf[5])
bf_cs_interact <- exp(anovaBFSCR@bayesFactor$bf[6])
bf_fullModel <- exp(anovaBFSCR@bayesFactor$bf[7])

# main effect US group: models [1] and [3] vs. null model and model [2]
bfIncGroupSCR <- (bf_usGroup + bf_usGroup_cs + bf_usGroup_interact + bf_fullModel) / 
  (bf_nullModel + bf_cs + bf_interact + bf_cs_interact); bfIncGroupSCR
# main effect CS type: models "main effect CS" & "main effects CS & group" vs.
#                      null model and "main effect group"
bfIncCsSCR <- (bf_cs + bf_usGroup_cs + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_interact + bf_usGroup_interact); bfIncCsSCR
# interaction: Full model vs. main-effects-only model
bfIncInteractSCR <- (bf_interact + bf_usGroup_interact + bf_cs_interact + bf_fullModel) / 
  (bf_nullModel + bf_usGroup + bf_cs + bf_usGroup_cs); bfIncInteractSCR

# quick & dirty graph of group x CS ANOVA on valence ratings
ezPlot(
  data = dataSCRLong_noNA[dataSCRLong_noNA$time == "allTr",],
  dv = SCR,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on SCR (difference scores) across groups
# delta [CS+av - CS+neu]
scrBothAvNeu_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
scrBothAvNeu_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                               dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                               dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                           paired = FALSE)
scrBothAvNeu_BF <- ttestBF(x = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "real"] -
                               dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "real"],
                           y = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "ima"] -
                               dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
scrBothAvMin_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
scrBothAvMin_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                               dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                               dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                           paired = FALSE)
scrBothAvMin_BF <- ttestBF(x = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "real"] -
                               dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "real"],
                           y = dataSCR_noNA$Av_allTr[dataSCR_noNA$usGroup == "ima"] -
                               dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
scrBothNeuMin_t <- t.test(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                          y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
scrBothNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"] -
                                dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"] -
                                dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                            paired = FALSE)
scrBothNeuMin_BF <- ttestBF(x = dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "real"] -
                                dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "real"],
                            y = dataSCR_noNA$Neu_allTr[dataSCR_noNA$usGroup == "ima"] -
                                dataSCR_noNA$Min_allTr[dataSCR_noNA$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided



# quick & dirty graph of group x CS ANOVA on SCR
ezPlot(
  data = dataSCRLong_noNA,
  dv = SCR,
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
  t = c(scrImaAvNeu_t$statistic, scrImaAvMin_t$statistic, scrImaNeuMin_t$statistic,
        scrRealAvNeu_t$statistic, scrRealAvMin_t$statistic, scrRealNeuMin_t$statistic,
        scrBothAvNeu_t$statistic, scrBothAvMin_t$statistic, scrBothNeuMin_t$statistic),
  df = c(scrImaAvNeu_t$parameter, scrImaAvMin_t$parameter, scrImaNeuMin_t$parameter,
         scrRealAvNeu_t$parameter, scrRealAvMin_t$parameter, scrRealNeuMin_t$parameter,
         scrBothAvNeu_t$parameter, scrBothAvMin_t$parameter, scrBothNeuMin_t$parameter), 
  p = c(scrImaAvNeu_t$p.value, scrImaAvMin_t$p.value, scrImaNeuMin_t$p.value,
        scrRealAvNeu_t$p.value, scrRealAvMin_t$p.value, scrRealNeuMin_t$p.value,
        scrBothAvNeu_t$p.value*3, scrBothAvMin_t$p.value*3, scrBothNeuMin_t$p.value*3),  # Bonferroni
  d = c(scrImaAvNeu_d$Cohens_d, scrImaAvMin_d$Cohens_d, scrImaNeuMin_d$Cohens_d,
        scrRealAvNeu_d$Cohens_d, scrRealAvMin_d$Cohens_d, scrRealNeuMin_d$Cohens_d,
        scrBothAvNeu_d$Cohens_d, scrBothAvMin_d$Cohens_d, scrBothNeuMin_d$Cohens_d),
  BF = c(exp(scrImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrImaAvMin_BF@bayesFactor[["bf"]][1]), exp(scrImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(scrRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrRealAvMin_BF@bayesFactor[["bf"]][1]), exp(scrRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(scrBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrBothAvMin_BF@bayesFactor[["bf"]][1]), exp(scrBothNeuMin_BF@bayesFactor[["bf"]][1]))
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

tableSCR <- flextable(tableData[1:9,])
tableSCR <- add_header_lines(tableSCR, top = TRUE, values = "SCR")
tableSCR <- align(tableSCR, align = "center")

save_as_docx(tableSCR, path = "Tables/tableSCR_raw.docx")



####################
### Plotting SCR ###
####################

csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-",
             expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")

dataSCRWithin <- dataSCR_noNA[,c("partInd","usGroup","Av_allTr","Neu_allTr","Min_allTr")]
# remove each participant's average from each single value
dataSCRWithin[,3:5] <- as.matrix(dataSCRWithin[,3:5]) -
  rowMeans(as.matrix(dataSCRWithin[,3:5])) 
# prepare data frame for bar plot with means from standard dataset and SE from
# dataset without betweem-subject variance
meanSCR <- data.frame(
  usGroup = factor(c(rep(1,3),rep(2,3)),
                   labels = c("Imagery-Based","Classical")),
  CS = factor(c(1,2,3,1,2,3),
              labels = csLabels[1:3]),
  mean = c(describe(dataSCR[dataSCR$usGroup == "ima", c(3,6,9)])$mean,
           describe(dataSCR[dataSCR$usGroup == "real", c(3,6,9)])$mean),
  se = c(describe(dataSCRWithin[dataSCRWithin$usGroup == "ima", 3:5])$se,
         describe(dataSCRWithin[dataSCRWithin$usGroup == "real", 3:5])$se)
)



plotFS <- 9
showSig <- TRUE

graphSCR <- ggplot(data = meanSCR, aes(x = usGroup, y = mean, fill = CS)) +
  theme_classic() +
  geom_col(aes(fill = CS), position = position_dodge(width = .9)) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = .1), position = position_dodge(width = .9)) +
  scale_x_discrete(aes(breaks = usGroup), name = "") +
  scale_y_continuous(name = "Mean SCR (normalized)") +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = usGroup, y = .2), colour = "black", size = (plotFS/.pt)-.5, fontface = "bold") +
  geom_text(aes(y = -0.01), label = csLabels, position = position_dodge(width = .9), colour = "black", size = (plotFS/.pt)-.5, fontface = "bold") + 
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"),
        plot.margin = unit(c(5,5,5,5), "mm"))

if (showSig == TRUE){
  graphSCR <- graphSCR +
    geom_segment(aes(x = 1.7, y = mean+se+.005, xend = 2.0, yend = mean+se+.005), data = meanSCR[4,]) +
    geom_text(aes(label = "***", x = 1.85, y = mean+se+.01), size = plotFS/2, data = meanSCR[4,]) +
    geom_segment(aes(x = 1.7, y = mean+se+.02, xend = 2.3, yend = mean+se+.02), data = meanSCR[4,]) +
    geom_text(aes(label = "***", x = 2.0, y = mean+se+.025), size = plotFS/2, data = meanSCR[4,])
}
graphSCR

# saving it
ggsave(filename = "Figures/Figure3_barPlot_SCR.eps",
       plot = graphSCR,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)

ggsave(filename = "Figures/Figure3_barPlot_SCR.pdf",
       plot = graphSCR,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)

