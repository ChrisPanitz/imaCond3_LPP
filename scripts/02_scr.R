# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: Main analyses on SCR (ANOVAs, pairwise comparisons, plotting)

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

# Rating data contains group membership
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# load SCR data from text file
importSCR <- read.csv(paste0(pathname, "/experimentData/imaCond3_scrmatrix_cs.txt"), sep = "")
importSCR <- merge(importRatings[,c("partCode","group")], importSCR, by = "partCode")
importSCR <- importSCR[order(importSCR$partCode),]

# create data frames in wide & long format for SCR
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


########################################
### Across groups - primary analyses ###
########################################

# descriptive statistics for SCR ratings across conditioning groups
describe(dataSCR)

# frequentist ANOVA on SCR across conditioning groups
anovaSCR <- ezANOVA(
  #data = dataSCRLong[dataSCRLong$time == "allTr",],
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
  #data = dataSCRLong[dataSCRLong$time == "allTr",],
  data = dataSCRLong[dataSCRLong$time == "allTr",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFSCR)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFSCR)

# quick graph of group x CS ANOVA on SCR
ezPlot(
  data = dataSCRLong[dataSCRLong$time == "allTr",],
  dv = SCR,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on SCR in imagery-based conditioning group
# CS+av vs CS+neu
scrAcrossAvNeu_t <- t.test(x = dataSCR$Av_allTr,
                        y = dataSCR$Neu_allTr,
                        alternative = "greater", paired = TRUE) # one-sided
scrAcrossAvNeu_d <- cohens_d(x = dataSCR$Av_allTr,
                          y = dataSCR$Neu_allTr,
                          paired = TRUE)
scrAcrossAvNeu_BF <- ttestBF(x = dataSCR$Av_allTr,
                          y = dataSCR$Neu_allTr,
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrAcrossAvMin_t <- t.test(x = dataSCR$Av_allTr,
                        y = dataSCR$Min_allTr,
                        alternative = "greater", paired = TRUE) # one-sided
scrAcrossAvMin_d <- cohens_d(x = dataSCR$Av_allTr,
                          y = dataSCR$Min_allTr,
                          paired = TRUE)
scrAcrossAvMin_BF <- ttestBF(x = dataSCR$Av_allTr,
                          y = dataSCR$Min_allTr,
                          nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrAcrossNeuMin_t <- t.test(x = dataSCR$Neu_allTr,
                         y = dataSCR$Min_allTr,
                         alternative = "two.sided", paired = TRUE) # two-sided
scrAcrossNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr,
                           y = dataSCR$Min_allTr,
                           paired = TRUE)
scrAcrossNeuMin_BF <- ttestBF(x = dataSCR$Neu_allTr,
                           y = dataSCR$Min_allTr,
                           nullInterval = NULL, paired = TRUE) # two-sided


# frequentist & bayesian t-tests on SCR (difference scores) between groups
# delta [CS+av - CS+neu]
scrBetweenAvNeu_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                           dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                           dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
scrBetweenAvNeu_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                           paired = FALSE)
scrBetweenAvNeu_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Neu_allTr[dataSCR$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
scrBetweenAvMin_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                           dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                           dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                         alternative = "two.sided", paired = FALSE) # two-sided
scrBetweenAvMin_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                           paired = FALSE)
scrBetweenAvMin_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Av_allTr[dataSCR$usGroup == "ima"] -
                             dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                           nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
scrBetweenNeuMin_t <- t.test(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"] -
                            dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                          y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"] -
                            dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
scrBetweenNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                            paired = FALSE)
scrBetweenNeuMin_BF <- ttestBF(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Neu_allTr[dataSCR$usGroup == "ima"] -
                              dataSCR$Min_allTr[dataSCR$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided



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
dataSCR <- dataSCR[!is.na(rowMeans(subset(dataSCR, select = Av_allTr:Min_2ndBl))),]
dataSCRLong <- dataSCRLong[!is.na(dataSCRLong$SCR),]

# frequentist ANOVA on SCR in classical conditioning group, including p. eta^2
# IV = CS; DV = SCR
anovaSCRReal <- ezANOVA(
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
  data = dataSCRLong[dataSCRLong$usGroup == "real" & dataSCRLong$time == "allTr",],
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
scrRealAvNeu_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
scrRealAvMin_t <- t.test(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                         y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                         alternative = "greater", paired = TRUE) # one-sided
scrRealAvMin_d <- cohens_d(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           paired = TRUE)
scrRealAvMin_BF <- ttestBF(x = dataSCR$Av_allTr[dataSCR$usGroup == "real"],
                           y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
scrRealNeuMin_t <- t.test(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                          y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                          alternative = "two.sided", paired = TRUE) # two-sided
scrRealNeuMin_d <- cohens_d(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            paired = TRUE)
scrRealNeuMin_BF <- ttestBF(x = dataSCR$Neu_allTr[dataSCR$usGroup == "real"],
                            y = dataSCR$Min_allTr[dataSCR$usGroup == "real"],
                            nullInterval = NULL, paired = TRUE) # two-sided



#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("across groups: CS+av vs CS+neu", "across groups: CS+av vs CS-", "across groups: CSneu vs CS-",
                     "imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-",
                     "between groups: delta CS+av / CS+neu", "between groups: delta CS+av / CS-", "between groups: delta CSneu / CS-"), 3),
  t = c(scrAcrossAvNeu_t$statistic, scrAcrossAvMin_t$statistic, scrAcrossNeuMin_t$statistic,
        scrImaAvNeu_t$statistic, scrImaAvMin_t$statistic, scrImaNeuMin_t$statistic,
        scrRealAvNeu_t$statistic, scrRealAvMin_t$statistic, scrRealNeuMin_t$statistic,
        scrBetweenAvNeu_t$statistic, scrBetweenAvMin_t$statistic, scrBetweenNeuMin_t$statistic),
  df = c(scrAcrossAvNeu_t$parameter, scrAcrossAvMin_t$parameter, scrAcrossNeuMin_t$parameter,
         scrImaAvNeu_t$parameter, scrImaAvMin_t$parameter, scrImaNeuMin_t$parameter,
         scrRealAvNeu_t$parameter, scrRealAvMin_t$parameter, scrRealNeuMin_t$parameter,
         scrBetweenAvNeu_t$parameter, scrBetweenAvMin_t$parameter, scrBetweenNeuMin_t$parameter), 
  p = c(scrAcrossAvNeu_t$p.value, scrAcrossAvMin_t$p.value, scrAcrossNeuMin_t$p.value,
        scrImaAvNeu_t$p.value, scrImaAvMin_t$p.value, scrImaNeuMin_t$p.value,
        scrRealAvNeu_t$p.value, scrRealAvMin_t$p.value, scrRealNeuMin_t$p.value,
        scrBetweenAvNeu_t$p.value*3, scrBetweenAvMin_t$p.value*3, scrBetweenNeuMin_t$p.value*3),  # Bonferroni
  d = c(scrAcrossAvNeu_d$Cohens_d, scrAcrossAvMin_d$Cohens_d, scrAcrossNeuMin_d$Cohens_d,
        scrImaAvNeu_d$Cohens_d, scrImaAvMin_d$Cohens_d, scrImaNeuMin_d$Cohens_d,
        scrRealAvNeu_d$Cohens_d, scrRealAvMin_d$Cohens_d, scrRealNeuMin_d$Cohens_d,
        scrBetweenAvNeu_d$Cohens_d, scrBetweenAvMin_d$Cohens_d, scrBetweenNeuMin_d$Cohens_d),
  BF = c(exp(scrAcrossAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrAcrossAvMin_BF@bayesFactor[["bf"]][1]), exp(scrAcrossNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(scrImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrImaAvMin_BF@bayesFactor[["bf"]][1]), exp(scrImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(scrRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrRealAvMin_BF@bayesFactor[["bf"]][1]), exp(scrRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(scrBetweenAvNeu_BF@bayesFactor[["bf"]][1]), exp(scrBetweenAvMin_BF@bayesFactor[["bf"]][1]), exp(scrBetweenNeuMin_BF@bayesFactor[["bf"]][1]))
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
tableData$BF <- format(tableData$BF, digits = 3)

tableSCR <- flextable(tableData[1:12,])
tableSCR <- add_header_lines(tableSCR, top = TRUE, values = "SCR")
tableSCR <- align(tableSCR, align = "center")

save_as_docx(tableSCR, path = paste0(pathname, "/tables/tableSCR_raw.docx"))



####################
### Plotting SCR ###
####################
# remove between-subject variance for plotting standard errors based on
# within-subject variance
dataSCRLongAll <- dataSCR[,c("partInd","usGroup","Av_allTr","Neu_allTr","Min_allTr")]
# remove each participant's average from each single value
dataSCRLongAll[dataSCRLongAll$usGroup == "ima",6:8] <- as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "ima",3:5])) + mean(as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "ima",3:5]))
dataSCRLongAll[dataSCRLongAll$usGroup == "real",6:8] <- as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "real",3:5])) + mean(as.matrix(dataSCRLongAll[dataSCRLongAll$usGroup == "real",3:5]))
names(dataSCRLongAll) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")
# into long format
dataSCRLongAll <- pivot_longer(data = dataSCRLongAll, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "scr")
dataSCRLongAll <- pivot_wider(data = dataSCRLongAll, names_from = "variance", values_from = "scr")
dataSCRLongAll$CS <- factor(dataSCRLongAll$CS, levels = c("Av","Neu","Min"))
levels(dataSCRLongAll$usGroup) <- c("Imagery-based","Classical")

# some general settings
plotFS <- 8
showSig <- TRUE
csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")

# plotting
graphSCR <- ggplot(data = dataSCRLongAll, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  geom_violin(alpha = .2, color = NA, bw = .1, position = position_dodge(0.8), width = 0.75) +
  geom_quasirandom(dodge.width = 0.8, size = .2, color = "gray70", width = 0.03) +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.2, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.35, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Mean SCR (normalized)", limits = c(-0.25, 0.6)) +
  geom_hline(yintercept = -0.2) +
  geom_text(x = 0.74, y = -0.25, label = csLabels[1], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(x = 1.0, y = -0.25, label = csLabels[2], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(x = 1.26, y = -0.25, label = csLabels[3], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(x = 1.74, y = -0.25, label = csLabels[1], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(x = 2.0, y = -0.25, label = csLabels[2], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(x = 2.26, y = -0.25, label = csLabels[3], colour = "black", size = (plotFS-1)/.pt) +
  geom_text(aes(label = usGroup, y = .6), colour = "black", size = ((plotFS-2)/.pt), fontface = "bold") +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),        axis.title.y = element_text(margin = margin(r = 5), size = plotFS),
        axis.text.y = element_text(margin = margin(r = 5), size = plotFS, color = "black"),
        axis.ticks.y = element_line(colour = "black"),
        plot.margin = unit(c(5,5,5,5), "mm"))

if (showSig == TRUE){
  graphSCR <- graphSCR +
    geom_segment(aes(x = 1.74, y = 0.46, xend = 2.0, yend = 0.46), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.87, y = 0.465), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 0.50, xend = 2.26, yend = 0.50), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 2.0, y = 0.505), size = plotFS/4, color = "gray20")
}
graphSCR

# saving it
ggsave(filename = "Figures/Figure3_plotSCR.eps",
       plot = graphSCR,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)

ggsave(filename = "Figures/Figure3_plotSCR.pdf",
       plot = graphSCR,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)