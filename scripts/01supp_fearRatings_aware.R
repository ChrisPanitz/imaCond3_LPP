# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: Supplementary analyses on fear ratings: repeat main analyses in contingency-aware participants

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
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# extract contingency aware participants
importRatings <- importRatings[importRatings$contAware == TRUE,]

# create data frames in wide & long format for fear ratings
dataFear <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resolve issue in which Bayes ANOVA crashes using alphanumeric codes as participant ID
  partCode = factor(importRatings$partCode),
  usGroup = factor(importRatings$group, levels = c("ima", "real")),
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



#######################################################
### Across groups - fear ratings - primary analyses ###
#######################################################

# descriptive statistics  for fear ratings across conditioning groups
describe(dataFear)

# frequentist ANOVA on fear ratings across conditioning groups
anovaFear <- ezANOVA(
  data = dataFearLong[dataFearLong$time == "Post",],
  dv = fear,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  type = 3,
  detailed = TRUE
); anovaFear$ANOVA$pEtaSq <- c(anovaFear$ANOVA$SSn[1] /
                                 (anovaFear$ANOVA$SSd[1]+anovaFear$ANOVA$SSn[1]),
                               anovaFear$ANOVA$SSn[2] /
                                 (anovaFear$ANOVA$SSd[2]+anovaFear$ANOVA$SSn[2]),
                               anovaFear$ANOVA$SSn[3] /
                                 (anovaFear$ANOVA$SSd[3]+anovaFear$ANOVA$SSn[3]),
                               anovaFear$ANOVA$SSn[4] /
                                 (anovaFear$ANOVA$SSd[4]+anovaFear$ANOVA$SSn[4])
); print(anovaFear)

# bayesian ANOVA on fear ratings across conditioning groups
set.seed(rngSeed); anovaBFFear <- anovaBF(
  formula = fear ~ usGroup*CS + partInd,
  data = dataFearLong[dataFearLong$time == "Post",],
  whichRandom = "partInd",
  whichModels = "all",
  iterations = 100000
); print(anovaBFFear)

# inclusion factors for bayesian ANOVA effects
bf_inclusion(anovaBFFear)

# quick graph of group x CS ANOVA on fear ratings
ezPlot(
  data = dataFearLong[dataFearLong$time == "Post",],
  dv = fear,
  wid = partInd,
  within = .(CS),
  between = .(usGroup),
  x = CS,
  split = usGroup
)  

# frequentist & bayesian t-tests on fear ratings across groups
# CS+av vs CS+neu
fearAcrossAvNeu_t <- t.test(x = dataFear$Av_Post,
                            y = dataFear$Neu_Post,
                            alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvNeu_d <- cohens_d(x = dataFear$Av_Post,
                              y = dataFear$Neu_Post,
                              paired = TRUE)
fearAcrossAvNeu_BF <- ttestBF(x = dataFear$Av_Post,
                              y = dataFear$Neu_Post,
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearAcrossAvMin_t <- t.test(x = dataFear$Av_Post,
                            y = dataFear$Min_Post,
                            alternative = "greater", paired = TRUE) # one-sided
fearAcrossAvMin_d <- cohens_d(x = dataFear$Av_Post,
                              y = dataFear$Min_Post,
                              paired = TRUE)
fearAcrossAvMin_BF <- ttestBF(x = dataFear$Av_Post,
                              y = dataFear$Min_Post,
                              nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearAcrossNeuMin_t <- t.test(x = dataFear$Neu_Post,
                             y = dataFear$Min_Post,
                             alternative = "two.sided", paired = TRUE) # two-sided
fearAcrossNeuMin_d <- cohens_d(x = dataFear$Neu_Post,
                               y = dataFear$Min_Post,
                               paired = TRUE)
fearAcrossNeuMin_BF <- ttestBF(x = dataFear$Neu_Post,
                               y = dataFear$Min_Post,
                               nullIntervall = NULL, paired = TRUE) # two-sided

# frequentist & bayesian t-tests on fear ratings (difference scores) across groups
# delta [CS+av - CS+neu]
fearBetweenAvNeu_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                               dataFear$Neu_Post[dataFear$usGroup == "real"],
                             y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                               dataFear$Neu_Post[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvNeu_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Post[dataFear$usGroup == "real"],
                               y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Post[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvNeu_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                 dataFear$Neu_Post[dataFear$usGroup == "real"],
                               y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                 dataFear$Neu_Post[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
fearBetweenAvMin_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                               dataFear$Min_Post[dataFear$usGroup == "real"],
                             y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                               dataFear$Min_Post[dataFear$usGroup == "ima"],
                             alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenAvMin_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                 dataFear$Min_Post[dataFear$usGroup == "real"],
                               y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Post[dataFear$usGroup == "ima"],
                               paired = FALSE)
fearBetweenAvMin_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                                 dataFear$Min_Post[dataFear$usGroup == "real"],
                               y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                                 dataFear$Min_Post[dataFear$usGroup == "ima"],
                               nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
fearBetweenNeuMin_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                dataFear$Min_Post[dataFear$usGroup == "real"],
                              y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                                dataFear$Min_Post[dataFear$usGroup == "ima"],
                              alternative = "two.sided", paired = FALSE) # two-sided
fearBetweenNeuMin_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                  dataFear$Min_Post[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Post[dataFear$usGroup == "ima"] -
                                  dataFear$Min_Post[dataFear$usGroup == "ima"],
                                paired = FALSE)
fearBetweenNeuMin_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                                  dataFear$Min_Post[dataFear$usGroup == "real"],
                                y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                                  dataFear$Min_Post[dataFear$usGroup == "ima"],
                                nullInterval = NULL, paired = FALSE) # two-sided



####################################################################
### Imagery-based conditioning - fear ratings - primary analyses ###
####################################################################

# descriptive statistics for fear ratings in imagery-based conditioning group
describe(dataFear[dataFear$usGroup == "ima",])

# frequentist ANOVA in imagery-based conditioning group, including p. eta^2
# IV = CS; DV = fear rating
anovaFearIma <- ezANOVA(
  data = dataFearLong[dataFearLong$usGroup == "ima" & dataFearLong$time == "Post",],
  dv = fear,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaFearIma$ANOVA$pEtaSq <-
  c(anovaFearIma$ANOVA$SSn[1] / (anovaFearIma$ANOVA$SSd[1]+anovaFearIma$ANOVA$SSn[1]),
    anovaFearIma$ANOVA$SSn[2] / (anovaFearIma$ANOVA$SSd[2]+anovaFearIma$ANOVA$SSn[2])
  ); print(anovaFearIma)

# bayesian ANOVA on fear ratings in imagery-based conditioning group
set.seed(rngSeed); anovaBFFearIma <- anovaBF(
  formula = fear ~ CS + partInd,
  data = dataFearLong[dataFearLong$usGroup == "ima" & dataFearLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFFearIma)

# frequentist & bayesian t-tests on fear ratings in imagery-based conditioning group
# CS+av vs CS+neu
fearImaAvNeu_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                         y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                         alternative = "greater", paired = TRUE) # one-sided
fearImaAvNeu_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                           y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                           paired = TRUE)
fearImaAvNeu_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                           y = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearImaAvMin_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                         y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                         alternative = "greater", paired = TRUE) # one-sided
fearImaAvMin_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                           y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                           paired = TRUE)
fearImaAvMin_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "ima"],
                           y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                           nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearImaNeuMin_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                          y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                          alternative = "two.sided", paired = TRUE) # two-sided
fearImaNeuMin_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                            y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                            paired = TRUE)
fearImaNeuMin_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "ima"],
                            y = dataFear$Min_Post[dataFear$usGroup == "ima"],
                            nullIntervall = NULL, paired = TRUE) # two-sided



################################################################
### Classical conditioning - fear ratings - primary analyses ###
################################################################

# descriptive statistics for fear ratings in classical conditioning group
describe(dataFear[dataFear$usGroup == "real",])

# frequentist ANOVA in classical conditioning group, including p. eta^2
# IV = CS; DV = fear rating
anovaFearReal <- ezANOVA(
  data = dataFearLong[dataFearLong$usGroup == "real" & dataFearLong$time == "Post",],
  dv = fear,
  wid = partInd,
  within = .(CS),
  type = 3,
  detailed = TRUE
); anovaFearReal$ANOVA$pEtaSq <-
  c(anovaFearReal$ANOVA$SSn[1] / (anovaFearReal$ANOVA$SSd[1]+anovaFearReal$ANOVA$SSn[1]),
    anovaFearReal$ANOVA$SSn[2] / (anovaFearReal$ANOVA$SSd[2]+anovaFearReal$ANOVA$SSn[2])
  ); print(anovaFearReal)

# bayesian ANOVA on fear ratings in classical conditioning group
set.seed(rngSeed); anovaBFFearReal <- anovaBF(
  formula = fear ~ CS + partInd,
  data = dataFearLong[dataFearLong$usGroup == "real" & dataFearLong$time == "Post",],
  whichRandom = "partInd",
  iterations = 100000
); print(anovaBFFearReal)

# frequentist & bayesian t-tests on fear ratings in classical conditioning group
# CS+av vs CS+neu
fearRealAvNeu_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                          y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                          alternative = "greater", paired = TRUE) # one-sided
fearRealAvNeu_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                            y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                            paired = TRUE)
fearRealAvNeu_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                            y = dataFear$Neu_Post[dataFear$usGroup == "real"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+av vs CS-
fearRealAvMin_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                          y = dataFear$Min_Post[dataFear$usGroup == "real"],
                          alternative = "greater", paired = TRUE) # one-sided
fearRealAvMin_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                            y = dataFear$Min_Post[dataFear$usGroup == "real"],
                            paired = TRUE)
fearRealAvMin_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"],
                            y = dataFear$Min_Post[dataFear$usGroup == "real"],
                            nullInterval = c(0, Inf), paired = TRUE) # one-sided x > y
# CS+neu vs CS-
fearRealNeuMin_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                           y = dataFear$Min_Post[dataFear$usGroup == "real"],
                           alternative = "two.sided", paired = TRUE) # two-sided
fearRealNeuMin_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                             y = dataFear$Min_Post[dataFear$usGroup == "real"],
                             paired = TRUE)
fearRealNeuMin_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "real"],
                             y = dataFear$Min_Post[dataFear$usGroup == "real"],
                             nullIntervall = NULL, paired = TRUE) # two-sided



#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("total sample: CS+av vs CS+neu", "total sample: CS+av vs CS-", "total sample: CSneu vs CS-",
                     "between groups: delta CS+av / CS+neu", "between groups: delta CS+av / CS-", "between groups: delta CSneu / CS-",
                     "imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-"), 5),
  t = c(fearAcrossAvNeu_t$statistic, fearAcrossAvMin_t$statistic, fearAcrossNeuMin_t$statistic,
        fearImaAvNeu_t$statistic, fearImaAvMin_t$statistic, fearImaNeuMin_t$statistic,
        fearRealAvNeu_t$statistic, fearRealAvMin_t$statistic, fearRealNeuMin_t$statistic, 
        fearBetweenAvNeu_t$statistic, fearBetweenAvMin_t$statistic, fearBetweenNeuMin_t$statistic),
  df = c(fearAcrossAvNeu_t$parameter, fearAcrossAvMin_t$parameter, fearAcrossNeuMin_t$parameter,
         fearImaAvNeu_t$parameter, fearImaAvMin_t$parameter, fearImaNeuMin_t$parameter,
         fearRealAvNeu_t$parameter, fearRealAvMin_t$parameter, fearRealNeuMin_t$parameter, 
         fearBetweenAvNeu_t$parameter, fearBetweenAvMin_t$parameter, fearBetweenNeuMin_t$parameter),
  p = c(fearAcrossAvNeu_t$p.value, fearAcrossAvMin_t$p.value, fearAcrossNeuMin_t$p.value,
        fearImaAvNeu_t$p.value, fearImaAvMin_t$p.value, fearImaNeuMin_t$p.value,
        fearRealAvNeu_t$p.value, fearRealAvMin_t$p.value, fearRealNeuMin_t$p.value,
        fearBetweenAvNeu_t$p.value*3, fearBetweenAvMin_t$p.value*3, fearBetweenNeuMin_t$p.value*3), # Bonferroni
  d = c(fearAcrossAvNeu_d$Cohens_d, fearAcrossAvMin_d$Cohens_d, fearAcrossNeuMin_d$Cohens_d,
        fearImaAvNeu_d$Cohens_d, fearImaAvMin_d$Cohens_d, fearImaNeuMin_d$Cohens_d,
        fearRealAvNeu_d$Cohens_d, fearRealAvMin_d$Cohens_d, fearRealNeuMin_d$Cohens_d, 
        fearBetweenAvNeu_d$Cohens_d, fearBetweenAvMin_d$Cohens_d, fearBetweenNeuMin_d$Cohens_d),
  BF = c(exp(fearAcrossAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearAcrossAvMin_BF@bayesFactor[["bf"]][1]), exp(fearAcrossNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(fearImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearImaAvMin_BF@bayesFactor[["bf"]][1]), exp(fearImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(fearRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearRealAvMin_BF@bayesFactor[["bf"]][1]), exp(fearRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(fearBetweenAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearBetweenAvMin_BF@bayesFactor[["bf"]][1]), exp(fearBetweenNeuMin_BF@bayesFactor[["bf"]][1]))
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

tableFear <- flextable(tableData[1:12,])
tableFear <- add_header_lines(tableFear, top = TRUE, values = "fear")

save_as_docx(tableFear, path = paste0(pathname, "/supplement/01s_tableFear_aware_raw.docx"))



#############################
### Plot for fear ratings ###
#############################
# remove between-subject variance for plotting standard errors based on
# within-subject variance
dataFearLongPost <- dataFear[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]
# remove each participant's average from each single value
dataFearLongPost[dataFearLongPost$usGroup == "ima",6:8] <- as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "ima",3:5]) -
  rowMeans(as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "ima",3:5])) + mean(as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "ima",3:5]))
dataFearLongPost[dataFearLongPost$usGroup == "real",6:8] <- as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "real",3:5]) -
  rowMeans(as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "real",3:5])) + mean(as.matrix(dataFearLongPost[dataFearLongPost$usGroup == "real",3:5]))
names(dataFearLongPost) <- c("partInd","usGroup","Av_btw","Neu_btw","Min_btw","Av_wth","Neu_wth","Min_wth")
# into long format
dataFearLongPost <- pivot_longer(data = dataFearLongPost, cols = Av_btw:Min_wth,
                                 names_to = c("CS","variance"), names_sep = "_", values_to = "fear")
dataFearLongPost <- pivot_wider(data = dataFearLongPost, names_from = "variance", values_from = "fear")
dataFearLongPost$CS <- factor(dataFearLongPost$CS, levels = c("Av","Neu","Min"))
levels(dataFearLongPost$usGroup) <- c("Imagery-based","Classical")

# some general settings
plotFS <- 8
showSig <- TRUE
csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")

# plotting
graphFear <- ggplot(data = dataFearLongPost, aes(x = usGroup, y = btw, fill = CS, color = CS)) +
  theme_classic() +
  stat_summary(aes(y = wth), fun.data = mean_se, geom = "errorbar", position=position_dodge(0.8), width = 0.1, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.8), width = 0.25, linewidth = 0.2) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_color_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  scale_y_continuous(name = "Fear rating (1-5)", limits = c(0.4,6), breaks = 1:5, oob = rescale_none, expand = c(0,0)) +
  geom_vline(xintercept = 0.41) +
  geom_rect(xmin = 0.4, xmax = 2.6, ymin = 0.35, ymax = 0.8, fill = "white", color = "white") +
  geom_hline(yintercept = 0.8) +
  geom_beeswarm(aes(color = CS), dodge.width = 0.8, cex = 0.6, size = .1, color = "gray70") +
  geom_violin(alpha = .2, color = NA, bw = .5, position = position_dodge(0.8), width = 0.75) +
  geom_text(x = 0.7, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.3, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
  geom_text(x = 1.7, y = 0.6, label = csLabels[1], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.0, y = 0.6, label = csLabels[2], colour = "black", size = plotFS/.pt) +
  geom_text(x = 2.3, y = 0.6, label = csLabels[3], colour = "black", size = plotFS/.pt) +
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
  graphFear <- graphFear +
    geom_segment(aes(x = 0.74, y = 5.1, xend = 1.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "**", x = 0.85, y =5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 0.74, y = 5.4, xend = 1.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.0, y = 5.45), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.1, xend = 2.0, yend = 5.1), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 1.85, y = 5.15), size = plotFS/4, color = "gray20") +
    geom_segment(aes(x = 1.74, y = 5.4, xend = 2.26, yend = 5.4), linewidth = 0.2, color = "gray20") +
    geom_text(aes(label = "***", x = 2.0, y = 5.45), size = plotFS/4, color = "gray20")
}
graphFear
