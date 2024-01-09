# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.2.2 (2022-10-31) -- "Innocent and Trusting"
# --- RStudio version: 2022.12.0
# --- script version: Mar 2022
# --- content: rating analyses of imagery-based conditioning data in Panitz & Mueller (2023)

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

# load rating data from text file
# (see imaCond3_allratings_readme.txt for more details)
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# create data frames in wide & long format for fear ratings
dataFear <- data.frame(
  partInd = factor(1:dim(importRatings)[1]), # could not resovle issue in which Bayes ANOVA crashes using alphanumeric codes as prticipant ID
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

# frequentist & bayesian t-tests on fear ratings (difference scores) across groups
# delta [CS+av - CS+neu]
fearBothAvNeu_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                            dataFear$Neu_Post[dataFear$usGroup == "real"],
                          y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                            dataFear$Neu_Post[dataFear$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
fearBothAvNeu_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                              dataFear$Neu_Post[dataFear$usGroup == "real"],
                            y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                              dataFear$Neu_Post[dataFear$usGroup == "ima"],
                            paired = FALSE)
fearBothAvNeu_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                              dataFear$Neu_Post[dataFear$usGroup == "real"],
                            y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                              dataFear$Neu_Post[dataFear$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+av - CS-]
fearBothAvMin_t <- t.test(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                            dataFear$Min_Post[dataFear$usGroup == "real"],
                          y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                            dataFear$Min_Post[dataFear$usGroup == "ima"],
                          alternative = "two.sided", paired = FALSE) # two-sided
fearBothAvMin_d <- cohens_d(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                              dataFear$Min_Post[dataFear$usGroup == "real"],
                            y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                              dataFear$Min_Post[dataFear$usGroup == "ima"],
                            paired = FALSE)
fearBothAvMin_BF <- ttestBF(x = dataFear$Av_Post[dataFear$usGroup == "real"] -
                              dataFear$Min_Post[dataFear$usGroup == "real"],
                            y = dataFear$Av_Post[dataFear$usGroup == "ima"] -
                              dataFear$Min_Post[dataFear$usGroup == "ima"],
                            nullInterval = NULL, paired = FALSE) # two-sided
# delta [CS+neu - CS-]
fearBothNeuMin_t <- t.test(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                             dataFear$Min_Post[dataFear$usGroup == "real"],
                           y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                             dataFear$Min_Post[dataFear$usGroup == "ima"],
                           alternative = "two.sided", paired = FALSE) # two-sided
fearBothNeuMin_d <- cohens_d(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                               dataFear$Min_Post[dataFear$usGroup == "real"],
                             y = dataFear$Neu_Post[dataFear$usGroup == "ima"] -
                               dataFear$Min_Post[dataFear$usGroup == "ima"],
                             paired = FALSE)
fearBothNeuMin_BF <- ttestBF(x = dataFear$Neu_Post[dataFear$usGroup == "real"] -
                               dataFear$Min_Post[dataFear$usGroup == "real"],
                             y = dataFear$Neu_Post[dataFear$usGroup == "ima"] - 
                               dataFear$Min_Post[dataFear$usGroup == "ima"],
                             nullInterval = NULL, paired = FALSE) # two-sided



#########################
### Table for t-tests ###
#########################

tableData <- data.frame(
  comparison = rep(c("imagery: CS+av vs CS+neu", "imagery: CS+av vs CS-", "imagery: CSneu vs CS-",
                     "classical: CS+av vs CS+neu", "classical: CS+av vs CS-", "classical: CSneu vs CS-",
                     "groups: delta CS+av / CS+neu", "groups: delta CS+av / CS-", "groups: delta CSneu / CS-"), 5),
  t = c(fearImaAvNeu_t$statistic, fearImaAvMin_t$statistic, fearImaNeuMin_t$statistic,
        fearRealAvNeu_t$statistic, fearRealAvMin_t$statistic, fearRealNeuMin_t$statistic,
        fearBothAvNeu_t$statistic, fearBothAvMin_t$statistic, fearBothNeuMin_t$statistic), 
  df = c(fearImaAvNeu_t$parameter, fearImaAvMin_t$parameter, fearImaNeuMin_t$parameter,
         fearRealAvNeu_t$parameter, fearRealAvMin_t$parameter, fearRealNeuMin_t$parameter,
         fearBothAvNeu_t$parameter, fearBothAvMin_t$parameter, fearBothNeuMin_t$parameter), 
  p = c(fearImaAvNeu_t$p.value, fearImaAvMin_t$p.value, fearImaNeuMin_t$p.value,
        fearRealAvNeu_t$p.value, fearRealAvMin_t$p.value, fearRealNeuMin_t$p.value,
        fearBothAvNeu_t$p.value*3, fearBothAvMin_t$p.value*3, fearBothNeuMin_t$p.value*3),  # Bonferroni
  d = c(fearImaAvNeu_d$Cohens_d, fearImaAvMin_d$Cohens_d, fearImaNeuMin_d$Cohens_d,
        fearRealAvNeu_d$Cohens_d, fearRealAvMin_d$Cohens_d, fearRealNeuMin_d$Cohens_d,
        fearBothAvNeu_d$Cohens_d, fearBothAvMin_d$Cohens_d, fearBothNeuMin_d$Cohens_d), 
  BF = c(exp(fearImaAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearImaAvMin_BF@bayesFactor[["bf"]][1]), exp(fearImaNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(fearRealAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearRealAvMin_BF@bayesFactor[["bf"]][1]), exp(fearRealNeuMin_BF@bayesFactor[["bf"]][1]),
         exp(fearBothAvNeu_BF@bayesFactor[["bf"]][1]), exp(fearBothAvMin_BF@bayesFactor[["bf"]][1]), exp(fearBothNeuMin_BF@bayesFactor[["bf"]][1]))
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

tableFear <- flextable(tableData[1:9,])
tableFear <- add_header_lines(tableFear, top = TRUE, values = "fear")

save_as_docx(tableFear, path = paste0(pathname, "/tables/tableFear_raw.docx"))



#############################
### Plot for fear ratings ###
#############################
# remove between-subject variance for plotting standard errors based on
# within-subject variance
dataFearWithin <- dataFear[,c("partInd","usGroup","Av_Post","Neu_Post","Min_Post")]
# remove each participant's average from each single value
dataFearWithin[,3:5] <- as.matrix(dataFearWithin[,3:5]) -
  rowMeans(as.matrix(dataFearWithin[,3:5])) 

# prepare data frame for bar plot with means from standard dataset and SE from
# dataset without between-subject variance
plotDataFear <- data.frame(
  usGroup = factor(c(rep("Imagery-Based",3),rep("Classical",3)),
                   levels = c("Imagery-Based","Classical")),
  CS = factor(c("CS+ av","CS+ neu","CS- ","CS+ av","CS+ neu","CS- "),
              levels = c("CS+ av","CS+ neu","CS- ")),
  mean = c(describe(dataFear[dataFear$usGroup == "ima", c(5,8,11)])$mean,
           describe(dataFear[dataFear$usGroup == "real", c(5,8,11)])$mean),
  se = c(describe(dataFearWithin[dataFearWithin$usGroup == "ima", 3:5])$se,
         describe(dataFearWithin[dataFearWithin$usGroup == "real", 3:5])$se)
)

# some general settings
plotFS <- 9
showSig <- TRUE
csLabels = c(expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-",
             expression(paste("CS+"[av])), expression(paste("CS+"[neu])), "CS-")

# bar graphs of group x CS effects on fear ratings
graphFear <- ggplot(data = plotDataFear, aes(x = usGroup, y = mean, fill = CS)) +
  theme_classic() +
  geom_col(aes(fill = CS), position = position_dodge(width = .9)) +
  scale_fill_discrete(type = scico(n = 3, palette = "davos", begin = .1, end = .7)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, width = .1), position = position_dodge(width = .9)) +
  scale_y_continuous(name = "Fear rating (1-5)", limits = c(0.5,5.2), oob = rescale_none, expand = c(0,0)) +
  #labs(title = "Fear") +
  geom_vline(xintercept = 0.41) +
  geom_rect(aes(xmin = 0.4, xmax = 2.6, ymin = 0.45, ymax = 1), fill = "white") +
  geom_hline(yintercept = 1) +
  geom_text(aes(y = 0.8), label = csLabels, position = position_dodge(.9), colour = "black", size = plotFS/.pt) + 
  geom_text(aes(label = usGroup, y = 5), colour = "black", size = plotFS/.pt, fontface = "bold") +
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
    geom_segment(aes(x = 0.7, y = mean+se+.1, xend = 1.0, yend = mean+se+.1), data = plotDataFear[1,]) +
    geom_text(aes(label = "***", x = 0.85, y = mean+se+.15), size = plotFS/2, data = plotDataFear[1,]) +
    geom_segment(aes(x = 0.7, y = mean+se+.4, xend = 1.3, yend = mean+se+.4), data = plotDataFear[1,]) +
    geom_text(aes(label = "**", x = 1.0, y = mean+se+.45), size = plotFS/2, data = plotDataFear[1,]) +
    geom_segment(aes(x = 1.7, y = mean+se+.1, xend = 2.0, yend = mean+se+.1), data = plotDataFear[4,]) +
    geom_text(aes(label = "***", x = 1.85, y = mean+se+.15), size = plotFS/2, data = plotDataFear[4,]) +
    geom_segment(aes(x = 1.7, y = mean+se+.4, xend = 2.3, yend = mean+se+.4), data = plotDataFear[4,]) +
    geom_text(aes(label = "***", x = 2.0, y = mean+se+.45), size = plotFS/2, data = plotDataFear[4,])
}
graphFear



# add margins to subplots
graphFear <- graphFear + theme(plot.title = element_blank(),
                               plot.margin = margin(5,5,5,5))

# saving it
ggsave(filename = paste0(pathname, "/figures/Figure2_barplotFear.eps"),
       plot = graphFear,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)

ggsave(filename = paste0(pathname, "/figures/Figure2_barplotFear.pdf"),
       plot = graphFear,
       width = 100,
       height = 70,
       units = "mm",
       dpi = 300
)
