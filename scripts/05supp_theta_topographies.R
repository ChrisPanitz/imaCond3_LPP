# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Mar 2024
# --- content: Plotting topographies for frontomedial Theta (across CS and pairwise comparisons within experimental groups)

###################
### preparing R ###
###################

# loading required packages
library(tidyr) # ver. 1.3.0
library(ggplot2) # ver. 3.4.2
library(ggpubr) # ver. 0.6.0
library(here) # ver. 1.0.1
library(eegUtils) # ver. 0.7.0


########################
### data preparation ###
########################

# load rating data from text file
pathname <- here()
importRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# load Theta data from text file
# (see imaCond3_theta_readme.txt for more details)
importTheta <- read.csv(paste0(pathname, "/experimentData/imaCond3_theta.txt"), sep = ",")
importTheta <- merge(importRatings[,c("partCode","group")], importTheta, by = "partCode")
importTheta <- importTheta[order(importTheta$partCode),]
importThetaLong <- pivot_longer(importTheta, cols = Fp1_acqTotal_csplusav:O2_acq2_csminus, 
                                names_to = c("channel", "time", "CS"), names_sep = "_", values_to = "theta")
# select all mean theta across all trials
importThetaLong <- importThetaLong[importThetaLong$time == "acqTotal", c("partCode","group","channel","CS","theta")]

topoAvgAVima <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "ima" & importThetaLong$CS == "csplusav", ], FUN = "mean")
topoAvgNEUima <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "ima" & importThetaLong$CS == "csplusneu", ], FUN = "mean")
topoAvgMINima <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "ima" & importThetaLong$CS == "csminus", ], FUN = "mean")

topoAvgAVreal <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "real" & importThetaLong$CS == "csplusav", ], FUN = "mean")
topoAvgNEUreal <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "real" & importThetaLong$CS == "csplusneu", ], FUN = "mean")
topoAvgMINreal <- aggregate(theta ~ channel, data = importThetaLong[importThetaLong$group == "real" & importThetaLong$CS == "csminus", ], FUN = "mean")

# Load channel locations and transform from 3D theta + radius into 2D x & y
# for plotting purposes
loadname <- paste0(pathname,"/channelLocations/chanLocs_biosemi64.txt")
chanLocs <- read.csv(loadname, sep = ";")
chanLocs$thetaRadian <- pi/180*chanLocs$theta
chanLocs$x <- chanLocs$radius*sin(chanLocs$thetaRadian)*200
chanLocs$y <- chanLocs$radius*cos(chanLocs$thetaRadian)*200
names(chanLocs) <- gsub("name","electrode",names(chanLocs))


# number of electrodes
nrChans = length(topoAvgAVima)

# Create data frame with factors lab, driving frequency & modulation function, 
# electrode name, x & y coordinates for plot & theta power for different conditions

dfTopos <- data.frame(
  electrode = topoAvgAVima$channel,
  ima_allCS = rowMeans(cbind(topoAvgAVima$theta, topoAvgNEUima$theta, topoAvgMINima$theta)),
  ima_DiffAvNeu = topoAvgAVima$theta - topoAvgNEUima$theta,
  ima_DiffAvMin = topoAvgAVima$theta - topoAvgMINima$theta,
  real_allCS = rowMeans(cbind(topoAvgAVreal$theta, topoAvgNEUreal$theta, topoAvgMINreal$theta)),
  real_DiffAvNeu = topoAvgAVreal$theta - topoAvgNEUreal$theta,
  real_DiffAvMin = topoAvgAVreal$theta - topoAvgMINreal$theta
)

dfTopos <- merge(chanLocs, dfTopos, by = "electrode")



# settings for topography plots
topoRes <- 200
chanCol <- "black"
nrColors <- 8

minLimAll <- min(c(dfTopos$ima_allCS, dfTopos$real_allCS))
maxLimAll <- max(c(dfTopos$ima_allCS, dfTopos$real_allCS))

minLimDiff <- min(c(dfTopos$ima_DiffAvNeu, dfTopos$ima_DiffAvMin, dfTopos$real_DiffAvNeu, dfTopos$real_DiffAvMin))
maxLimDiff <- max(c(dfTopos$ima_DiffAvNeu, dfTopos$ima_DiffAvMin, dfTopos$real_DiffAvNeu, dfTopos$real_DiffAvMin))
absLimDiff <- max(abs(c(minLimDiff, maxLimDiff)))

### create topography plots

topoAllIma <- topoplot(data = dfTopos,
                    contour = FALSE, interp_limit = "head", highlights = "Fz",
                    grid_res = topoRes, quantity = "ima_allCS", 
                    limits = c(-maxLimAll, maxLimAll), method = "gam")
topoAllIma$guides$fill$title <- "Power Density (µV^2/Hz)"
topoAllIma$guides$fill$barheight <- unit(9, "lines")


topoAllReal <- topoplot(data = dfTopos,
                       contour = FALSE, interp_limit = "head", highlights = "Fz",
                       grid_res = topoRes, quantity = "real_allCS", 
                       limits = c(-maxLimAll, maxLimAll), method = "gam")
topoAllReal$guides$fill$title <- "Power Density (µV^2/Hz)"
topoAllReal$guides$fill$barheight <- unit(9, "lines")

topoAllImaAvNeu <- topoplot(data = dfTopos,
                       contour = FALSE, interp_limit = "head", highlights = "Fz",
                       grid_res = topoRes, quantity = "ima_DiffAvNeu", 
                       limits = c(-absLimDiff, absLimDiff), method = "gam")
topoAllImaAvNeu$guides$fill$title <- "Diff. in Power Density (µV^2/Hz)"
topoAllImaAvNeu$guides$fill$barheight <- unit(12, "lines")

topoAllImaAvMin <- topoplot(data = dfTopos,
                            contour = FALSE, interp_limit = "head", highlights = "Fz",
                            grid_res = topoRes, quantity = "ima_DiffAvMin", 
                            limits = c(-absLimDiff, absLimDiff), method = "gam")
topoAllImaAvMin$guides$fill$title <- "Diff. in Power Density (µV^2/Hz)"
topoAllImaAvMin$guides$fill$barheight <- unit(12, "lines")

topoAllRealAvNeu <- topoplot(data = dfTopos,
                            contour = FALSE, interp_limit = "head", highlights = "Fz",
                            grid_res = topoRes, quantity = "real_DiffAvNeu", 
                            limits = c(-absLimDiff, absLimDiff), method = "gam")
topoAllRealAvNeu$guides$fill$title <- "Diff. in Power Density (µV^2/Hz)"
topoAllRealAvNeu$guides$fill$barheight <- unit(12, "lines")

topoAllRealAvMin <- topoplot(data = dfTopos,
                            contour = FALSE, interp_limit = "head", highlights = "Fz",
                            grid_res = topoRes, quantity = "real_DiffAvMin", 
                            limits = c(-absLimDiff, absLimDiff), method = "gam")
topoAllRealAvMin$guides$fill$title <- "Diff. in Power Density (µV^2/Hz)"
topoAllRealAvMin$guides$fill$barheight <- unit(12, "lines")

# add margins
topoAllIma <- topoAllIma + theme(plot.margin = unit(c(5,5,5,5), "mm"))
topoAllReal <- topoAllReal + theme(plot.margin = unit(c(5,5,5,5), "mm"))
topoAllImaAvNeu <- topoAllImaAvNeu + theme(plot.margin = unit(c(5,5,5,5), "mm"))
topoAllImaAvMin <- topoAllImaAvMin + theme(plot.margin = unit(c(5,5,5,5), "mm"))
topoAllRealAvNeu <- topoAllRealAvNeu + theme(plot.margin = unit(c(5,5,5,5), "mm"))
topoAllRealAvMin <- topoAllRealAvMin + theme(plot.margin = unit(c(5,5,5,5), "mm"))


allTopos <- ggarrange(topoAllIma, topoAllImaAvNeu, topoAllImaAvMin,
                      topoAllReal, topoAllRealAvNeu, topoAllRealAvMin,
                      labels = c("Imagery-based: Average across all CS",
                                 "   Imagery-based: CS+av vs CS+neu   ",
                                 "    Imagery-based: CS+av vs CS-     ",
                                 "  Classical: Average across all CS  ",
                                 "     Classical: CS+av vs CS+neu     ",
                                 "      Classical: CS+av vs CS-       "),
                      hjust = -0.2, nrow = 2, ncol = 3)
allTopos

# save it
ggsave(filename = paste0(pathname, "/supplement/05s_theta_topographies.pdf"),
       plot = allTopos,
       width = 400,
       height = 200,
       units = "mm"
)
