# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# --- RStudio version: 1.3.1093
# --- script version: Jan 2023
# --- content: demographics (age, gender, & handedness)

################
### packages ###
################

# (installing and) loading required packages
# will install newest package version, not necessarily the version originally used!
if(!is.element("psych",installed.packages()[,1])) {install.packages("psych")}
  library(psych)
if(!is.element("here",installed.packages()[,1])) {install.packages("here")}
library(here)



##################################
### loading and reporting data ###
##################################

# load rating data from text file
# (see imaCond3_demogrphicsAndRatings_readme.txt for more details)
pathname <- here()
dfRatings <- read.csv(paste0(pathname, "/experimentData/imaCond3_demographicsAndRatings.txt"), sep=",")

# descriptives for age across and separate for conditioning groups
describe(dfRatings$age)
describeBy(dfRatings$age, group = dfRatings$group)

# frequencies for self-reported gender across and separate for conditioning groups
table(dfRatings$gender)
table(dfRatings$gender, dfRatings$group)

# frequencies for self-reported handedness across and separate for conditioning groups
table(dfRatings$handedness)
table(dfRatings$handedness, dfRatings$group)