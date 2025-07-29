# --- author: Christian Panitz
# --- encoding: en_US.UTF-8
# --- R version: 4.3.1 (2023-06-16) -- "Beagle Scouts"
# --- RStudio version: 2023.06.0
# --- script version: Jul 2025
# --- content: demographics (age, gender, handedness, BIS sensitivity)

################
### packages ###
################

# loading required packages
library(psych) # ver. 2.3.9
library(here) # ver. 1.0.1
library(effectsize)
library(BayesFactor)


##################################
### loading and reporting data ###
##################################

# load rating data from text file
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

# Average ZKPQ Neuroticism-Anxiety scores, separate for gender
describeBy(dfRatings$bis, group = dfRatings$gender)

t.test(age ~ group, data = sosciData)
cohens_d(age ~ group, data = sosciData)

t.test(zkpq_nanx ~ group, data = sosciData)
cohens_d(zkpq_nanx ~ group, data = sosciData)