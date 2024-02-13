# readme experimentData

folder contains various txt files with preprocessed single-participant data to be read into R


*******************************************
*** imaCond3_demographicsAndRatings.txt ***
*******************************************

Demographic data (age, gender, handedness) and rating data collected for each participant
First row = variable names
Each following row = 1 participant

Variables:
partCode: alphanumeric participant ID
gender: factor with 2 levels ("female", "male")
age: age in years
handedness: factor with 2 levels ("right", "left")
group: conditioning group, 2 levels ("ima" = imagery-based conditioning group, "real" = classical conditioning group)
cs_permut: values 1-3; permutation set for assignment of face stimulus to CS
cue_permut: values 1-3; permutation set for assignment of geometric cues to imagery (-99: not assigned in class. cond. group)
contAware: contingency awareness, 2 levels ("TRUE", "FALSE")

CS ratings (unpleas_cs* - dis_cs*): values: 1-5
- names: [rating type + CS type + rating index]
- 5 different scales: unpleas(antness), arous(al), anx(iety), ang(er), dis(gust)
- 3 CS: csplus_av(ersive), csplus_neu(tral), csminus
- 4 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition

Unpleasantness ratings for the imagery associated with the different cues; values 0-10 (99 = no imagery associated; only relevant for imagery-based cond. group, therefore NA in classical conditioning group)
- names: unpleas_cue[cue type + rating index]
- 3 cues: cue_av(ersive), cue_neu(tral), cue_minus
- 3 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition

add_learn_trials (only relevant in imagery-based cond. group):
# of additional trials participants needed to correctly report cue-imagery associations
(before the conditioning procedure, every participant had the script read to them and went through two rounds of having presented the contingencies on the screen. Conditioning would only start if they could report the correct contingencies. If they were not able to do so, they had to go additional rounds of instructions, i.e., additional learning trials, until they could report the contingencies correctly)

Unpleasantness ratings for the real US; alues 0-10 (only relevant for class. cond. group ,therefore NA in imagery-based conditioning group)
- names: unpleas_us[cue type + rating index]
- 2 US: us_av(ersive), us_neu(tral)
- 3 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition

Contingency ratings at the end of the experiment in classical conditioning group; values 0-3 ("How often was the [CS] followed by the [US]?" 0 = never, 1 = sometimes, 2 = often, 4 = always)
- names: cont_us[us type + cs type]
- 2 US: us_av(ersive), us_neu(tral)
- 3 CS: csplus_av(ersive), csplus_neu(tral), csminus

Contingency ratings at the end of the experiment in imagery-based conditioning group; values in % ("How often was the [CS] followed by the [Cue]?")
- names: cont_cue[us type + cs type]
- 3 Cues: cue_av(ersive), us_neu(tral), cue_nothing
- 3 CS: csplus_av(ersive), csplus_neu(tral), csminus


*********************************
*** imaCond3_ibimatrix_cs.txt ***
*********************************

Preprocessed IBI responses to CS (-1 to 7 s; in 500 ms bins; baseline-corrected, averaged across trials)
First row = variable names
Each following row = 1 participant

Variables:
partCode: alphanumeric participant ID

IBI values (ibi_*)
names: ibi_[phase + CS + bin]
phase: Hab = habituation; Akq1 = first half of acquisition, Akq2 = second half of acquisition, Akq = all acquisition trials
CS: 51 = CS+aversive, 52 = CS+neutral, 53 = CS-
bin: bl1 (-1000 to -500 ms rel. to CS), bl2 (-500 to 0 ms), 1 through 14 = post-CS bins with size 500 ms



*********************************
*** imaCond3_scrmatrix_cs.txt ***
*********************************

Preprocessed SCR to CS (baseline-corrected peak from 1 to 5 s post-CS; averaged across trials)
First row = variable names
Each following row = 1 participant

Variables:
partCode: alphanumeric participant ID

SCR values (scr_*)
names: ibi_[phase + CS + normalization]
phase: Akqall = all acquisition trials; Akq1 = first half of acquisition, Akq2 = second half of acquisition
CS: 51 = CS+aversive, 52 = CS+neutral, 53 = CS-
normalization: norm = normalized according to Lykken & Venables (1971), raw = not normalized



**************************
*** imaCond3_theta.txt ***
**************************

Preprocessed Theta power to CS (power density from 4 to 8 Hz in time window 0 to 2 s after CS)
First row = variable names
Each following row = 1 participant

Variables:
partCode: alphanumeric participant ID

All other columns
names: [channel + phase + CS]
channel : channel name according to 10-20 system
phase: acqTotal = all acquisition trials; acq1 = first half of acquisition, acq2 = second half of acquisition
CS: csplusav(ersive), csplusneu(tral), csminus



**************************
*** subfolder: erpData ***
**************************

Each file in the subfolder erpData contains the ERP of one participant in one condition
columns = samples (-200 to 1000 ms at 1024 Hz)
rows = channels (channel order as in the channelLocations.txt file)

file names:
[partCode]_IMAKON03_Average_[phase]_EKP_P9P10_[CS].dat
with
partCode: alphanumeric participant ID
phase: Akqu1 = first half of acquisition, Akqu2 = second half of acquisition, Akquges = all acquisition trials
CS: S51 = CS+aversive, S52 = CS+neutral, S53 = CS-
