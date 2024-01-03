Data collected by presentation for each participant
First row = variable names
Each following row = 1 participant

Variables:
vpcode: alphanumeric participant ID
group: values 0-1 (0 = imagery-based conditioning group, 1 = classical conditioning group)
cs_permut: values 1-3; permutation set for assignment of face stimulus to CS
cue_permut: values 1-3; permutation set for assignment of geometric cues to imagery (-99: not assigned in class. cond. group)

CS ratings (unpleas_cs* - dis_cs*): values: 1-5
- names: [rating type + CS type + rating index]
- 5 different ratings (unpleasantness, arousal, anxiety, anger, disgust)
- 3 CS: csplus_av(ersive), csplus_neu(tral), csminus
- 4 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition

Unpleasantness ratings for the imagery associated with the different cues (only relevant for imagery-based cond. group)
- names: unpleas_cue[cue type + rating index]
- 3 cues: cue_av(ersive), cue_neu(tral), cue_minus
- 3 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition

add_learn_trials (only relevant in imagery-based cond. group):
# of additional trials participants needed to correctly report cue-imagery associations
(before the conditioning procedure, every participants had the script read to them and went through two rounds of having presented the contingencies on the screen. Conditioning would only start if they could report the correct contingencies. If they were not able to do so, they had to go additional rounds of instructions, i.e., additional learning trials, until they could report the contingencies correctly)

Unpleasantness ratings for the real US (only relevant for class. cond. group)
- names: unpleas_us[cue type + rating index]
- 2 US: cue_av(ersive), cue_neu(tral)
- 3 rating indices: 1 = before habituation; 2 = between habituation and acquisition; 3 = between 1st & 2nd half of acquisition; 4 = after acquisition


Ratings:
vpcode group cs_permut cue_permut unpleas_csplus_av_1 unpleas_csplus_av_2 unpleas_csplus_av_3 unpleas_csplus_av_4 unpleas_csplus_neu_1 unpleas_csplus_neu_2 unpleas_csplus_neu_3 unpleas_csplus_neu_4 unpleas_csminus_1 unpleas_csminus_2 unpleas_csminus_3 unpleas_csminus_4 arous_csplus_av_1 arous_csplus_av_2 arous_csplus_av_3 arous_csplus_av_4 arous_csplus_neu_1 arous_csplus_neu_2 arous_csplus_neu_3 arous_csplus_neu_4 arous_csminus_1 arous_csminus_2 arous_csminus_3 arous_csminus_4 anx_csplus_av_1 anx_csplus_av_2 anx_csplus_av_3 anx_csplus_av_4 anx_csplus_neu_1 anx_csplus_neu_2 anx_csplus_neu_3 anx_csplus_neu_4 anx_csminus_1 anx_csminus_2 anx_csminus_3 anx_csminus_4 ang_csplus_av_1 ang_csplus_av_2 ang_csplus_av_3 ang_csplus_av_4 ang_csplus_neu_1 ang_csplus_neu_2 ang_csplus_neu_3 ang_csplus_neu_4 ang_csminus_1 ang_csminus_2 ang_csminus_3 ang_csminus_4 dis_csplus_av_1 dis_csplus_av_2 dis_csplus_av_3 dis_csplus_av_4 dis_csplus_neu_1 dis_csplus_neu_2 dis_csplus_neu_3 dis_csplus_neu_4 dis_csminus_1 dis_csminus_2 dis_csminus_3 dis_csminus_4 unpleas_cue_av_1 unpleas_cue_av_2 unpleas_cue_av_3 unpleas_cue_neu_1 unpleas_cue_neu_2 unpleas_cue_neu_3 unpleas_cue_minus_1 unpleas_cue_minus_2 unpleas_cue_minus_3 add_learn_trials unpleas_us_av_1 unpleas_us_av_2 unpleas_us_av_3 unpleas_us_neu_1 unpleas_us_neu_2 unpleas_us_neu_3 cont_usav_csplusav cont_usav_csplusneu cont_usav_csminus cont_usneu_csplusav cont_usneu_csplusneu cont_usneu_csminus
EK08ER20 0 2 2 3 1 1 2 3 2 2 3 4 3 2 3 1 2 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 4 2 2 3 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 99 3 4 1 2 1 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
EL07CH15 1 2 -99 3 1 4 5 2 1 1 1 4 4 2 1 2 2 5 5 3 1 2 1 4 4 2 2 1 1 5 4 1 2 2 2 3 2 1 2 1 1 5 4 1 1 2 1 1 2 1 1 1 1 3 3 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 7 9 10 0 0 0 2 1 0 0 3 0 
ER07AS29 0 3 3 3 3 4 4 4 3 2 2 3 3 2 2 2 2 3 3 3 3 4 4 1 1 1 1 2 2 3 3 3 2 2 1 1 1 1 1 1 1 3 3 2 2 1 1 1 1 1 1 2 1 2 2 3 3 1 1 3 2 2 2 5 5 6 0 2 2 0 0 0 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
LT07EL10 0 1 3 3 2 3 2 3 3 1 2 2 2 2 2 3 2 4 1 3 2 3 1 3 2 1 1 2 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 3 1 1 1 1 1 1 1 1 1 3 3 1 0 0 0 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
AY05LI06 0 2 3 3 3 4 4 2 2 2 3 3 4 3 2 3 3 4 4 3 3 3 3 2 4 4 3 3 1 3 2 1 1 1 1 2 1 2 1 1 1 5 4 1 1 1 2 1 2 3 1 2 1 4 3 1 1 2 4 1 1 1 2 9 10 9 99 99 7 0 1 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
CH08AN26 0 3 3 3 3 4 4 4 2 3 2 3 3 2 1 1 3 4 5 3 1 2 3 2 2 3 2 1 2 3 4 2 2 1 2 1 1 1 1 1 1 2 4 1 1 1 2 1 1 1 1 1 1 2 3 1 1 1 1 1 1 1 1 3 6 8 1 4 4 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
ST06RT16 1 1 -99 4 5 5 5 3 4 2 2 2 2 3 3 4 4 4 4 1 3 1 3 3 3 1 2 2 2 3 3 2 1 1 1 3 2 1 2 2 3 5 5 2 1 1 2 1 1 2 1 4 4 4 4 3 2 2 3 2 1 2 2 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 10 9 1 3 5 2 0 0 0 2 0 
OP05IG02 0 2 3 3 3 3 3 2 2 2 2 4 3 3 3 1 1 1 1 2 1 2 1 1 1 1 1 1 1 2 2 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 3 2 1 1 6 6 5 0 0 0 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
LD08LF01 1 3 -99 3 2 3 4 4 5 3 2 3 4 2 2 1 3 5 3 3 4 3 2 2 3 2 2 1 1 2 2 3 3 1 1 2 2 1 1 1 1 3 2 2 2 1 1 1 1 1 1 2 1 1 2 4 4 3 2 1 3 2 2 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 9 9 0 0 0 2 0 0 0 3 0 
IN09LF17 1 2 -99 3 3 3 4 3 2 2 2 4 4 4 3 2 1 4 4 1 1 2 1 3 2 1 2 1 1 2 2 1 1 1 1 2 2 1 1 1 1 1 4 1 1 1 1 2 2 1 1 2 2 2 1 2 2 2 1 3 3 2 2 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 9 10 1 0 0 2 0 0 0 2 0 
KE05RD15 1 2 -99 3 3 4 4 4 1 3 3 4 4 4 4 2 2 3 3 3 4 3 3 4 4 3 4 1 2 2 2 1 1 2 1 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 10 10 10 1 2 1 1 1 3 1 2 1 
AL05IM24 0 1 2 3 4 3 3 4 4 4 4 3 2 2 2 1 3 3 3 1 3 3 3 2 2 2 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 1 3 3 3 3 3 1 1 1 1 2 2 2 2 3 4 4 3 1 1 1 1 8 7 8 4 2 3 0 0 99 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
KE06ER07 1 3 -99 2 2 4 4 3 4 3 2 4 3 3 2 2 3 3 4 3 3 2 3 3 2 1 1 1 1 2 3 2 2 2 1 1 1 1 1 1 1 4 4 1 1 1 1 1 1 1 1 1 1 3 4 1 2 1 1 1 2 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 8 8 9 0 0 0 2 0 0 0 2 0 
LZ05AN21 1 2 -99 3 4 4 5 4 3 2 2 4 4 3 2 2 3 5 5 2 3 3 2 4 4 2 1 2 2 4 4 2 2 1 1 1 2 2 1 1 3 3 4 4 2 2 1 3 3 3 3 3 4 3 4 2 3 2 2 5 5 3 3 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 7 8 7 0 0 0 2 0 0 0 2 0 
EN05AN11 1 3 -99 2 2 5 5 3 4 2 2 4 4 2 1 1 1 5 5 2 3 3 2 3 3 1 1 1 1 3 4 1 1 1 1 1 1 1 1 1 1 3 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 9 10 0 2 2 2 0 0 0 2 0 
EN07LF06 0 3 3 5 3 2 3 3 4 3 3 3 3 2 3 2 1 1 4 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 8 4 3 1 1 1 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
LE06KO07 0 1 3 4 5 1 5 2 3 3 3 3 4 4 3 1 1 2 1 2 2 2 2 1 2 2 2 2 3 4 4 1 2 2 3 4 3 3 2 2 3 2 3 1 1 1 1 2 2 2 2 3 5 4 5 2 2 3 3 4 4 3 3 3 5 5 1 2 4 0 2 2 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
ER07LD07 1 3 -99 4 3 5 4 4 4 3 3 3 4 3 2 1 2 5 5 1 2 4 3 1 2 2 2 3 1 3 3 3 3 2 2 1 2 2 1 2 1 2 2 2 3 2 1 1 1 1 1 2 1 1 1 2 2 1 1 2 2 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 8 8 9 1 2 0 2 0 1 0 2 0 
TZ03UL29 1 1 -99 3 4 4 4 4 4 4 3 2 2 1 1 4 3 4 3 2 3 3 1 2 1 1 1 2 3 3 3 2 2 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 4 2 3 4 4 2 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 8 7 7 3 0 0 2 0 0 0 3 0 
LL06ER17 0 3 2 1 1 4 3 3 3 4 4 3 3 3 3 3 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 3 2 1 1 1 1 1 1 2 2 2 3 4 3 1 2 1 1 1 1 3 2 1 1 2 1 1 1 1 1 99 2 1 0 1 99 0 1 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
HA04IH18 1 3 -99 3 3 5 5 4 5 2 1 3 4 2 1 1 2 5 5 3 3 3 1 2 2 1 1 1 1 4 5 3 4 1 1 1 1 1 1 1 1 5 5 3 3 1 1 1 1 1 1 1 1 5 5 2 2 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 1 10 10 4 0 0 2 0 0 0 3 0 
NE08IM07 0 3 2 3 3 4 4 3 3 3 3 4 5 2 1 2 1 2 2 1 3 2 4 3 1 1 2 1 1 4 1 2 4 3 4 5 5 2 1 1 1 2 5 1 1 1 4 1 1 1 1 1 3 3 5 1 4 2 4 5 5 1 1 3 3 6 0 7 1 99 1 8 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
CH05ER11 0 1 1 4 4 5 5 2 2 3 3 1 1 1 1 3 4 4 4 4 3 3 2 4 4 2 2 4 3 4 5 1 1 2 2 1 1 1 1 3 4 5 5 1 1 3 3 1 1 1 1 4 4 5 5 2 2 3 2 1 1 1 1 9 9 10 2 3 3 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
ER06RD12 1 1 -99 4 4 4 5 4 4 4 3 2 2 1 1 1 3 5 5 2 4 3 2 1 1 1 2 3 4 4 4 3 3 1 1 1 1 1 1 2 3 3 3 4 4 2 3 1 1 1 1 1 1 3 3 1 3 3 2 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 8 9 9 0 0 0 2 0 0 0 2 0 
SS08NN01 0 2 1 3 3 3 3 4 2 1 2 4 4 4 4 1 1 1 1 3 2 1 1 1 1 1 1 1 1 2 2 2 1 1 1 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 3 2 2 2 9 9 5 0 0 2 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
UN06AS18 0 1 2 4 4 4 4 4 4 3 3 3 3 3 2 1 1 3 3 1 1 1 1 1 1 1 1 2 1 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 6 6 3 3 2 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
AN06US04 1 1 -99 5 4 4 4 5 4 3 3 3 3 3 3 4 4 4 5 3 2 3 5 3 2 1 1 4 2 2 2 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 10 10 8 3 0 0 2 0 0 0 2 0 
IL06IN05 1 1 -99 4 4 3 3 3 3 3 3 2 2 2 3 4 4 4 5 1 1 2 4 3 2 1 2 4 3 3 4 1 1 1 3 1 1 1 1 3 1 2 4 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 10 10 0 1 2 2 0 0 1 2 0 
EZ06EN31 0 2 1 3 3 4 3 3 2 3 3 4 4 2 2 2 2 3 2 2 3 3 4 3 3 2 2 2 2 2 2 3 1 2 1 3 2 2 1 2 2 3 3 2 1 2 1 3 2 1 1 2 1 2 2 1 1 1 1 2 2 1 1 9 5 4 2 1 1 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
AU07NS21 1 2 -99 3 3 3 3 2 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 2 1 1 1 1 2 3 1 1 2 2 2 1 2 1 1 1 1 1 1 1 1 1 3 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 7 10 4 5 4 2 0 1 1 2 0 
EN09RD15 0 1 2 2 3 3 3 1 3 4 4 4 5 5 4 1 1 2 3 2 1 3 4 3 3 4 3 1 1 1 2 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 5 4 4 3 1 1 1 2 1 1 1 1 4 4 4 1 5 10 6 0 5 2 0 5 2 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
IL09ER08 1 2 -99 3 3 4 5 2 2 2 1 3 2 2 1 1 1 5 4 1 1 3 3 1 1 1 1 1 1 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 10 10 0 0 0 2 0 0 0 2 0 
TA08RT19 1 3 -99 3 3 4 5 3 3 3 3 3 3 3 3 3 3 2 2 3 1 2 1 3 1 1 1 3 3 1 1 3 1 1 1 3 1 1 1 3 3 4 4 3 1 1 1 3 1 1 1 3 3 1 1 3 1 1 1 3 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 5 6 8 3 4 5 1 0 0 0 1 0 
AB09ER21 0 3 2 4 5 5 5 3 2 3 3 2 2 2 2 3 4 5 4 2 1 1 1 1 1 1 1 2 3 4 2 1 1 1 1 1 1 1 1 2 3 4 3 1 1 1 1 1 1 1 1 3 4 5 4 2 1 1 1 1 1 1 1 8 9 8 4 1 0 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
IN05IM17 0 3 1 3 3 5 5 4 4 3 3 3 2 3 3 2 3 4 4 3 3 4 3 1 2 2 2 2 3 4 5 3 3 3 2 2 2 1 3 1 2 3 4 1 2 2 3 1 2 2 2 2 3 3 4 2 2 2 3 2 2 2 1 8 8 9 1 0 3 0 0 1 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
BE06ND04 0 1 1 3 4 4 4 3 3 2 2 4 4 3 3 1 1 1 1 1 1 2 3 1 1 2 1 3 3 3 3 1 1 1 1 2 2 1 1 1 2 3 3 1 1 1 1 1 1 1 1 1 2 3 2 1 1 1 1 3 3 2 2 5 5 5 0 1 0 0 1 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
IX07AS13 1 2 -99 4 5 5 5 2 2 2 2 3 4 4 4 2 2 3 3 2 1 1 1 2 1 1 1 3 4 5 5 1 1 1 1 4 4 2 2 1 3 5 5 1 1 1 1 3 5 2 2 2 1 5 5 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 10 10 1 0 0 2 0 0 0 2 0 
UN05WE18 0 2 2 4 3 4 4 3 2 2 2 4 4 3 4 2 3 3 3 1 2 2 3 3 2 3 4 3 3 3 4 2 1 1 1 3 3 3 4 1 1 1 1 1 1 1 1 1 1 1 1 2 3 2 3 1 2 2 2 3 3 2 3 9 9 8 3 1 3 99 0 2 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
NN05CH10 1 3 -99 2 2 1 3 3 4 3 2 2 2 1 1 2 3 4 4 1 3 2 2 1 1 1 1 1 1 3 3 2 2 2 1 1 1 1 1 1 1 4 4 1 2 1 1 1 1 1 1 1 1 2 2 4 3 2 2 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 6 7 7 0 0 0 2 0 0 0 2 0 
ER05EL02 0 1 1 4 4 4 3 3 3 2 2 2 1 1 2 4 4 3 3 2 2 2 2 2 3 1 1 2 3 3 2 3 1 1 1 2 1 1 1 1 2 2 4 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 3 4 6 1 2 0 0 0 1 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
VA06ID03 1 1 -99 2 3 4 3 4 4 2 4 3 4 2 3 1 1 4 4 3 2 2 4 2 3 1 1 1 1 3 2 2 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 8 7 7 0 0 0 2 3 0 0 3 0 
NN09AS30 0 2 1 3 3 4 3 3 3 2 3 4 4 4 3 1 1 1 1 3 2 1 2 2 2 2 1 1 2 3 2 1 1 1 1 3 3 3 4 1 2 1 3 1 1 1 1 3 3 3 2 3 2 3 3 1 2 1 1 3 3 3 3 6 3 2 0 1 0 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
MY05EN23 1 3 -99 4 3 4 4 3 4 2 2 3 3 2 2 2 2 3 4 2 2 2 3 1 2 1 1 2 1 2 2 1 2 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 8 9 0 0 0 2 0 0 0 2 0 
NN09LD03 1 2 -99 4 4 4 3 2 3 3 2 3 2 3 3 3 3 3 1 2 3 3 1 2 2 3 1 4 3 2 2 2 2 2 1 2 3 2 2 4 4 4 4 1 3 4 3 3 2 4 4 3 3 3 2 1 2 2 1 1 1 2 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 6 3 3 1 0 0 2 0 0 0 2 0 
CK06RT15 0 3 1 3 4 4 5 4 5 3 4 3 2 2 2 3 3 4 3 3 3 3 3 3 3 3 1 4 3 3 4 3 3 3 4 3 2 1 1 3 4 5 5 3 4 4 3 1 3 1 1 4 4 4 4 4 4 3 4 2 2 1 1 10 9 10 2 4 2 0 0 0 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
HM08LE07 1 1 -99 4 3 3 3 3 3 2 2 3 3 3 3 1 1 3 3 1 2 3 3 2 3 3 3 2 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 3 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 9 8 7 0 0 0 2 0 0 0 2 0 
IG08AR03 1 1 -99 5 5 2 4 4 3 3 2 4 2 1 1 2 2 3 4 1 1 2 3 1 3 1 1 3 2 2 3 1 1 1 1 1 1 1 1 4 4 3 4 3 3 1 2 2 1 1 1 5 5 4 4 4 4 3 2 1 2 1 1 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 6 7 6 0 0 0 2 0 0 0 2 0 
JA05RL22 0 2 3 3 3 4 3 4 4 3 3 4 4 3 3 2 1 1 1 2 1 1 1 1 1 1 1 3 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 2 2 2 0 99 99 0 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 -99 
