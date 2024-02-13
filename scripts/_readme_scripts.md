# readme scripts

R scripts for statistical analyses and plotting

00_demographics.R
- compute mean age and distribution of gender and handedness, across and within conditioning groups

01_fearRatings.R
- compares mean fear ratings between CS types at the end of acquisition
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots bar plot for comparing fear ratings between all Conditioning Group x CS Type cells; saves it in the "figures" folder
- saves results from t-tests in "tables" folder

01supp_fearRatings_aware.R
- supplementary analyses: repeats main analyses on fear ratings (01_fearRatings.R) only in contingency-aware participants
- saves line plots and ANOVA/t-test output in "supplement" folder

01supp_fearRatings_timeFactor.R
- supplementary analyses: repeats main analyses on fear ratings (01_fearRatings.R) adding the factor Time (pre-acquisition, mid-acquisition, post-acquisition)
- CS Type x Time ANOVAs within conditioning groups + Conditioning Group x CS Type x Time ANOVA across all participants
- follow-up t-tests within groups compare CS responses separately for Time levels
- follow-up t-tests across groups compare differential CS responses separately for Time levels
- saves line plots and ANOVA/t-test output in "supplement" folder

01supp_otherRatings.R
compares mean unpleasantness, arousal, anger, and disgust ratings between CS types at the end of acquisition; reported in the supplement
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots bar plots for comparing all ratings between all Conditioning Group x CS Type cells; saves them in the "supplement folder"

02_scr.R
- compares mean normalized (Lykken & Vanables, 1971) SCR between CS types averaged across all acquisition trials
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots bar plot for comparing SCRs between all Conditioning Group x CS Type cells; saves it in the "figures" folder
- saves results from t-tests in "tables" folder

02supp_scr_aware.R
- supplementary analyses: repeats main analyses on SCRs (02_scr.R) only in contingency-aware participants
- saves line plots and ANOVA/t-test output in "supplement" folder

02supp_scr_timeFactor.R
- supplementary analyses: repeats main analyses on SCRs (02_scr.R) adding the factor Time (1st half of trials, 2nd half of trials)
- CS Type x Time ANOVAs within conditioning groups + Conditioning Group x CS Type x Time ANOVA across all participants
- follow-up t-tests within groups compare CS responses separately for Time levels
- follow-up t-tests across groups compare differential CS responses separately for Time levels
- saves line plots and ANOVA/t-test output in "supplement" folder

03_ibi.R
- compares mean IBI change (rel. to baseline; mean IBI 2 to 5 s post-CS) between CS types averaged across all acquisition trials
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots IBI time courses and bar plots for means, comparing IBI between all Conditioning Group x CS Type cells; saves them in the "figures" folder
- saves results from t-tests in "tables" folder

03b_ibi_exploratoryTimeWin
- repeats main analyses on IBI (03_ibi.R) using exploratory time window for averaging (4 to 7 s)

03supp_ibi_aware.R
- supplementary analyses: repeats main analyses on IBI (03_ibi.R) only in contingency-aware participants
- saves line plots and ANOVA/t-test output in "supplement" folder

03supp_ibi_timeFactor.R
- supplementary analyses: repeats main analyses on IBI (03_ibi.R) adding the factor Time (1st half of trials, 2nd half of trials)
- CS Type x Time ANOVAs within conditioning groups + Conditioning Group x CS Type x Time ANOVA across all participants
- follow-up t-tests within groups compare CS responses separately for Time levels
- follow-up t-tests across groups compare differential CS responses separately for Time levels
- saves line plots and ANOVA/t-test output in "supplement" folder

04_LPP.R
- compares mean LPP amplitude at Pz (300 to 700 ms post-CS) between CS types averaged across all acquisition trials
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots ERPs at Pz and bar plots for mean amplitudes, comparing LPP between all Conditioning Group x CS Type cells; saves them in "figures" folder
- saves results from t-tests in "tables" folder

04b_LPP_exploratoryTimeWin
- repeats main analyses on LPP (04_LPP.R) using exploratory time window (300 to 1000 ms) and channel selection (Pz, POz, P1, P2, P3, P4, PO3, PO4) for averaging
- saves plots and ANOVA/t-test output in "supplement" folder

04supp_LPP_aware.R
- supplementary analyses: repeats main analyses on LPP (04_LPP.R) only in contingency-aware participants
- saves line plots and ANOVA/t-test output in "supplement" folder

04supp_LPP_timeFactor.R
- supplementary analyses: repeats main analyses on LPP (04_ibi.R) adding the factor Time (1st half of trials, 2nd half of trials)
- CS Type x Time ANOVAs within conditioning groups + Conditioning Group x CS Type x Time ANOVA across all participants
- follow-up t-tests within groups compare CS responses separately for Time levels
- follow-up t-tests across groups compare differential CS responses separately for Time levels
- saves line plots and ANOVA/t-test output in "supplement" folder

05supp_theta.R
- supplementary analyses: compares mean theta power density (4 to 8 Hz) at Fz (0 to 2 s post-CS) between CS types averaged across all acquisition trials
- computes one-factorial ANOVAs with factor CS Type (CS+av, CS+neu, CS-) within conditioning groups + follow-up t-tests for pairwise CS comparisons
- computes two-factorial ANOVA with factors Conditioning Group (imagery-based, classical) and CS Type (CS+av, CS+neu, CS-) + follow-up t-tests for comparisons of differential CS responses between groups
- plots line plots for mean amplitudes, comparing LPP between all Conditioning Group x CS Type cells; saves them in "supplement" folder

05supp_theta_timeFactor.R
- supplementary analyses: repeats analyses on theta (04supp_theta.R) adding the factor Time (1st half of trials, 2nd half of trials)
- CS Type x Time ANOVAs within conditioning groups + Conditioning Group x CS Type x Time ANOVA across all participants
- follow-up t-tests within groups compare CS responses separately for Time levels
- follow-up t-tests across groups compare differential CS responses separately for Time levels
- saves line plots and ANOVA/t-test output in "supplement" folder

05supp_theta_topographies.R
- plots topographies for power density in theta band (4 to 8 Hz) separately for the two conditioning groups: average across all CS, [difference CS+av - CS+neu], and difference [CS+av - CS-]
- saves them in the "supplement folder"