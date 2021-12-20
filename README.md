# BolkanStoneEtAl
<pre>
 contact: Scott Bolkan sbolkan AT princeton DOT edu 													\
 code to analyze data/generate plots in Bolkan, Stone et al (2022) 											\
 download data from 10.6084/m9.figshare.17299142 or https://figshare.com/s/84695a0cd8cf37a446b9 							\
 see https://github.com/irisstone/glmhmm for general use application of GLM-HMMs to similarly 								\
 structured decision-making behavioral dat, as well as additional publication plots.									\
 
 This repo analyzes data and makes plots for 														\
 (1)  Figure 1      - DMS laser while mice navigate a virtual corridor 											\
 (2)  Figure 2      - cross-task comparison of choice accuracy and motor performance (y-velocity, 							\
 			x-position, view angle, distance) during VR-based									 	\
                      	accumulation of evidence, no distractors, and permanent cues tasks). 								\
 (3)  Figure 3      - Cross-task effects of laser on choice bias in DMS indirect/direct/noOpsin 							\
 			illumination groups, as well as effects of laser on choice bias in 								\
			NAc indirect/direct/noOpsin groups during AoE 											\
 (4)  ExtData Fig1  - inhibition validation (silicon optrode recording in A2a-Cre/D1R-Cre mice 								\
 			expressing NpHR while ambulating on a running wheel)  										\
 (5)  ExtData Fig2  - 					\
 (6)  ExtData Fig3  - 					\
 (7)  ExtData Fig4  - 					\
 (8)  ExtData Fig5  - 					\
 (9)  ExtData Fig6  - 					\
 (10) ExtData Fig7  - 					\
 (11) ExtData Fig8  - 					\
 (12) ExtData Fig9  - 					\
 (13) ExtData Fig12 - 					\
 (14) ExtData Fig14 - 					\
 (15) ExtData Fig15 - 					\
 
 see https://github.com/irisstone/glmhmm for plots of 	\ 
 (1)  Figure 4      -					\
 (2)  Figure 5      -					\
 (3)  Figure 6      -					\
 (4)  Figure 7      -					\
 (5)  ExtData Fig10 -					\
 (6)  ExtData Fig11 - 					\
 (7)  ExtData Fig13 - 					\
 
 
STEP-BY-STEP									 								\
(1) clone/download this repo and add to matlab path 												\
(2) download data on figshare from link at top 													\
	(a) zipped Data folder (12.8GB), after download and unzip (21.9GB),
		take stock of Data folder location 												\
	(c) Contents of Data folder (see BolkanStoneEtAl/Code/utility/logExplanation.m for 							\
		more extensive explanation of behavior logs) : 											\
		(i)   DMS_AoE_GLMHMM_states.mat (DMS D1R and D2R concatenated behavioral logs 							\
			with associated GLM-HMM stateIdx and stateProb)										\
		(ii)  OffOn_TasksOrGroup_all.mat (DMS D1R,D2R,NoOpsin and NAc D1R,D2R,NoOpsin 							\
			concatentated behavioral logs in mice performing the evidence 								\
			accumulation (AoE), no distrators (nd) or permanent cues (pc) tasks) 							\
		(iii) rtPPdata.mat (real-time conditioned place preference, Ethovision output 							\
			measures, DMS D1R,D2R,NoOpsin laser) 											\
		(iv)  runningWheelSpikes.mat (summary of 110 single-units in DMS from D1R-Cre/A2a-Cre 						\
			mice expressing DIO-NpHR, laser sweeps while freely moving on running wheel) 						\
		(v)   taskShapingData.mat (summary of AoE, N.D., and P.C. task shaping performance)						\
		(vi)  virtualCorridor.mat (DMS D1R,D2R,NoOpsin concatenated behavioral logs in mice 						\
			navigating a virtual corridor)      									         	\
		(vii) 4 .xls docs summarizing stereological quantification of D1R/D2R in situ 							\
			hybridization expression in A2a/D2R/D1R-Cre mouse lines							          	\	
(3) in MATLAB open .m file .../BolkanStoneEtAl/Code/utility/globalParams.m                                                            		\
	(a) change globalParams.dataPath to match location of downloaded data                                                         		\
(4) run e.g. figure1_script in matlab command line                                                                                    		\
	(a) will generate all figure1 data plots and save in your matlab workspace the 								\
		analyzed variables used to generate each plot  											\
(5) otherwise, edit figure1_script from matlab command line                                                                           		\
	(a) each script is bracketed in sections that can be run in steps to                                                          		\
		(i)   load pre-processed source data                                                                                  		\
		(ii)  analyze data in commented sections                                                                              		\
		(iii) plot data of given analysis                                                                                     		\			
