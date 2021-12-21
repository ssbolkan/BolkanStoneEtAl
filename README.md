# BolkanStoneEtAl
<pre>
 ***contact: Scott Bolkan sbolkan AT princeton DOT edu*** 											\
 --code to analyze data/generate plots in Bolkan, Stone et al (2022) 										\
 --download data from 10.6084/m9.figshare.17299142 or https://figshare.com/s/84695a0cd8cf37a446b9 						\
 --see https://github.com/irisstone/glmhmm for general use application of GLM-HMMs to similarly 						\
 structured decision-making behavioral data, as well as additional publication plots.								\
 
 This repo analyzes data and makes plots for 													\
 (1)  Figure 1      - motor effects of DMS laser while mice navigate a virtual corridor 							\
 (2)  Figure 2      - cross-task comparison of choice accuracy and motor performance (y-velocity, 						\
 			x-position, view angle, distance) during VR-based								 	\
                      	accumulation of evidence, no distractors, and permanent cues tasks). 							\
 (3)  Figure 3      - Cross-task effects of laser on choice bias in DMS indirect/direct/noOpsin 						\
 			illumination groups, as well as effects of laser on choice bias in 							\
			NAc indirect/direct/noOpsin groups during AoE 										\
 (4)  ExtData Fig1  - inhibition validation (silicon optrode recording in A2a-Cre/D1R-Cre mice 							\
 			expressing NpHR while ambulating on a running wheel)  									\
			(related to Figure 1a-c, ExtData Fig3)											\
 (5)  ExtData Fig2  - stereological quantification of striatal D1R/D2R fluoroescent in situ hybridization data					\
 			(related to Figure 1a-c, Figure 2, ExtData Fig1, and ExtData Fig3)							\
 (6)  ExtData Fig3  - raster plots of all significantly laser inhibited single-units in A2a- or D1R-cre mice					\
 			(related to Figure 1a-c, ExtData Fig1)											\
 (7)  ExtData Fig4  - additional summary analysis of effects of DMS laser on motor variables during VR corridor					\
 			(related to Figure 1d-j)												\
 (8)  ExtData Fig5  - real-time conditioned place preference, NpHR in DMS of A2a-/D2R- or D1R-Cre mice or NoOpsin laser control			\
 (9)  ExtData Fig6  - additional summary analysis of effects VR decision-making tasks on motor variables					\
  			(related to Figure 2)													\
 (10) ExtData Fig7  - indiv. & group psychometric performance during laser on/off trials across VR tasks and DMS/NAc laser groups		\
   			(related to Figure 3)													\
 (11) ExtData Fig8  - additional comparisons laser-induced choice bias across VR tasks and DMS/NAc laser groups					\
    			(related to Figure 3)													\
 (12) ExtData Fig9  - additional summary analysis of effects of DMS laser on motor variables across VR decision-making tasks			\
  			(related to Figure 2 and Figure 3)											\
 (13) ExtData Fig12 - additional analysis task engagement indicators across GLM-HMM states							\
 			(related to Figure 5-7)													\
 (14) ExtData Fig14 - additional summary analysis of motor variables across GLM-HMM states (laser off and laser on-off seperately)		\
 			(related to Figure 5-7, ExtData Fig9)											\
 (15) ExtData Fig15 - summary of behavioral performance across VR decision-making task shaping							\
 
 see https://github.com/irisstone/glmhmm for plots of:											 	\ 
 (1)  Figure 4      - Bernoulli GLM model, weights, and psychometric performance								\
 (2)  Figure 5      - comparison of GLM vs. GLM-HMM model performance										\
 (3)  Figure 6      - GLM-HMM weights, psychometrics, and descriptive state analyses								\
 (4)  Figure 7      - GLM-HMM state transition analyses												\
 (5)  ExtData Fig10 - GLM-HMM model validation analyses												\
 (6)  ExtData Fig11 - summary of within-session GLM-HMM state occupancy for all individual mice							\
 (7)  ExtData Fig13 - additional model validation and simulation analyses									\
 
 
STEP-BY-STEP									 								\
(1) clone/download this repo and add to matlab path 												\
(2) download data on figshare from link at top 													\
	(a) zipped Data folder (12.8GB), after download and unzip (21.9GB),									\
		***take stock of Data folder location*** 											\
	(b) see <yr local path>/BolkanStoneEtAl/Code/utility/logExplanation.m for extensive							\
		explanation of content of behavior logs											 	\
	(c) see <yr local path>/BolkanStoneEtAl/Code/utility/SpikeDataAll_Explanation.m for 							\
		extensive explanation of contents of spike data summary									 	\	
	(d) Contents of Data folder: 														\
		(i)   DMS_AoE_GLMHMM_states.mat (DMS D1R and D2R concatenated behavioral logs 							\
			with associated GLM-HMM stateIdx and stateProb)										\
		(ii)  OffOn_TasksOrGroup_all.mat (DMS D1R,D2R,NoOpsin and NAc D1R,D2R,NoOpsin 							\
			concatentated behavioral logs in mice performing the evidence 								\
			accumulation (_aoe), no distrators (_nd) or permanent cues (_pc) tasks) 							\
		(iii) rtPPdata.mat (real-time conditioned place preference, Ethovision output 							\
			measures, DMS D1R,D2R,NoOpsin laser) 											\
		(iv)  runningWheelSpikes.mat (summary of 110 single-units in DMS from D1R-Cre/A2a-Cre 						\
			mice expressing DIO-NpHR, laser sweeps while freely moving on running wheel) 						\
		(v)   taskShapingData.mat (summary of AoE, N.D., and P.C. task shaping performance)						\
		(vi)  virtualCorridor.mat (DMS D1R,D2R,NoOpsin concatenated behavioral logs in mice 						\
			navigating a virtual corridor)      									         	\
		(vii) 4 .xls docs summarizing stereological quantification of D1R/D2R in situ 							\
			hybridization expression in A2a/D2R/D1R-Cre mouse lines							          	\	
(3) in MATLAB open/edit .m file <yr local path>/BolkanStoneEtAl/Code/utility/globalParams.m                                                            		\
	(a) change globalParams.dataPath to match location of downloaded data                                                         		\
(4) run e.g. figure1_script in matlab command line                                                                                    		\
	(a) will generate all figure1 data plots and save in your matlab workspace the 								\
		analyzed variables used to generate each plot  											\
(5) otherwise, edit figure1_script from matlab command line                                                                           		\
	(a) each script is bracketed in sections that can be run in steps to                                                          		\
		(i)   load pre-processed source data                                                                                  		\
		(ii)  analyze data in commented sections                                                                              		\
		(iii) plot data of given analysis                                                                                     		\			
