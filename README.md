# BolkanStoneEtAl
<pre>
 code to analyze data/generate plots in Bolkan, Stone et al (2022) \

(1) clone/download this repo and add to matlab path \
(2) download data from figshare \
	(a) sbolkan/BolkanStoneEtAl  zipped Data folder (12.8GB) \
	(b) after download and unzip (21.9GB) - take stock of Data folder location \
	(c) Contents of Data folder (see BolkanStoneEtAl/Code/utility/logExplanation.m for more extensive explanation of behavior logs) : \
		(i)   DMS_AoE_GLMHMM_states.mat (DMS D1R and D2R concatenated behavioral logs with associated GLM-HMM stateIdx and stateProb)
		(ii)  OffOn_TasksOrGroup_all.mat (DMS D1R,D2R,NoOpsin and NAc D1R,D2R,NoOpsin concatentated behavioral logs \ 
		      in mice performing the evidence accumulation (AoE), no distrators (nd) or permanent cues (pc) tasks
		(iii) rtPPdata.mat (real-time conditioned place preference, Ethovision output measures, DMS D1R,D2R,NoOpsin laser)
		(iv)  runningWheelSpikes.mat (summary of 110 single-units in DMS from D1R-Cre/A2a-Cre mice expressing DIO-NpHR, laser sweeps \ 
		      while freely moving on running wheel)
		(v)   taskShapingData (summary of AoE, N.D., and P.C. performance during task shaping)
		(vi)  virtualCorridor (DMS D1R,D2R,NoOpsin concatenated behavioral logs in mice navigating a virtual corridor)
		(vii) 4 .xls docs summarizing stereological quantification of D1R/D2R in situ hybridization expression \	
(3) in MATLAB open .m file .../BolkanStoneEtAl/Code/utility/globalParams.m \
	(a) change globalParams.dataPath to match location of downloaded data \
(4) run e.g. figure1_script in matlab command line \
	(a) will generate all figure1 data plots and save in your matlab workspace the analyzed variables used to generate each plot \
(5) otherwise, edit figure1_script from matlab command line \
	(a) each script is bracketed in sections that can be run in steps to \
		(i)   load pre-processed source data \
		(ii)  analyze data in commented sections \
		(iii) plot data of given analysis \		
(6) for additional GLM-HMM model validation and analysis/plots in Bolkan,Stone, et al see \
https://github.com/irisstone/glmhmm \
which also includes more general use applications to similarly structured data.


