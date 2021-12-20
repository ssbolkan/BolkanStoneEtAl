# BolkanStoneEtAl
 code to analyze data/generate plots in Bolkan, Stone et al (2022)

(1) clone/download repo and add to matlab path
	(a) 
(2) download data from figshare 
	(a) 
	(b) 
	(c) (12.8GB), 
	(d) unzip (21.9GB) and take stock of Data folder location 
(3) in MATLAB open .m file .../BolkanStoneEtAl/Code/utility/globalParams.m
	(a) change globalParams.dataPath to match location of downloaded data
(4) run e.g. figure1_script in matlab command line 
	(a) will generate all figure1 data plots and save in your matlab workspace the analyzed variables used to generate each plot
(5) otherwise, edit figure1_script from matlab command line
	(a) each script is bracketed in sections that can be run in steps to
		(i)   load pre-processed source data
		(ii)  analyze data in commented sections
		(iii) plot data of given analysis

