<2024.12.23.Mon version>
% PLEASE check the latest update version

################### FILE LOCATIONS START ####################

 ~/programs/
	/MATLAB 						(v.2023b) 
	/fsl 							(v. 6.0.7.10)

 ~/github/
 	/canlab							(pulled on 24.12.23)
 		/CanlabCore
 		/CanlabPrivate
	/cocoan							(pulled on 24.12.23)
		/cocoanCORE
		/humanfmri_preproc_bids_v2
		/humanfmri_preproc_bids
		/surface_preprocessing
	/spm12							(from ~/nas01/resources/spm12)
	
################### FILE LOCATIONS END ####################



# DESCRIPTIONS

Cocoan specific modifications in "~/github/spm12/spm_defaults.m"
	defaults.stats.fmri.hpf = 128           ->  180;
	defaults.stats.fmri.cvi = 'AR(1)'       -> 'none';
	defaults.mask.thresh    = 0.8           -> -Inf;
	defaults.stats.maxmem   = 2^30          -> 2^33;
	
	

	

