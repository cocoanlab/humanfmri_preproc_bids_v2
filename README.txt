<2024.12.23.Mon version>
% PLEASE check the latest update version

1. Programs
 1-1. MATLAB (v.2023b) 
 1-2. FSL (v. 6.0.7.10)

2. Essentials
 2-1. Canlab Private
 2-2. canlab
	1) CanlabCore
	2) Lindquist_Dyanmic_Correlation
	3) MediationToolbox
 2-3. cocoan
	1) cocoanCORE
	2) humanfmri_preproc_bids
	3) surface_preprocessing
 2-4. masksprivate
	1) Masks_private
 2-5. yanchogosu
	1) yanchogosu_toolbox

3. SPM12
 : cocoan setting in "nas01/resources/spm12"
 ( List                      original ver. -> cocoan ver.)
 - defaults.stats.fmri.hpf = 128           ->  180;
 - defaults.stats.fmri.cvi = 'AR(1)'       -> 'none';
 - defaults.mask.thresh    = 0,8           -> -Inf;
 - defaults.stats.maxmem   = 2^30          -> 2^33;
