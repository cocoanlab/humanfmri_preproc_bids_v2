function out = make_param

%% Specify Inputs
out.preproc_dir                                 = '';
% out.project                                   = '';
out.subjects                                    = '';

out.dcm_funcstr                                 = {'HEAT', 'REST'};
out.dcm_fmapstr                                 = {'_distortion_corr_'};
out.dcm_T1str                                   = {'_T1_'};
% dcm_*str          : Specifiers that can indicate "func", "fmap", "T1"
%                     from the scanner. refer to "$dcmscnr_dir"
%                     e.g., If my resting-state run EPIs are stored as
%                                 "series_36_33_ISO_ME_func_REST01"
%                           and task run EPIS as   
%                                 "series_36_33_ISO_ME_func_TASK01"
%                           I can specify "$dcm_funcstr" as {'_TASK', '_REST'}

% Specify Runs
% out.func_run_nums{1}                            = 1:8;
% out.func_tasks{1}                               = 

%% 1. Disdaq
out.('disdaq').t_disdaq                         = 5; % in sec
out.('disdaq').output_prefix                    = 'c_';
out.('disdaq').TR                               = 1; % in sec
out.('disdaq').echos                            = [13 29.3 45.58]; % check your own sequence

%% 2. QC (outlier detection)
% 2-1. Mean Image (QC)
out.('save_mean_img').input_prefix              = 'c_';
out.('save_mean_img').output_prefix             = 'mean_beforepreproc_';
out.('save_mean_img').ref_echonum               = 2; % middle echo as ref.

% 2-2. Save Implicit Mask
out.('save_implicit_mask').input_prefix         = 'c_'; % disdaq ver.
out.('save_implicit_mask').output_postfix       = '_bet';
out.('save_implicit_mask').ref_echonum          = 1;    % Using 1st Echo for BET

% 1-3. Spike Detect Before Preprocessing
out.('spike_id').input_prefix                   = 'c_'; % disdaq ver
out.('spike_id').input_postfix                  = '';
out.('spike_id').mask_prefix                    = '';
out.('spike_id').mask_postfix                   = '_betmask';
out.('spike_id').ref_echonum                    = 2; % Using Echo 2 for SpikeDetection...

%% 2. Slice Time Correction
out.('slice_timing_correction').input_prefix    = 'c_';
out.('slice_timing_correction').output_prefix   = 'a';

% ** ONLY FOR MULTI-ECHO **  %
% The length of TR in our current ME seuqnce is 1 secs. According to
% tedana community, they recommend to use the parameters of one of echos.
% So We acquire the realgnment parameters of Second echo images. 

%% 3. Motion Correction (Realignment)
out.('motion_correction').input_prefix          = 'c_';     % before slice timing correction.
out.('motion_correction').output_prefix         = 'r';
out.('motion_correction').ref_imgnum            = 1;        % not to be used if use_sbref is true.  
out.('motion_correction').use_sbref             = true; 
out.('motion_correction').ref_echonum           = 2;
out.('motion_correction').ref_prefix            = 'ref_';

% 3-2. Actual Realign
% same name? 'motion correction'?
out.('motion_correction').input_prefix          = 'ac_';  
out.('motion_correction').output_prefix         = 'r';
out.('motion_correction').ref_echonum           = 2;

%% 4. TEDANA
out.('tedana').input_prefix                     = 'rac_';
out.('tedana').mask_postfix                     = '_betmask';
out.('tedana').output_prefix                    = 't';

%% 5. Distotion Correction
% Please chekc the FSL directory
out.('distotion_correction').epi_enc_dir        = 'ap';
out.('distotion_correction').input_prefix       = 'trac_';
out.('distotion_correction').output_prefix      = 'dc';
out.('distotion_correction').ap_indicator       = '*to_ap*';
out.('distotion_correction').pa_indicator       = '*_64ch_pa_*';

%% 6. co-registration and normalization
% co-reg: T1 to EPI 
out.('co_registration').t1_indicator            = 'anat*T1*';
out.('co_registration').ref_prefix              = 'mean_dctrac_';
out.('co_registration').input_prefix            = 'dctrac_';
out.('co_registration').output_prefix           = 'w';

%% 7. Smoothing
out.('smoothing').input_prefix                  = 'wdctrac_'; 
out.('smoothing').output_prefix                 = 's'; 
out.('smoothing').fwhm                          = 6; 

end