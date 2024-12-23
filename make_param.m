function out = make_param

% Specify Inputs
out.datdir          = '';
out.project         = '';
out.subjects        = '';

out.dcm_funcstr
out.dcm_fmapstr
out.T1str

% Specify Runss
% Disdaq
out.t_disdaq        = '';
out.output_prefix   = 'c_';
out.TR
out.echos

% 1-1. Mean Image (QC)
out.input_prefix     = 'c_'
out.output_prefix    = 'mean_beforepreproc_';
out.ref_echonum      = 2; % middle echo as ref.

% 1-2. Save Implicit Mask
out.input_prefix     = 'c_'; % disdaq ver.
out.output_postfix   = '_bet';
out.ref_echonum      = 1;    % Using 1st Echo for BET

% 1-3. Spike Detect Before Preprocessing
out.input_prefix     = 'c_'; % disdaq ver
out.input_postfix    = '';
out.mask_prefix      = '';
out.mask_postfix     = '_betmask';
out.ref_echonum      = 2; % Using Echo 2 for SpikeDetection...

% 2. Slice Time Correction
out.input_prefix
out.output_prefix

% 3. Motion Correction (Realignment)
out.input_prefix       = 'c_';  % before slice timing correction.
out.output_prefix      = 'r';
out.ref_imgnum         = 1; % not be used if use_sbref is true.  
out.use_sbref          = true; 
out.ref_echonum        = 2;
out.ref_prefix         = 'ref_';

%% 3-2. Actual Realign 
input_prefix       = 'ac_';  
output_prefix      = 'r';
ref_echonum        = 2; % 
out.('save_implicit_mask').use_sbref = true;
out.ref_imgnum = 1;
out.ref_echonum = 2;

out.slice_timing
out.distortion 
out.coregistration
out.normalization
out.smoothing_kernel


end