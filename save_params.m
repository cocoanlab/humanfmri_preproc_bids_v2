clear;clc
addpath(genpath('~/Programs/essentials')); 
addpath(genpath('~/github/humanfmri_preproc_bids_v2'));

params.preproc_dir      = '';      % char
params.subjects         = {};      % cell
params.intput_prefix    = '';      % char
params.input_postfix    = '';      % char
params.output_prefix    = '';      % char
params.output_postfix   = '';      % char
params.mask_prefix      = '';      % char
params.mask_postfix     = '';      % char

params.t_disdaq         = NaN; 
params.TR               = NaN;
params.echos            = [NaN NaN]; % in milisec
params.ref_echonum      = NaN;        % middle echo as ref.

params.add_dirs.fsl_dir         = ''; % till /bin
params.add_dirs.dcm2niix_dir    = ''; 
params.add_dirs.tedana_dir      = '';

params.epi_enc_dir    = 'ap';
parmas.ap_indicator   = '*to_ap*';
params.pa_indicator   = '*_64ch_pa_*';
params.t1_indicator = '*T1*';



%% Specify required program directories.
add_dirs.fsl_dir             = 
add_dirs.dcm2niix_dir        = ''; 
add_dirs.tedana_dir          = '/home/cocoan/Programs/fsl/bin'; % 

int_dirs = fieldnames(add_dirs);
for dir_i = 1:numel(int_dirs)
    addstr = add_dirs.(int_dirs{dir_i});
    if (~isempty(addstr)) && (~contains(getenv('PATH'), addstr))
        setenv('PATH', [getenv('PATH'), ':', addstr])
    end
end 