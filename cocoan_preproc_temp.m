function resparams = cocoan_preproc_temp(int_process, params)


% SET environmennt for FSL (and tedana)
param.add_dirs.fsl_dir             = '/home/cocoan/Programs/fsl/bin'; % till /bin
param.add_dirs.dcm2bids_dir        = ''; 
param.add_dirs.tedana_dir          = '/home/cocoan/Programs/fsl/bin'; % 

int_dirs = fieldnames(param.add_dirs);
for dir_i = 1:numel(int_dirs)
    addstr = param.add_dirs.(int_dirs{dir_i});
    if (~isempty(addstr)) && (~contains(getenv('PATH'), addstr))
        setenv('PATH', [getenv('PATH'), ':', addstr])
    end
end 
setenv('FSLOUTPUTTYPE', 'NIFTI');


% SET directory 
addpath(genpath('~/Programs/essentials'));
preproc_dir    = params.preproc_dir; % char
subjects       = params.subjects;    % cell


% SET parameters
if isempty(params.default_prefix)    % default prefix for specifying EPI images.
                                     % If you have followed the cocoan pipeline,
                                     % it's good to leave empty.
    funcstr    = 'func';
else
    funcstr    = params.default_prefix;
end




if ischar(int_process)
    process_is = int_process;
else
    error('First Input should be character specifying the preprocessing step')
end


%% RUN preproces step 

switch lower(process_is)
    
  % =================================================================== %
    %                                                                     %
    %                    PREPROCESSEING STEP                              %
    %                                                                     %
    % =================================================================== %
    
    case {'disdaq'}  
    % 1) Discard some TRs, 2) write save mean images, 
    % and 3) save implicit mask
    
	% Required params
    %   -params.t_disdaq (in seconds.)
    %   -params.TR   (in seconds.)
    
    % Optional params.
    %   -params.output_prefix (default: 'c_') will be appended at the resulting output
           
    case {'save_implicit_mask'}
        
    case {'spike_id'}
    % Outlier detection (spike) based on the statistics from RMSSD and
    % Mahalanobis distance
    % It requires functions in preprocess/canlab                
    

    case {'slice_timing_correction' , 'st_correction'}
    % It can be ignored if only short TR ( less than 1 secs)   .
    
    % ** ONLY FOR MULTI-ECHO **  %
    % The length of TR in our current ME seuqnce is 1 secs. According to
    % tedana community, they recommend to use the parameters of one of echos.
    % So We acquire the realgnment parameters of Second echo images. 
                            
    case {'motion_correction'}
    % The consensus is to do 1) estimate realigment parameter using
    % before-slice timing correction images and 2) alignment using after-slice
    % timing correction images. 
    % by JJ
        
    case {'distotion_correction'}
            
           
        
    case {'co-registration','co-regi', 'normalizatoin','spatial_normalization'}
	% 1) co-reg: T1 to EPI ( cannot understand what actually this part is
	% doing...) and 2) then normalization
    
    case {'smoothing'} 
    % temporal smoothing
    
    % =================================================================== %
    %                                                                     %
    %                    MULTI-echo -related STEP                         %
    %                                                                     %
    % =================================================================== %    
    
    case {'tedana','multi-echo'} % ** ONLY FOR MULTI-ECHO **  %
    % 1) Writing optimally combined images using three echos by estimating t2*
    % 2) Denoising ME-ICA using optimally combined images    
    
    % =================================================================== %
    %                                                                     %
    %                    QC-related STEP                                  %
    %                                                                     %
    % =================================================================== %    
    case {'save_mean_img'}
        
        
    case {'est_nuisance'}
    %
        
end





