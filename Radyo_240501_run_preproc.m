%% You can start from scratch at this code.
% Things to be checked.
% 1) At which level should be compiled as a function.
% 2) Work with BIDs?
% 3) assert(numel(dcm_funcdirs) == numel(raw_funcdirs), ...
%        'Check desired functional runs and runs from actual scanners');
% 4) Motage code...need further revision. where is "canlab_preproc_show_montage".?
clear;clc
addpath(genpath('~/Programs/essentials')); 
% contains SPM, cocoan, canlab
addpath(genpath('~/github/humanfmri_preproc_bids_v2'));

%% Specify required program directories.
add_dirs.fsl_dir             = '/home/cocoan/Programs/fsl/bin'; % till /bin
add_dirs.dcm2niix_dir        = ''; 
add_dirs.tedana_dir          = '/home/cocoan/Programs/fsl/bin'; % 

int_dirs = fieldnames(add_dirs);
for dir_i = 1:numel(int_dirs)
    addstr = add_dirs.(int_dirs{dir_i});
    if (~isempty(addstr)) && (~contains(getenv('PATH'), addstr))
        setenv('PATH', [getenv('PATH'), ':', addstr])
    end
end 

%% Specify inputs. 1/2
% datdir            : char. where your processed outcome will be stored.
% prjname           : char. your project name. Every output will be in "$datdir/$prjname"
% subjects          : cell (1 x n_subs). 
%                     THIS WILL BE YOUR FUTURE SUBJECT-LEVEL FOLDER NAME.
%                       e.g.,{"sub-003", "sub-007"}                      
datdir              = '/home/cocoan/nas02/data'; 
prjname             = 'Radyo';
subjects            = {'sub-radyo-004'};

%% Specify inputs. 2/2
% dcm_*str          : Specifiers that can indicate "func", "fmap", "T1"
%                     from the scanner. refer to "$dcmscnr_dir"
%                     e.g., If my resting-state run EPIs are stored as
%                                 "series_36_33_ISO_ME_func_REST01"
%                           and task run EPIS as   
%                                 "series_36_33_ISO_ME_func_TASK01"
%                           I can specify "$dcm_funcstr" as {'_TASK', '_REST'}
dcm_funcstr         = {'HEAT', 'REST'};
dcm_fmapstr         = {'_distortion_corr_'};
dcm_T1str           = {'_T1_'};

%% Setting Directory 
% Should manually move dicom files to "$dcmscnr_dir"
datadir             = fullfile(datdir, prjname); 
data_imaging_dir    = fullfile(datadir, 'Imaging');
preproc_dir         = fullfile(data_imaging_dir, 'preprocessed');
dcmscnr_dir         = fullfile(data_imaging_dir, 'dicom_from_scanner');
raw_dir             = fullfile(data_imaging_dir, 'raw');

if ~exist(dcmscnr_dir, 'dir')
    cellfun(@mkdir, fullfile(dcmscnr_dir, subjects))
    fprintf("\n\nMove your original dicom files to \n%s\n\n", dcmscnr_dir);
    fprintf("It should be look like...\n");
    fprintf("~/dicom_from_scanner/%s/%s\n", subjects{1}, 'series_28_33_$SCAN_NAME1')
end

%% Specify Runs & Make Corresponding directories.
clear func_run_num func_tasks
func_run_nums{1}     = 1:8;      func_tasks{1} = [repmat({'heat'}, 1, 7), {'rest'}];
% func_run_nums{2}     = [1:7, 1];      func_tasks{2} = [repmat({'heat'}, 1, 7), {'rest'}];
% one cell per subject.
%   e.g.)
%   func_run_nums{1} = 1:2               => subjects{1} has 2 runs.
%   func_tasks{1}    = {'Heat', 'Caps'}  => subjects{1} has 'heat' and 'caps' runs.

humanfmri_a1_make_directories(subjects, data_imaging_dir, func_run_nums, func_tasks);
% This results in 'raw' folder under "$data_imaging_dir".


%% Move "dicom_from_scanner" to "~/raw/dicom/"
%  This is to preserve dicom files in structured form.
%  moving dicom files from "dicom_from_scanner" to "~/raw/dicom/".

for sub_i = 1:numel(subjects)
    int_sub = subjects{sub_i};
    % set dirs in "dicom_from_scanner"
    dcm_dirs       = sort_ycgosu(fullfile(dcmscnr_dir, int_sub));
    dcm_funcdirs   = dcm_dirs(contains(dcm_dirs, dcm_funcstr));
    dcm_anatdirs   = dcm_dirs(contains(dcm_dirs, dcm_T1str));
    dcm_fmapdirs   = dcm_dirs(contains(dcm_dirs, dcm_fmapstr));

    % set dirs in "~/raw/dicom"
    raw_funcdirs   = sort_ycgosu(fullfile(raw_dir, int_sub, 'dicom', '*task*'));
    raw_anatdir    = sort_ycgosu(fullfile(raw_dir, int_sub, 'dicom', '*anat*'), 'char'); % 
    raw_fmapdir    = sort_ycgosu(fullfile(raw_dir, int_sub, 'dicom', '*fmap*'), 'char'); % 

    assert(numel(dcm_funcdirs) == numel(raw_funcdirs), ...
        'Check desired functional runs and runs from actual scanners');

    fprintf('====== COPYING functional images of %s to ~/raw/dicom  ======\n\n', int_sub)
    for func_i = 1:numel(dcm_funcdirs)
        [~, scan_runname] = fileparts(dcm_funcdirs{func_i});
        copyfile(dcm_funcdirs{func_i}, fullfile(raw_funcdirs{func_i}, scan_runname));
    end

    fprintf('====== COPYING anatomical images of %s to ~/raw/dicom ======\n\n', int_sub)
    for anat_i = 1:numel(dcm_anatdirs)
        [~, scan_runname] = fileparts(dcm_anatdirs{anat_i});
        copyfile(dcm_anatdirs{anat_i}, fullfile(raw_anatdir, scan_runname));
    end

    fprintf('====== COPYING fieldmap images of %s to ~/raw/dicom ======\n\n', int_sub)
    for fmap_i = 1:numel(dcm_fmapdirs)
        [~, scan_runname] = fileparts(dcm_fmapdirs{fmap_i});
        copyfile(dcm_fmapdirs{fmap_i}, fullfile(raw_fmapdir, scan_runname));
    end
end

%% DICOM 2 Nifti using dcm2niix 
% see the examples: https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#General_Usage
% ".nii"
int_modalities = {'func', 'anat', 'fmap'}; % 'dwi';

for sub_i = 1:numel(subjects)
     subj_is        = subjects{sub_i};
     preproc_subdir = fullfile(preproc_dir, subj_is);

     for mod_i = 1:numel(int_modalities)
         int_modstr = int_modalities{mod_i};

         inputdirs  = sort_ycgosu(fullfile(data_imaging_dir, 'raw' , subj_is, ...
             'dicom', ['*', int_modstr, '*']));
         modoutdir  = fullfile(preproc_subdir, int_modstr);

         for dcm_i = 1:numel(inputdirs)
             inputdir = inputdirs{dcm_i};
             [~, taskmod_indicator] = fileparts(inputdir);
             if contains(modoutdir, taskmod_indicator)
                outputdir = modoutdir;
             else
                outputdir = fullfile(modoutdir, taskmod_indicator);
             end

             if ~exist(outputdir, 'dir'), mkdir(outputdir); end
             system(sprintf('dcm2niix -w 0 -o %s %s', outputdir, inputdir))
         end
     end
end
fprintf('DONE dcm2niix\n')

%% START MAIN PREPROCESSING...
t_disdaq = 7;   % in sec
TR       = 0.83; % in sec
echos    = [13 29.3 45.58]; %

dirnames2make = {'mean_func', 'qc_images', 'covariates'};
for sub_i = 1:numel(subjects)
    for dir_i = 1:numel(dirnames2make)
        dir2make = fullfile(preproc_dir, subjects{sub_i}, dirnames2make{dir_i});
        if ~exist(dir2make, 'dir'), mkdir(dir2make); end
    end
end

%% 1. Disdaq. (prefix: c_)
int_process = 'disdaq';
setenv('FSLOUTPUTTYPE', 'NIFTI');
n_disdaq = ceil(t_disdaq / TR);
output_prefix   = 'c_';

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis = sort_ycgosu(fullfile(func_runs{run_i}, 'func*nii'));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        for e_i = 1:numel(run_niis)
            input_nii = run_niis{e_i};
            [folder_is, file_is, ext] = fileparts(input_nii);
            n_vol = numel(spm_vol(input_nii));
            output_nii = fullfile(folder_is, [output_prefix, file_is, ext]);

            if ~exist(output_nii, 'file')
                system(sprintf('fslroi %s %s %d %d', ...
                    input_nii, output_nii, n_disdaq, n_vol - n_disdaq));
            end
        end
    end
end
                        


%% 1 - 1. Write mean image and save mean image (QC)
setenv('FSLOUTPUTTYPE', 'NIFTI');
input_prefix    = 'c_';
output_prefix   = 'mean_beforepreproc_';
ref_echonum     = 2; % middle echo as ref.

for sub_i = 1:numel(subjects)
    subj_is      = subjects{sub_i};
    mean_funcdir = fullfile(preproc_dir, subj_is, 'mean_func');
    qc_dir       = fullfile(preproc_dir, subj_is, 'qc_images');

    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
        end

        input_nii  = run_niis{contains(run_niis, sprintf('_e%d', ref_echonum))};
        [~, file_is, ext] = fileparts(input_nii);
        output_nii = fullfile(mean_funcdir, [output_prefix, file_is ext]);
        saveto     =  fullfile(qc_dir, [output_prefix, file_is, '.png']);
        cmd        = sprintf('fslmaths %s -Tmean %s', input_nii, output_nii);
        if ~exist(output_nii,'file') && ~exist(saveto, 'file')
            system(cmd)
            draw_montage(output_nii, saveto) 
        end
    end
end

%% 1 - 2. save implicit mask. 
setenv('FSLOUTPUTTYPE', 'NIFTI');
input_prefix    = '';
output_postfix  = '_bet';
ref_echonum     = 1; % Using Echo 1 for BET.

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        input_nii  = run_niis{contains(run_niis, sprintf('_e%d', ref_echonum))};
        [folder_is, file_is, ext] = fileparts(input_nii);
        output_nii = fullfile(folder_is, [file_is, output_postfix, ext]);
        if ~exist(output_nii, 'file')
            system(sprintf('bet %s %s -f 0.3 -n -R', input_nii, output_nii)); 
        end
        
        mask_nii   = fullfile(folder_is, [file_is, output_postfix, 'mask', ext]);
        if ~exist(mask_nii, 'file')
            system(sprintf('fslmaths %s -bin %s',    output_nii, mask_nii)); 
        end        
    end
end

%% 1 - 3. SpikeDetect Before Preproc.
input_prefix    = 'c_';  input_postfix   = '';
mask_prefix     = '';    mask_postfix    = '_betmask';
ref_echonum     = 2; % Using Echo 2 for SpikeDetection...

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end
        input_nii  = run_niis{contains(run_niis, sprintf('_e%d.nii', ref_echonum))};
        spk_mat    = fullfile(preproc_dir, subj_is,'covariates', sprintf('spike_covariates_%s.mat', taskstr));

        if ~exist(spk_mat, 'file')
            [folder_is, file_is, ext] = fileparts(input_nii);
            mask_nii   = sort_ycgosu(fullfile(folder_is, ['*', mask_postfix, ext]), 'char');
            dat = fmri_data(input_nii, mask_nii); 
            dat = spk_calc_save(dat);
            
            spike_covariates.dat   = dat.covariates;
            spike_covariates.files = input_nii;
           
            save(spk_mat,'spike_covariates');
        end
    end
end

%% 2. Slice time correction 
% The length of TR in our current ME seuqnce is 1 secs. According to
% tedana community, they recommend to use the parameters of one of echos.
% So We acquire the realgnment parameters of Second echo images 
input_prefix    = 'c_';
output_prefix   = 'a'; 

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
            run_jsons = sort_ycgosu(fullfile(func_runs{run_i}, '*json'));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end
        
        for e_i = 1:numel(run_niis)
            input_nii = run_niis{e_i};
            [folder_is, file_is, ext] = fileparts(input_nii);
            output_nii = fullfile(folder_is, [output_prefix, file_is, ext]);
            if exist(output_nii, 'file'), continue; end

            input_json = fullfile(run_jsons{e_i});
            fid = fopen(input_json);raw = fread(fid, inf);str = char(raw');
            fclose(fid);json_read = jsondecode(str);

            slice_time = json_read.SliceTiming .* 1000; % sec -> msec
            mbf = json_read.MultibandAccelerationFactor;
            
            slice_timing_job = []; 
            slice_timing_job{1}.spm.temporal.st.scans{1} = spm_select('expand', {input_nii});
            % 1. nslices
            Vfist_vol = spm_vol([input_nii, ',1']);
            numSlices = Vfist_vol(1).dim(3);
            slice_timing_job{1}.spm.temporal.st.nslices = numSlices;
    
            % 2.TR
            slice_timing_job{1}.spm.temporal.st.tr = TR;
    
            % 3. Acqui time
            slice_timing_job{1}.spm.temporal.st.ta = TR - TR*mbf / numSlices;
    
            % 4. slice order
            slice_timing_job{1}.spm.temporal.st.so = slice_time;
            slice_timing_job{1}.spm.temporal.st.refslice = min(slice_time);
            slice_timing_job{1}.spm.temporal.st.prefix = output_prefix;
            % saving slice time correction job
            
            % RUN 
            spm('defaults', 'fmri');    spm_jobman('initcfg');
            spm_jobman('run', slice_timing_job)
        end
    end
end
fprintf('======RUNNING ST CORRECTION-DONE=====\n');
%% 3-1. Motion correction (realignment)
% The consensus is to do 1) estimate realigment parameter using
% before-slice timing correction images and 2) alignment using after-slice
% timing correction images.
%
% Bit confusing... using 1) param for 2) or using param seperately in 2) ?
% For this code, "using 1) param for 2)"

% 1) Estimate Params
setenv('FSLOUTPUTTYPE', 'NIFTI');
input_prefix       = 'c_';  % before slice timing correction.
output_prefix      = 'r';
ref_imgnum         = 1; % not be used if use_sbref is true.  
use_sbref          = true; 
ref_echonum        = 2;
ref_prefix         = 'ref_';

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        input_nii = run_niis{contains(run_niis, sprintf('_e%d.nii', ref_echonum))};
        [folder_is, file_is, ext] = fileparts(input_nii);
        ref_nii = fullfile(folder_is, [ref_prefix file_is ext]);
        if ~exist(ref_nii, 'file')
            if use_sbref
                sbref_nii = sort_ycgosu(fullfile(strrep(folder_is, '_bold', '_sbref'), ...
                    sprintf('*_e%d.nii', ref_echonum)), 'char');
                copyfile(sbref_nii, ref_nii);
            else
                system(sprintf('fslroi %s %s %d %d', input_nii, ref_nii, ref_imgnum-1, 1))
            end
        end
        
        output_nii = fullfile(folder_is, [output_prefix, file_is, ext]);
        if ~exist(output_nii, 'file')
            system(sprintf('mcflirt -in %s -reffile %s -o %s -plots -mats', ...
                input_nii, ref_nii, output_nii));
    
            mvmts = importdata([output_nii, '.par']);
            saveto = fullfile(preproc_dir, subj_is,'qc_images', ['qc_mvmt_' taskstr '.png']); 
            mvmt_calc_save(mvmts, saveto, 'software', 'FSL');        
        end
    end
end


%% 3-2. Actual Realign 
input_prefix       = 'ac_';  
output_prefix      = 'r';
ref_echonum        = 2; % 

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));
    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
        end

        for e_i = 1:numel(run_niis)
            input_nii = run_niis{e_i};
            [folder_is, file_is, ext] = fileparts(input_nii);
            output_nii = fullfile(folder_is, [output_prefix, file_is, ext]);
            if ~exist(output_nii, 'file')
                [~, tempdir] = system('mktemp -d');
                tempdir = strtrim(tempdir);
                system(sprintf('fslsplit %s %s -t', input_nii, fullfile(tempdir, 'images')));
                temp_inlist  = filenames(fullfile(tempdir, 'images*.nii'));
                temp_outlist = strrep(temp_inlist, fullfile(tempdir, 'images'), fullfile(tempdir, 'rimages'));
                mc_mats      = sort_ycgosu(fullfile(folder_is, '*mat', 'MAT*'));
                for vol_i = 1:numel(temp_inlist)
                    system(sprintf('flirt -in %s -ref %s -applyxfm -init %s -out %s', ...
                        temp_inlist{vol_i}, ref_nii, mc_mats{vol_i}, temp_outlist{vol_i}));
                end
                temp_outlist_cat = strcat(temp_outlist, {' '});
                temp_outlist_cat = cat(2, temp_outlist_cat{:});
                system(sprintf('fslmerge -t %s %s', output_nii, temp_outlist_cat));
                pause(1);
                system(sprintf('rm -r %s', tempdir));
            end
        end
    end
end

%% 4. RUN TEDANA 
% 1) Writing optimally combined images using three echos by estimating t2*
% 2) Denoising ME-ICA using optimally combined images
input_prefix   = 'rac_';
mask_postfix   = '_betmask';
output_prefix  = 't';

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end
        
        [folder_is, file_is] = fileparts(run_niis); 
        folder_is = unique(folder_is);
        file_is   = unique(cellfun(@(x) x(1:end-3), file_is, 'UniformOutput', false));
        folder_is = folder_is{:}; file_is = file_is{:}; % 
        
        mask_nii  = sort_ycgosu(fullfile(folder_is, ['*' mask_postfix '.nii']), 'char');
        assert(numel(run_niis) == numel(echos), 'Check your inputs and number of echos')
        n_echos = numel(run_niis);
        tedana_outdir = fullfile(folder_is, 'tedana_output');

        cmdstr1 = repmat('%s ', 1, n_echos);
        cmdstr2 = repmat('%.2f ', 1, n_echos);
        cmdstr  = sprintf(['tedana -d ', cmdstr1, '-e ' cmdstr2, '--out-dir %s --mask %s' ...
            ' --fittype curvefit --maxrestart 10 --maxit 500 --tedpca kic'], ...
            run_niis{:}, echos, tedana_outdir, mask_nii);

        output_nii = fullfile(folder_is, [output_prefix, file_is, ext]);
        if ~exist(output_nii, 'file')
            print_header('Running tedana ', [subj_is, ' ', taskstr]);
            system(cmdstr)
            gzoutput_nii = fullfile(tedana_outdir, 'desc-denoised_bold.nii.gz');
            % refer to "desc-tedana_registry.json" for the file you want.
            gunzip(gzoutput_nii)
            copyfile(fullfile(tedana_outdir, 'desc-denoised_bold.nii'), output_nii);
        end   
    end
end

%% 5. Distotion correction 
setenv('FSLOUTPUTTYPE', 'NIFTI');
epi_enc_dir    = 'ap';
input_prefix   = 'trac_';
output_prefix  = 'dc';
ap_indicator   = '*to_ap*';
pa_indicator   = '*_64ch_pa_*';

current_paths  = split(getenv('PATH'), ':');
fsl_dir        = current_paths{contains(current_paths, 'fsl')};

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    fmap_dir = fullfile(preproc_dir, subj_is, 'fmap');
    qc_dir   = fullfile(preproc_dir, subj_is, 'qc_images');
    fmap_combdir = fullfile(fmap_dir, 'distortion_combined');
    if ~exist(fmap_combdir, 'dir'), mkdir(fmap_combdir); end
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));
    ap_nii = sort_ycgosu(fullfile(fmap_dir, [ap_indicator 'nii']), 'char');
    pa_nii = sort_ycgosu(fullfile(fmap_dir, [pa_indicator 'nii']), 'char');
    ap_json = sort_ycgosu(fullfile(fmap_dir, [ap_indicator 'json']), 'char');
    pa_json = sort_ycgosu(fullfile(fmap_dir, [pa_indicator 'json']), 'char');
    dc_param = fullfile(fmap_dir, ['dc_param_', epi_enc_dir, '.txt']);

    topup_out      = fullfile(fmap_combdir, 'topup_out');
    topup_fieldout = fullfile(fmap_combdir, 'topup_fieldout');
    topup_unwarped = fullfile(fmap_combdir, 'topup_unwarped');
    
    if ~exist([topup_unwarped, '.nii'], 'file')
        % calculating field map START...
        dc_comb_nii = fullfile(fmap_combdir, 'dc_combined.nii');
        system(['fslmerge -t ', dc_comb_nii, ' ', ap_nii, ' ', pa_nii]);
        n_ap = numel(spm_vol(ap_nii));  n_pa = numel(spm_vol(pa_nii));    
        distort_json = {ap_json, pa_json};
        rdtime       = NaN(1,2);
        for i = 1:2
            fid = fopen(distort_json{i});   raw = fread(fid, inf);  str = char(raw');
            fclose(fid);  json_read = jsondecode(str);
            rdtime(1,i) = json_read.TotalReadoutTime;
        end
    
        ap_rdtime = rdtime(1);      pa_rdtime = rdtime(2);
        fileID    = fopen(dc_param, 'w');
        dc_param_dat = [repmat([0 -1 0 ap_rdtime], n_ap, 1); repmat([0 1 0 pa_rdtime], n_pa, 1)]; % in case of 'AP'
        fprintf(fileID, repmat([repmat('%.4f\t', 1, size(dc_param_dat, 2)), '\n'], 1, size(dc_param_dat, 1)), dc_param_dat');
        fclose(fileID);
    
        topup_config   = fullfile(fileparts(fsl_dir), 'src/fsl-topup/flirtsch/b02b0.cnf');
        system(['topup --imain=', dc_comb_nii, ' --datain=', dc_param, ' --config=', topup_config, ' --out=', topup_out, ...
            ' --fout=', topup_fieldout, ' --iout=', topup_unwarped]);
        % calculating field map DONE...
        % Save images of field map - NOT YET... where is "canlab_preproc_show_montage"
        % tu_unwarped_png{1} = fullfile(preproc_dir, subj_is,'qc_images', 'topup_unwarped_dir-ap_epi.png');
        % tu_unwarped_png{2} = fullfile(preproc_dir, subj_is,'qc_images', 'topup_unwarped_dir-pa_epi.png');
        %
        % for top_i = 1:numel(tu_unwarped_png)
        %     tu_before_list = cellstr(strcat(dc_comb_nii, ',', num2str([2*top_i-1;2*top_i])));
        %     tu_after_list = cellstr(strcat([topup_unwarped '.nii'], ',', num2str([2*top_i-1;2*top_i])));
        %     draw_montage([tu_before_list; tu_after_list], tu_unwarped_png{top_i});
        %     drawnow;
        % end
        % close all;
    end
    

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        for ii = 1:numel(run_niis)
            input_nii = run_niis{ii};
            [folder_is, file_is, ext] = fileparts(input_nii);
            output_nii     = fullfile(folder_is, [output_prefix, file_is, ext]);
            meanoutput_nii = fullfile(folder_is, ['mean_', output_prefix, file_is, ext]);
            meanoutput_png = fullfile(qc_dir, ['mean_', output_prefix, file_is, '.png']);

            cmd = sprintf(['applytopup --imain="%s" --inindex=1 --topup="%s" --datain="%s"' ...
                ' --method=jac --interp=spline --out="%s"'],  ...
                input_nii, topup_out, dc_param, output_nii);
            if ~exist(output_nii, 'file') || ~exist(meanoutput_nii, 'file')
                system(cmd);
                system(sprintf('fslmaths %s -abs %s', output_nii, output_nii)); % removing negative values...?
                system(sprintf('fslmaths %s -Tmean %s', output_nii, meanoutput_nii));
                draw_montage(meanoutput_nii, meanoutput_png);
            end
        end
    end
end


%% 8. co-registration and normalization
% co-reg: T1 to EPI 
t1_indicator = 'anat*T1*';
ref_prefix = 'mean_dctrac_'; %
input_prefix = 'dctrac_'; %
output_prefix = 'w';

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));
    anat_dir = fullfile(preproc_dir, subj_is, 'anat');
    anat_nii = sort_ycgosu(fullfile(anat_dir, [t1_indicator 'nii']), 'char');

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
            ref_nii   = sort_ycgosu(fullfile(func_runs{run_i}, [ref_prefix, 'func*nii']), 'char');
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        input_nii = run_niis{:};
        [folder_is, file_is, ext] = fileparts(input_nii);
        output_nii     = fullfile(folder_is, [output_prefix, file_is, ext]);
        if ~exist(output_nii, 'file')
            % coreg
            matlabbatch = [];
            print_header('co-registeration', [subj_is, '    ', taskstr]);
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ref_nii ',1']};
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {anat_nii};
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'coreg-';
            spm('defaults','fmri');
            spm_jobman('initcfg');
            spm_jobman('run', matlabbatch);

            % segmentation & normalization
            matlabbatch = [];
            coreg_anat_nii = sort_ycgosu(fullfile(anat_dir, 'coreg-*nii'), 'char');
            print_header('normalization', [subj_is, '    ', taskstr]);
            for j = 1:6
                matlabbatch{1}.spm.spatial.preproc.tissue(j).tpm{1} = [which('TPM.nii') ',' num2str(j)];
            end
            matlabbatch{1}.spm.spatial.preproc.channel.vols{1} = coreg_anat_nii;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'subj';
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
    
            [~, anat_filename] = fileparts(matlabbatch{1}.spm.spatial.preproc.channel.vols{1});
            deformation_nii = fullfile(anat_dir, ['y_' anat_filename '.nii']);
            matlabbatch{2}.spm.spatial.normalise.write.subj.def = {deformation_nii};
            matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {input_nii};
            matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 7;
            matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = output_prefix;
    
            spm('defaults','fmri');
            spm_jobman('initcfg');
            spm_jobman('run', matlabbatch);

            close all;
        end
    end
end

%% 9. Smoothing
input_prefix   = 'wdctrac_'; %
output_prefix  = 's'; %
fwhm = 6; 

for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subj_is, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis  = sort_ycgosu(fullfile(func_runs{run_i}, [input_prefix, 'func*nii']));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
            continue
        end

        input_nii = run_niis{:};
        [folder_is, file_is, ext] = fileparts(input_nii);
        output_nii     = fullfile(folder_is, [output_prefix, file_is, ext]);

        if ~exist(output_nii, 'file')
            matlabbatch = {};
            matlabbatch{1}.spm.spatial.smooth.prefix = output_prefix;
            matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type; 0 = same as before
            matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit mask; 0 = no
            matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(fwhm, 1, 3); % override whatever the defaults were with this
            matlabbatch{1}.spm.spatial.smooth.data = {input_nii};
    
            spm('defaults','fmri');
            spm_jobman('initcfg');
            spm_jobman('run', matlabbatch);
        end
    end
end





