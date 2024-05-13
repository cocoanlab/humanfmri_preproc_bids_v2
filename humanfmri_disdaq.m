function humanfmri_disdaq()

n_disdaq = round(t_disdaq / TR);
disdaq_prefix   = 'c_';

fprintf('======THROWING FIRST %02d IMAGES OUT=======\n', n_disdaq)
for sub_i = 1:numel(subjects)
    subj_is = subjects{sub_i};
    func_runs = sort_ycgosu(fullfile(preproc_dir, subjects{sub_i}, ...
        'func', '*bold*'));

    for run_i = 1:numel(func_runs)
        [~, taskstr] = fileparts(func_runs{run_i});
        try
            run_niis = sort_ycgosu(fullfile(func_runs{run_i}, 'func*nii'));
        catch
            fprintf('-----%s----%s MISSING------- CONTINUE\n', subj_is, taskstr);
        end

        for echo_i = 1:numel(run_niis)
            input_nii = run_niis{echo_i};
            [foldername_is, filename_is, ext] = fileparts(input_nii);
            n_vol = numel(spm_vol(input_nii));
            output_nii = fullfile(foldername_is, [disdaq_prefix, filename_is, ext]);

            if ~exist(output_nii, 'file')
                system(sprintf('fslroi %s %s %d %d', ...
                    input_nii, output_nii, n_disdaq, n_vol - n_disdaq));
            end
        end
    end
end

fprintf('======THROWING FIRST %02d IMAGES OUT-DONE==\n', n_disdaq)

end
