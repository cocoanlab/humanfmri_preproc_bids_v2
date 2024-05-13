function obj = spk_calc_save(obj)

[foldername_is, filename_is] = fileparts(obj.fullpath(1, :));
subdir        = fileparts(fileparts(foldername_is));
obj.images_per_session = size(obj.dat,2);

% spike identification
diary(fullfile(subdir,'qc_images', ['qc_diary_' filename_is '.txt']));
obj = preprocess(obj, 'outliers', 'plot');  % Spike detect and globals by slice
subplot(5, 1, 5);
obj = preprocess(obj, 'outliers_rmssd', 'plot');  % RMSSD Spike detect
diary off;
sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);

qcspikefilename = fullfile(subdir,'qc_images', ['qc_spike_plot_' filename_is '.png']); 
saveas(gcf,qcspikefilename);

close all;

end