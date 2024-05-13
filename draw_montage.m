function draw_montage(output_nii, saveto)

spm_image('init', output_nii);
cluster_orthviews_montage(8, 'axial', [], 'onerow');
fh = findobj('Tag', 'Graphics'); set(fh, 'Visible', 'off');
fh = findobj('Tag', 'montage_axial');
figure(fh);
scn_export_papersetup(300);
saveas(fh, saveto);

close all;

end