clear;
close all;
pwd2 = fileparts(matlab.desktop.editor.getActiveFilename);
addpath([pwd2 '\functions']);


params.load_dir = 'C:\Users\ys2605\Desktop\stuff\data_others\motion_correction_hb';

params.load_fname = 'ch1_m1m-nm_dob11202023_iue111723-spasap5-jrgeco_NA1_920nm150mw_z250_zoom1-512x102_135hz_stim_00003.tif';

params.num_planes = 1;

max_frames = 1e10;

params.load_tif_format = 'int16';  % for scanimage or leave empty for auto

[Y, params] = f_load_mov(params);

[~, save_fname, ~] = fileparts(params.load_fname);

for n_pl = 1:numel(Y)
    [d1, d2, T] = size(Y{n_pl});
    f_save_mov_YS(Y{n_pl}(:,:,1:min([T, max_frames])), sprintf('%s\\%s_pl%d.h5', params.load_dir, save_fname, n_pl), '/mov');
end

