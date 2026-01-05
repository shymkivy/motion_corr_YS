addpath('C:\Users\ys2605\Desktop\stuff\motion_corr_YS\functions');

fpath =  'F:\AC_data\caiman_data_echo\movies\';
fname = 'M4264_im8_A1_4cont8_9_5_24_mpl5_pl3.h5';

start = 20000;
dur = 10000;

[~, fname2, ext] = fileparts(fname);

Y = h5read([fpath, fname], '/mov');

Y2 = Y(:,:,start:(start+dur));

f_save_mov_YS(Y2, [fpath, fname2, '_', num2str(start), ext], '/mov');


