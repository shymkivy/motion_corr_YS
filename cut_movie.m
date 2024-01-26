addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\general_functions');

fpath =  'F:\AC_data\caiman_data_missmatch\movies\';
fname = 'M10_im1_A2_ammn1_5_31_20_cut.h5';

start = 5000;
dur = 1000;

[~, fname2, ext] = fileparts(fname);

Y = h5read([fpath, fname], '/mov');

Y2 = Y(:,:,start:(start+dur));

f_save_mov_YS(Y2, [fpath, fname2, '_', num2str(start), ext], '/mov');


