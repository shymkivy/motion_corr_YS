%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiffs
%           2: Tiff stack
%           3: H5 stack
%       Crop synch pulses
%       bidi shift correct
%       custom moco
%       Save as H5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear;
close all;

%%
pwd2 = fileparts(matlab.desktop.editor.getActiveFilename);
addpath([pwd2 '\functions']);

%%

params.load_dir = 'D:\data\AC\2p\2020\M8_5_21_20_missmatch';  % 
params.save_dir = 'F:\AC_data\caiman_data_missmatch';
%save_dir = 'C:\Users\rylab_dataPC\Desktop\Yuriy\DD_data\proc_data';
%save_dir = 'E:\data\Auditory\caiman_out_multiplane';
%save_dir = 'J:\mouse\backup\2018\caiman_out_dLGN';
%save_dir = 'L:\data\Auditory\caiman_out';

params.load_fname = 'A2_ammn-001';  % can be a dir if is a list of tiffs
params.save_fname = 'M8_im1_A2_ammn1_5_21_20';

params.im_target_fname = ''; % cuts file with target for moco, or empty string for none

params.num_planes = 1; % enter number of planes
params.do_moco = 1;
params.do_bidi = 0;
params.moco_zero_edge = 1;

%%
f_preprocess_mov(params);

%%