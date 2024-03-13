%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiff folde
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

% params.load_dir = 'F:\AC_data\caiman_data_missmatch\movies';  % 
% params.save_dir = 'F:\AC_data\caiman_data_missmatch';
% 
% params.load_fname = 'M10_im1_A2_ammn1_5_31_20_cut_5000.h5';  % can be a dir if is a Prairie list of tiffs
% params.fname = 'M8_im1_A2_ammn1_5_21_20';

params.load_dir = 'C:\Users\ys2605\Desktop\stuff\data_others\motion_correction_hb';  % 
params.save_dir = 'C:\Users\ys2605\Desktop\stuff\data_others\motion_correction_hb';

params.load_fname = 'ch1_m1m-nm_dob11202023_iue111723-spasap5-jrgeco_NA1_920nm150mw_z250_zoom1-512x102_135hz_stim_00003_pl1.h5';  % can be a dir if is a Prairie list of tiffs
params.fname = 'ch1_mc';

params.load_tif_format = 'int16';  % for scanimage or leave empty for auto

params.im_target_fname = ''; % cuts mat file with target for moco, string or cell of strings, or *_h5cutsdata.mat file

params.num_planes = 1; % enter number of planes

params.do_moco = 1;
params.do_bidi = 0;
params.align_pulse_crop_method = 0;

params.moco_rigid_method = 12; % 0=one iteration method; other are described inside f_preprocess_mov.m
params.moco_zero_edge = 0;
%%
f_preprocess_mov(params);

%%