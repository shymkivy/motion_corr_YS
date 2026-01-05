%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiff folder
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
addpath(genpath([pwd2 '\functions']));

%params.suite2P_matlab_path = 'C:\Users\ys2605\Desktop\stuff\libs\Suite2P_matlab';

%%

mouse_tag = 'L';
date_tag = '10_28_25';

%params.load_dir = 'F:\VR\10_28_25\L';  % 
params.load_dir = ['D:\VR\', date_tag, '\', mouse_tag];
params.save_dir = ['F:\VR\data_proc\', mouse_tag];

params.load_fname = [mouse_tag, '-001'];  % can be a dir if is a Prairie list of tiffs
params.save_fname = [mouse_tag, '_', date_tag, '_cut'];

params.im_target_fname = ''; % cuts mat file with target for moco, string or cell of strings, or *_h5cutsdata.mat file

params.align_pulse_crop_method = 1;         % 0=no cuts; 1=auto; 2=manual
params.do_moco = 1;
params.do_bidi = 0;

params.manually_split_planes = 0;       % if file has multiple sequential planes not specified by prairie view, indicate number of planes

params.moco_rigid_method = 27; % 0=one iteration method; other are described inside f_preprocess_mov.m
params.moco_zero_edge = 0;      
%%
params = f_set_params(params);

%% load
Y = f_load_mov([params.load_dir '\\' params.load_fname], params);
%Y{1} = Y{1}(:,:,1:10000);

%% compute cuts (cut out alignment pulses)
if params.align_pulse_crop_method
    [Y, cuts_data, params] = f_compute_alignment_pulse_cuts(Y, params);
end

%% bidi fix
if params.do_bidi
    [Y, cuts_data, params] = f_compute_bidi_wrap(Y, cuts_data, params);
end

%% moco
if params.do_moco
    [Y, cuts_data, params] = f_compute_moco_wrap(Y, cuts_data, params);
end

%% save
f_save_proc_outputs(Y, cuts_data,  params);
disp('Done');
