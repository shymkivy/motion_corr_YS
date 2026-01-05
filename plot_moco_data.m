clear;
close all;

pwd2 = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath([pwd2 '\functions']));

%%
fpath = 'F:\VR\data_proc\L\preprocessing\L_10_21_25_h5cutsdata.mat';

data_load = load(fpath);

cuts_data = data_load.cuts_data;
params = data_load.params;

%%
f_mc_plot_cuts_data(cuts_data, params.save_fname);

fprintf('smooth params:\n')
disp(params.params_moco.smooth_std);


figure()
plot(cuts_data{1}.vid_cuts_trace)