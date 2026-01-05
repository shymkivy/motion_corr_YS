function f_preprocess_mov(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiffs
%           2: Tiff stack
%           3: H5 stack
%       Crop pulses
%       Save as H5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('params', 'var'); params = struct(); end
params = f_set_params(params);

%% load
Y = f_load_mov([params.load_dir '\\' params.load_fname], params);
Y{1} = Y{1}(:,:,1:10000);

%% compute cuts
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

end