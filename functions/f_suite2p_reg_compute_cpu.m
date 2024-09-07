function [dsall, input_frame, ops1] = f_suite2p_reg_compute_cpu(data, input_frame, ops)

%% register Y

% suite2p_matlab_path = 'C:\Users\ys2605\Desktop\stuff\libs\Suite2P_matlab';
% addpath(suite2p_matlab_path);
% addpath([suite2p_matlab_path '\preRegistration']);
% addpath([suite2p_matlab_path '\registration']);
% addpath([suite2p_matlab_path '\utils']);

[d1, d2, T] = size(data);

if ~exist('ops', 'var')
    ops = struct;
end

ops.splitFOV = [1 1];
ops.NiterPrealign = 20;
ops.kriging = 1; % subpix align??
ops.useGPU = 0;
ops.planesToProcess = 1;
ops.alignAcrossPlanes = 0;
ops.nplanes = 1;
ops.smooth_time_space = [0 0 0];

[Ly, Lx, T] = size(data);
ops.Ly = Ly;
ops.Lx = Lx;

[xFOVs, yFOVs] = get_xyFOVs(ops);

ops1 = cell(1, 1);
ops1{1} = ops;


if ~exist('input_frame', 'var')
    input_frame = [];
end
if ~isempty(input_frame)
    ops1{1,1}.mimg = input_frame;
else
    % get 500 random frames
    samp_frames = randsample(T, 1000);
    ops1{1,1} = alignIterative(single(data(:,:,samp_frames)), ops);
    input_frame = ops1{1,1}.mimg;
end

ops1{1,1}.DS          = [];
ops1{1,1}.CorrFrame   = [];
ops1{1,1}.Corr_z      = [];
ops1{1,1}.mimg1       = zeros(ops1{1,1}.Ly, ops1{1,1}.Lx);
ops1{1}.Nframes(1) = 0;

data = reshape(data, d1,d2,T,1);

[dsall, ops1] = rigidOffsets_YS(data, 1, 1, 1, ops, ops1);
% [dsall, ops1] = rigidOffsets(data, 1, 1, 1, ops, ops1);

%dreg = rigidMovie(data, ops1, dsall, yFOVs, xFOVs);

% f_save_tif_stack2_YS(Y_prereg, [save_path '\movie_pre_reg']);
% f_save_tif_stack2_YS(dreg, [save_path '\movie_post_reg']);


end
