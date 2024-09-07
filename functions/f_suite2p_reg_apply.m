function dreg = f_suite2p_reg_apply(data, dsall)

%% register Y

% suite2p_matlab_path = 'C:\Users\ys2605\Desktop\stuff\libs\Suite2P_matlab';
% addpath(suite2p_matlab_path);
% addpath([suite2p_matlab_path '\preRegistration']);
% addpath([suite2p_matlab_path '\registration']);
% addpath([suite2p_matlab_path '\utils']);

[d1, d2, T] = size(data);

ops.splitFOV = [1 1];
ops.NiterPrealign = 20;
ops.kriging = 1; % subpix align??
ops.useGPU = 1;
ops.planesToProcess = 1;
ops.alignAcrossPlanes = 0;
ops.nplanes = 1;
ops.smooth_time_space = [];


[Ly, Lx, T] = size(data);
ops.Ly = Ly;
ops.Lx = Lx;

[xFOVs, yFOVs] = get_xyFOVs(ops);

ops1 = cell(1, 1);
ops1{1} = ops;

ops1{1,1}.DS          = [];
ops1{1,1}.CorrFrame   = [];
ops1{1,1}.mimg1       = zeros(ops1{1,1}.Ly, ops1{1,1}.Lx);
ops1{1}.Nframes(1) = 0;

data = reshape(data, d1,d2,T,1);

dreg = rigidMovie(data, ops1, dsall, yFOVs, xFOVs);


end
