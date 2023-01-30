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

%% error check
if ~isfield(params, 'load_dir'); error('Input load dir'); end
if ~numel(params.load_dir); error('Input load dir'); end
if ~isfield(params, 'save_dir'); error('Input save dir'); end
if ~numel(params.save_dir); error('Input save dir'); end

%% default params
% loading
if ~isfield(params, 'num_planes'); params.num_planes = 1; end
if ~isfield(params, 'use_prairie_mpl_tags'); params.use_prairie_mpl_tags = 1; end
if ~isfield(params, 'mpl_tags'); params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; end % multiplane data tags in prairie
if ~isfield(params, 'prairie_chan_tag'); params.prairie_chan_tag = 'Ch2'; end
if ~isfield(params, 'h5_movie_tag'); params.h5_movie_tag = '/mov'; end

% saving
if ~isfield(params, 'save_all_steps'); params.save_all_steps = 0; end                       % during save
if ~isfield(params, 'save_indiv_h5info'); params.save_indiv_h5info = 1; end                 % during save
if ~isfield(params, 'trim_output_num_frames'); params.trim_output_num_frames = 0; end       % 0 or number of frames to save in separate output
if ~isfield(params, 'moco_rigid_method'); params.moco_rigid_method = 1; end                 % which moco rigid method to use
if ~isfield(params, 'moco_nonrigid_method'); params.moco_nonrigid_method = 1; end           % which moco nonrigid method to use
if ~isfield(params, 'overwrite_moco_rigid'); params.overwrite_moco_rigid = 1; end           % overwrite if file already exists
if ~isfield(params, 'overwrite_moco_nonrigid'); params.overwrite_moco_nonrigid = 1; end     % overwrite if file already exists

% processing
if ~isfield(params, 'align_pulse_crop_method'); params.align_pulse_crop_method = 0; end     % crop light synch pulses from movie 0=nothing; 1=auto; 2=manual
if ~isfield(params, 'do_bidi'); params.do_bidi = 0; end
if ~isfield(params, 'do_moco'); params.do_moco = 1; end
if ~isfield(params, 'moco_zero_edge'); params.moco_zero_edge = 1; end
if ~isfield(params, 'do_nonrigid'); params.do_nonrigid = 0; end
if ~isfield(params, 'im_target_fname'); params.im_target_fname = ''; end
if ~isfield(params, 'block_size'); params.block_size = 1000; end

save_dir = params.save_dir;
save_fname = params.save_fname;
num_planes = params.num_planes;
save_all_steps = params.save_all_steps;
save_indiv_h5info = params.save_indiv_h5info;
do_moco = params.do_moco;
do_bidi = params.do_bidi;

%% bidirectional shift fix params
if ~isfield(params, 'params_bidi'); params.params_bidi = struct(); end
params_bidi = params.params_bidi;

if ~isfield(params_bidi, 'smooth_std'); params_bidi.smooth_std = [0.5 0.5 1]; end       % smoothing video before bidi [m, n, T]
if ~isfield(params_bidi, 'reg_lambda'); params_bidi.reg_lambda = [2 0.5]; end           % regularizations for step 1 and 2
if ~isfield(params_bidi, 'fix_range'); params_bidi.fix_range = -15:15; end              % maximal movement
if ~isfield(params_bidi, 'num_iterations'); params_bidi.num_iterations = 1; end
if ~isfield(params_bidi, 'use_planes'); params_bidi.use_planes = [1 3]; end             % for multiplane, whick planes to use
if ~isfield(params_bidi, 'plot_stuff'); params_bidi.plot_stuff = 0; end

%% moco params
if ~isfield(params, 'params_moco'); params.params_moco = struct(); end 
params_moco = params.params_moco;

if ~isfield(params_moco, 'high_val_cut_thresh'); params_moco.high_val_cut_thresh = 0.01; end    % to reduce super bright signals
if ~isfield(params_moco, 'reg_lambda_base'); params_moco.reg_lambda_base = [1 .2]; end          % some regularization
if ~isfield(params_moco, 'medfilt'); params_moco.medfilt = 0; end                               % median filter motion output
if ~isfield(params_moco, 'plot_stuff'); params_moco.plot_stuff = 0; end

%params_moco.image_target = [];

params_moco.im_target_fname = [params.im_target_fname];

% descriptions of different moco methods
if params.moco_rigid_method == 1 % regular multiplane
    params_moco.num_iterations = 2;
    
    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 1;...
                              0.5 0.5 0.5];

    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 2 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [0 .2;... % 1
                              2 .2;...
                              2 .5;...
                              2 .5];
elseif params.moco_rigid_method == 21 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [1 .2;... % 1
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 22 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 2;...
                              0.5 0.5 2];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 23 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              4 1;...
                              4 1];
elseif params.moco_rigid_method == 24 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 2; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [5 .5;...
                              2 .2;...
                              2 .5;...
                              4 1;...
                              4 1];
                          
elseif params.moco_rigid_method == 25 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 6; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 2;...
                              0.5 0.5 1;...
                              0.5 0.5 0;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [.1 .01;...
                              .1 .01;...
                              .1 .01;...
                              .1 .01];
                                                
elseif params.moco_rigid_method == 3 % multiplane super noisy; dream/chrmine
    
    params_moco.num_iterations = 4; % 

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 32 % even more super noisy; dream/chrmine
    
    params_moco.num_iterations = 3; % 

    params_moco.smooth_std = [1 1 12;...
                              1 1 7;... 
                              1 1 5;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
end

% list of nonrigid methods                                     
if params.moco_nonrigid_method == 1
    params_moco.nonrigid_smooth_std = [0.5 0.5 6];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 60;
    params_moco.nonrigid_block_overlap = 40;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3];
    
elseif params.moco_nonrigid_method == 2
    params_moco.nonrigid_smooth_std = [0.5 0.5 1];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 50;
    params_moco.nonrigid_block_overlap = 30;
    params_moco.nonrigid_block_smooth = [0.5 0.5 1];

elseif params.moco_nonrigid_method == 3
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 40;
    params_moco.nonrigid_block_overlap = 30;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]
elseif params.moco_nonrigid_method == 4
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 30;
    params_moco.nonrigid_block_overlap = 15;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]
end

%%
save_dir_movie = [save_dir '\movies'];
save_dir_cuts = [save_dir '\preprocessing'];

params_moco.save_fname = params.save_fname;
params_moco.save_dir = save_dir_movie;

params.params_bidi = params_bidi;
params.params_moco = params_moco;

disp(save_fname);
fprintf('Pulse corp method = %d\n', params.align_pulse_crop_method);

if params.do_nonrigid
    fprintf('Moco nonrigind method = %d\n', params.moco_nonrigid_method);
else
    fprintf('Moco rigid method = %d\n', params.moco_rigid_method);
end


%%
if ~exist(save_dir_movie, 'dir'); mkdir(save_dir_movie); end
if ~exist(save_dir_cuts, 'dir'); mkdir(save_dir_cuts); end
if ~exist([save_dir_movie '\ave_proj'], 'dir'); mkdir([save_dir_movie '\ave_proj']); end

%%
colors1 = parula(params.num_planes);
proc_steps = '_cut';

if ~isempty(params_moco.im_target_fname)
    cuts_target_fname = [save_dir_cuts '\' params_moco.im_target_fname '_h5cutsdata.mat'];
    if exist(cuts_target_fname, 'file')
        moco_init_load = load(cuts_target_fname);
    end
end

%% load cuts data
cuts_fname = [save_dir_cuts '\' save_fname '_h5cutsdata.mat'];
if exist(cuts_fname, 'file')
    load_data = load(cuts_fname);
    cuts_data = load_data.cuts_data;
else
    cuts_data = cell(num_planes,1);
end

%% load
[~, ~, ext1] = fileparts(params.load_fname);

load_path = [params.load_dir '\' params.load_fname];

if ~numel(ext1)
    if exist(load_path, 'dir') % is a directory
        load_type = 1; 
    else
        error('provide correct file name, with extension, or directory')
    end
else
    if sum(strcmpi(ext1, {'.h5', '.hdf5'}))
        load_type = 3; 
    elseif sum(strcmpi(ext1, {'.tif', '.tiff'}))
        load_type = 2; 
    else
        error('Only accepts tiff, tif, h5, hdf5 or directory with tifs')
    end
end

% load types
% 1 = Prairie tiffs
% 2 = tiff stack (needs file name)
% 3 = h5 stack (needs file name)


Y = cell(num_planes,1);

if load_type == 1
    use_mpl_tags = 0;
    if num_planes > 1
        if params.use_prairie_mpl_tags
            use_mpl_tags = 1;
        end
    end
    if use_mpl_tags
        for n_pl = 1:num_planes
            Y{n_pl} = f_collect_prairie_tiffs4(load_path, params.mpl_tags{n_pl});
        end
    else
        Y_full = f_collect_prairie_tiffs4(load_path, params.prairie_chan_tag);
    end
elseif load_type == 2
    Y_full = bigread3(load_path, 1);
elseif load_type == 3
    Y_full = h5read(load_path, params.h5_movie_tag);
end

if num_planes > 1
    last_time = size(Y_full,3);
    params.ave_trace_full = squeeze(mean(mean(Y_full, 1),2));
    figure; plot(params.ave_trace_full);
    title('Full ave trace');
    for n_pl = 1:num_planes
        ind_mpl = n_pl:num_planes:last_time;
        Y{n_pl} = Y_full(:,:,ind_mpl);
    end
else
    Y{1} = Y_full;
end
clear Y_full;


%% compute cuts
if ~isfield(cuts_data{1}, 'vid_cuts_trace')
    cuts_data = cell(num_planes,1);
    [~, ~, T] = size(Y{1});
    vid_cuts_trace_all = true(T,1);
    for n_pl = 1:num_planes
        cuts_data{n_pl} = params;
        if num_planes>1
            cuts_data{n_pl}.title_tag = sprintf('_mpl%d_pl%d', num_planes, n_pl);
        else
            cuts_data{n_pl}.title_tag = '';
        end
        cuts_data{n_pl}.ave_trace = squeeze(mean(mean(Y{n_pl}, 1),2));
        cuts_data{n_pl} = f_compute_align_cuts(cuts_data{n_pl});
        vid_cuts_trace_all = vid_cuts_trace_all.*cuts_data{n_pl}.vid_cuts_trace;
    end
    for n_pl = 1:num_planes
        cuts_data{n_pl}.vid_cuts_trace = vid_cuts_trace_all;
    end
    save(cuts_fname, 'params', 'cuts_data');
end

%% apply cuts

for n_pl = 1:num_planes
    Y{n_pl}(:,:,~cuts_data{n_pl}.vid_cuts_trace) = [];
    if save_all_steps
        f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov');
    end
end

%% bidi fix

if do_bidi
    Y_pre_corr = Y;
    % Y = Y_pre_corr;
    if ~isfield(cuts_data{1}, 'bidi_out')
        % compute
        for n_pl = 1:num_planes
            fprintf('%s %s\n', save_fname, cuts_data{n_pl}.title_tag);
            params_bidi.title_tag = cuts_data{n_pl}.title_tag;
            [~, cuts_data{n_pl}.bidi_out] = f_fix_bidi_shifts3(Y{n_pl}, params_bidi);
        end
        save(cuts_fname, 'params', 'cuts_data');
    end
    
    bidi_shifts_all = cell(num_planes,1);
    for n_pl = 1:num_planes
        bidi_shifts_all{n_pl} = sum(cuts_data{n_pl}.bidi_out.best_shifts,2);
    end
    bidi_shifts_all = cat(2,bidi_shifts_all{:});
    % only use top 3
    if isfield(params_bidi, 'use_planes')
        mean_tag = num2str(params_bidi.use_planes(1):min([params_bidi.use_planes(2) num_planes]));
        mean_bidi_shifts = round(mean(bidi_shifts_all(:,params_bidi.use_planes(1):min([params_bidi.use_planes(2) num_planes])),2));
    else
        mean_tag = 'all';
        mean_bidi_shifts = round(mean(bidi_shifts_all(:,1:min([3 num_planes])),2));
    end
    
    figure; hold on;
    for n_pl = 1:num_planes
        plot(bidi_shifts_all(:,n_pl), 'color', colors1(n_pl, :))
    end
    plot(mean_bidi_shifts, 'k');
    title(['computed bidi shifts all planes; black-average pl ' mean_tag])
    

    % apply
    for n_pl = 1:num_planes
        tic;
        [d1, d2, T] = size(Y{n_pl});
        if cuts_data{n_pl}.bidi_out.params.do_interp
            
            num_blocks = ceil(T/params.block_size);
            deg_per_fov = 180 * cuts_data{n_pl}.bidi_out.params.laser_open_frac;
            y0 = 1:d1;
            deg0 = linspace(-deg_per_fov/2, deg_per_fov/2, d2);
            x0 = sin(deg0/360*2*pi);
            x0n = x0 - min(x0);
            x0n = x0n/max(x0n)*(d2-1)+1;
            x_coords = 1:d2;
            
            start1 = 1;
            
            fprintf('Applying nonrigid reg in blocks; block #/%d:', num_blocks)
            for n_bl = 1:num_blocks
                fprintf('..%d', n_bl);
                end1 = min((start1 + params.block_size-1), T);
                block_size2 = (end1 - start1+1);
                z0 = 1:block_size2;
                Y_temp = single(Y{n_pl}(:,:,start1:end1));
                [X_corr, Y0, Z0] = meshgrid(x_coords, y0, z0);
                [X_real, ~, ~] = meshgrid(x0n, y0, z0);
                Y_temp = interp3(X_corr, Y0, Z0, Y_temp, X_real, Y0, Z0);
                Y_temp = f_bidi_apply_shift(Y_temp, bidi_shifts_all(start1:end1,:));
                Y_temp = interp3(X_real, Y0, Z0, Y_temp, X_corr, Y0, Z0);
                Y{n_pl}(:,:,start1:end1) = uint16(Y_temp);
                start1 = end1 + 1;
            end
        else
            Y_temp = f_bidi_apply_shift(Y_temp, bidi_shifts_all);
        end
        fprintf('\nDone with nonrigid apply; compute time = %.1f\n', toc);
        
    end
    proc_steps = [proc_steps '_bidi'];
    
    if save_all_steps
        for n_pl = 1:num_planes
            f_save_mov_YS(Y{n_pl}(:,:,1:min(25000, size(Y{n_pl},3))), [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov');
        end
    end
end

%% moco

if do_moco
    Y_pre_corr = Y;
    %Y = Y_pre_corr;
    
    % correct movie to itself
    for n_pl = 1:num_planes
        if ~isfield(cuts_data{n_pl}, 'dsall') || params.overwrite_moco_rigid
            fprintf('%s %s\n', save_fname, cuts_data{n_pl}.title_tag);
            [~, mc_out] = f_mc_rigid(Y{n_pl}, params_moco);
            cuts_data{n_pl}.dsall = mc_out.dsall;
        end
    end
    
    save(cuts_fname, 'params', 'cuts_data');
    
    % add all coorection and use median actoss planes
    [~, dsall1_all, dsall1_all_mf] = f_mc_dsall_proc(cuts_data);
    
    if params_moco.medfilt
        dsall1_use = dsall1_all_mf;
    else
        dsall1_use = dsall1_all;
    end
    f_mc_plot_cuts_data(cuts_data, save_fname);
    
    %Y2 = Y_pre_moco;
    for n_pl = 1:num_planes
        Y{n_pl} = uint16(f_suite2p_reg_apply(Y{n_pl}, dsall1_use));
        %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
    end
    
    % apply global offset to while movie
    for n_pl = 1:num_planes
        Y_temp = single(Y{n_pl});
        cuts_data{n_pl}.image_target = mean(Y_temp,3);
        cuts_data{n_pl}.image_target_std = std(Y_temp,0,3);
    end

    ds_base_all = zeros(num_planes, 2);
    for n_pl = 1:num_planes
        cuts_data{n_pl}.ds_base = [0 0];
        if ~isempty(params_moco.moco_target_fname)
            if ~isempty(moco_init_load.cuts_data{n_pl}.image_target)
                cuts_data{n_pl}.image_target_external = moco_init_load.cuts_data{n_pl}.image_target;
                cuts_data{n_pl}.ds_base = f_suite2p_reg_compute(cuts_data{n_pl}.image_target, cuts_data{n_pl}.image_target_external, params_moco.reg_lambda_base);
            end
        end
        ds_base_all(n_pl, :) = cuts_data{n_pl}.ds_base;
    end
            
    save(cuts_fname, 'params', 'cuts_data');

    if ~isempty(params_moco.im_target_fname)
        figure; plot(ds_base_all); title('correction to external input database');
    end
    
    dsall1_use2 = ones(size(dsall1_use));
    for n_pl = 1:num_planes
        Y{n_pl} = uint16(f_suite2p_reg_apply(Y{n_pl}, dsall1_use2.*ds_base_all(n_pl, :)));
        %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
    end
     
    if params.moco_zero_edge
        for n_pl = 1:num_planes
            Y{n_pl} = f_mc_zero_edges(Y{n_pl}, dsall1_use, ds_base_all(n_pl, :));
        end
    end
    
    proc_steps = [proc_steps '_moco'];
    if save_all_steps
        for n_pl = 1:num_planes
            f_save_mov_YS(Y{n_pl}(:,:,1:min(25000, size(Y{n_pl},3))), [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov')
        end
    end
    
    if params.do_nonrigid
        Y_pre_corr = Y;
        %Y = Y_pre_corr;
        
        for n_pl = 1:num_planes
            if ~isfield(cuts_data{n_pl}, 'nr_corr_data') || params.overwrite_moco_nonrigid
                [~, mc_out2] = f_mc_nonrigid(Y{n_pl}, params_moco);
                cuts_data{n_pl}.nr_corr_data = mc_out2.nr_corr;
            end
        end

        f_mc_plot_nrcuts_data(cuts_data, save_fname)
        
        for n_pl = 1:num_planes
            Y{n_pl} = f_mc_apply_nonrigid_corr(Y{n_pl}, cuts_data{n_pl}.nr_corr_data);
        end
        proc_steps = [proc_steps '_nonrigid'];
        if 1%save_all_steps
            for n_pl = 1:num_planes
                f_save_mov_YS(Y{n_pl}(:,:,1:min(25000, size(Y{n_pl},3))), [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov')
            end
        end
    end
    
end

%% save
params.params_bidi = params_bidi;
params.params_moco = params_moco;

for n_pl = 1:num_planes
    params.save_mov_path = [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag '.h5'];
    params.cuts_data = cuts_data{n_pl};
    
    f_save_mov_YS(Y{n_pl}, params.save_mov_path, '/mov');
    
    if params.trim_output_num_frames
        params.save_mov_path_trim = [save_dir_movie '\' save_fname cuts_data{n_pl}.title_tag proc_steps '_trim.hdf5'];
        f_save_mov_YS(Y{n_pl}(:,:,1:round(params.trim_output_num_frames)), params.save_mov_path_trim, '/mov');  
    end
    
    if save_indiv_h5info
        save([save_dir '\' save_fname cuts_data{n_pl}.title_tag '_h5cutsinfo.mat'], 'params');
    end
    
    tmp_fig = figure; imagesc(squeeze(mean(Y{n_pl},3)));
    title([save_fname ' Ave prjection ' cuts_data{n_pl}.title_tag], 'Interpreter', 'none');
    axis tight equal;
    saveas(tmp_fig,[save_dir_movie '\ave_proj\' save_fname cuts_data{n_pl}.title_tag '_ave_proj.fig']);
end

disp('Done')

end