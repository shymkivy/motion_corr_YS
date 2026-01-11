function params = f_set_params(params)

%% error check
if ~isfield(params, 'load_dir'); error('Input load dir'); end
if ~numel(params.load_dir); error('Input load dir'); end
if ~isfield(params, 'save_dir'); error('Input save dir'); end
if ~numel(params.save_dir); error('Input save dir'); end

%% default params
% loading
%if ~isfield(params, 'num_planes'); params.num_planes = 1; end
%if ~isfield(params, 'use_prairie_mpl_tags'); params.use_prairie_mpl_tags = 1; end
%if ~isfield(params, 'prairie_mpl_tags'); params.prairie_mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; end % multiplane data tags in prairie
%if ~isfield(params, 'prairie_multipage_tif'); params.prairie_multipage_tif = true;  end
%if ~isfield(params, 'prairie_chan_tag'); params.prairie_chan_tag = 'Ch2'; end
if ~isfield(params, 'h5_movie_tag'); params.h5_movie_tag = '/mov'; end

% saving
if ~isfield(params, 'save_all_steps'); params.save_all_steps = 0; end                       % during save
if ~isfield(params, 'save_all_steps_frames'); params.save_all_steps_frames = 10000; end                         % number of frames to save in all steps movies
if ~isfield(params, 'save_in_parts'); params.save_in_parts = 0; end                         % also save same movie in smaller parts to be able to load with imagej
if ~isfield(params, 'save_in_parts_size'); params.save_in_parts_size = 30000; end           % num frames per part
if ~isfield(params, 'save_indiv_h5info'); params.save_indiv_h5info = 1; end                 % during save
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

%% bidirectional shift fix params
if ~isfield(params, 'params_bidi'); params.params_bidi = struct(); end
params_bidi = params.params_bidi;

if ~isfield(params_bidi, 'smooth_std'); params_bidi.smooth_std = [0.5 0.5 1]; end       % smoothing video before bidi [m, n, T]
if ~isfield(params_bidi, 'reg_lambda'); params_bidi.reg_lambda = [2 0.5]; end           % regularizations for step 1 and 2
if ~isfield(params_bidi, 'fix_range'); params_bidi.fix_range = -15:15; end              % maximal movement
if ~isfield(params_bidi, 'num_iterations'); params_bidi.num_iterations = 1; end
if ~isfield(params_bidi, 'use_planes'); params_bidi.use_planes = [1 3]; end             % for multiplane, whick planes to use
if ~isfield(params_bidi, 'plot_stuff'); params_bidi.plot_stuff = 0; end

params.params_bidi = params_bidi;

%% moco params
if ~isfield(params, 'params_moco'); params.params_moco = struct(); end 
params_moco = params.params_moco;

if ~isfield(params_moco, 'high_val_cut_thresh'); params_moco.high_val_cut_thresh = 0.01; end    % to reduce super bright signals
if ~isfield(params_moco, 'reg_lambda_base'); params_moco.reg_lambda_base = 1; end          % some regularization
if ~isfield(params_moco, 'medfilt'); params_moco.medfilt = 0; end                               % median filter motion output
if ~isfield(params_moco, 'plot_stuff'); params_moco.plot_stuff = 0; end

%params_moco.image_target = [];

%params_moco.suite2P_matlab_path = params.suite2P_matlab_path;
params_moco.im_target_fname = [params.im_target_fname];

params.params_moco = params_moco;

params.params_moco = f_moco_rigid_params(params.params_moco, params.moco_rigid_method);
params.params_moco = f_moco_nonrigid_params(params.params_moco, params.moco_nonrigid_method);
%%
params.proc_steps = '';

params.save_dir_movie = [params.save_dir '\movies'];
params.save_dir_cuts = [params.save_dir '\preprocessing'];

params.params_moco.save_fname = params.save_fname;
params.params_moco.save_dir = params.save_dir_movie;

if ~exist(params.save_dir_movie, 'dir'); mkdir(params.save_dir_movie); end
if ~exist(params.save_dir_cuts, 'dir'); mkdir(params.save_dir_cuts); end
if ~exist([params.save_dir_movie '\ave_proj'], 'dir'); mkdir([params.save_dir_movie '\ave_proj']); end

if ~isfield(params, 'load_fname')
    if isfield(params, 'dset_name')
        params.load_fname = [params.dset_name(1:end-1) '-00' params.dset_name(end)];
    end
end

end