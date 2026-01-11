%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This runs the process_mov pipeline in a loop through files listed in
%   xlsx file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clear;
close all; 

pwd2 = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath([pwd2 '\functions']));

% %%
% load_dir = {'H:\data\AC\2021\',...
%             'H:\data\AC\2022\'};
%          
load_dir = {'D:\VR\'};%...
            %'I:\mouse\auditory\2018'};

%save_dir = {'F:\AC_data\caiman_data_echo\'};%,...
save_dir = {'D:\VR\data_proc\'};%,...
%save_dir = {'F:\AC_data\caiman_data_dream\'};

params.dset_table_fpath = 'F:\VR\data_proc\VR_data.xlsx';

params.limit.dset_name =        '';
params.limit.experiment =       '';
params.limit.mouse_id =         'RL';
params.limit.mouse_tag =        '';
params.limit.FOV_num =          '';

params.save_all_steps = 0;

%%
AC_data = f_s0_parse_tab_data(params);

mouse_id_all = unique(AC_data.mouse_id, 'stable');

%% set default params
AC_data.do_moco(isnan(AC_data.do_moco)) = 1;
AC_data.do_bidi(isnan(AC_data.do_bidi)) = 0;
AC_data.mc_zero_edge(isnan(AC_data.mc_zero_edge)) = 1;
AC_data.mc_rigid_met(isnan(AC_data.mc_rigid_met)) = 1;
AC_data.mc_do_nonrigid(isnan(AC_data.mc_do_nonrigid)) = 0;
%%

fprintf('Running %d dsets total...\n', size(AC_data,1))

for n_ms = 1:numel(mouse_id_all)
    AC_data2 = AC_data(strcmpi(AC_data.mouse_id, mouse_id_all{n_ms}),:);
    fprintf('Mouse id %s; %d dsets...\n', mouse_id_all{n_ms}, size(AC_data2,1));
    % check if folder exists
    
    % set moco target as first in list
    idx2 = logical(sum(AC_data2.im_num == unique(AC_data2.mc_to_dset)',2));
    AC_data2 = [AC_data2(idx2,:); AC_data2(~idx2,:)];
    
    for n_dset = 1:size(AC_data2,1)
        do_s0 = true;

        load_path = sprintf('%s\\%s\\',AC_data2.dset_name{n_dset}, AC_data2.mouse_id{n_dset});      % adjust here for relative loading structure

        if isstring(load_dir)
            load_dir = {load_dir};
        end
        
        fold_exist = false(numel(load_dir),1);
        for n_data_dir = 1:numel(load_dir)
            load_dir2 = load_dir{n_data_dir};
            if exist([load_dir2 '\' load_path], 'dir')
                fold_exist(n_data_dir) = 1;
            end
        end
        
        if ~sum(fold_exist)
            do_s0 = 0;
            warning(['Data directory does not exist: ' load_path])
        else
            params.load_dir = [load_dir{fold_exist} '\' load_path];
        end
        
        if do_s0
            cdset = AC_data2(n_dset,:);

            params.save_fname = sprintf('%s_%s_im%d', cdset.mouse_id{1}, cdset.dset_name{1}, cdset.im_num);
            
            if ~iscell(save_dir)
                save_dir = {save_dir}
            end

            % check it output already exists
            num_match = 0;
            if iscell(save_dir)
                for n_dir = 1:numel(save_dir)
                    dir_list = dir([save_dir{n_dir} '\' cdset.mouse_id{1} '\movies\*.h5']);
                    dir_names = {dir_list.name};
                    for n_file = 1:numel(dir_names)
                        pat1 = strfind(dir_names{n_file}, params.save_fname);
                        if ~isempty(pat1)
                            num_match = num_match + 1;
                        end
                    end
                end
            end

            if ~num_match
                params.num_planes = cdset.mpl;
                params.do_moco = cdset.do_moco;
                params.do_nonrigid = cdset.mc_do_nonrigid;
                params.moco_zero_edge = cdset.mc_zero_edge;
                params.do_bidi = cdset.do_bidi;
                params.moco_rigid_method = cdset.mc_rigid_met;
                params.moco_nonrigid_method = cdset.mc_nonrigid_met;
                params.dset_name = cdset.dset_name{1};
                params.save_dir = [save_dir{1} '\' cdset.mouse_id{1}];
                if or(cdset.align_pulse_crop_method == 0, cdset.align_pulse_crop_method == 2)
                    params.align_pulse_crop_method = cdset.align_pulse_crop_method;
                else
                    params.align_pulse_crop_method = 1; % default is 1 = auto; 2 = manual; 0 = full movie
                end
                
                params.im_target_fname = '';
                if params.do_moco
                    if ~isempty(cdset.mc_to_dset)
                        if ~isnan(cdset.mc_to_dset)
                            if cdset.im_num ~= cdset.mc_to_dset
                                source_dset = cdset.mc_to_dset;
                                source_dset_idx = AC_data2.im_num == source_dset;
                                fname_dset1 = sprintf('%s_im%d_%s_%s', AC_data2.mouse_id{source_dset_idx}, AC_data2.im_num(source_dset_idx), AC_data2.dset_name{source_dset_idx}, AC_data2.mouse_tag{source_dset_idx});
                                params.im_target_fname = fname_dset1;
                                fprintf('Using target registration im from %s\n', params.im_target_fname);
                            end
                        end
                    end
                end
                
                params.load_fname = cdset.recording_name{1};
                params = f_set_params(params);
                
                % load
                Y = f_load_mov([params.load_dir '\\' params.load_fname], params);

                % compute cuts (cut out alignment pulses)
                [Y, cuts_data, params] = f_compute_alignment_pulse_cuts(Y, params);

                % bidi fix
                if params.do_bidi
                    [Y, cuts_data, params] = f_compute_bidi_wrap(Y, cuts_data, params);
                end

                % moco
                if params.do_moco
                    [Y, cuts_data, params] = f_compute_moco_wrap(Y, cuts_data, params);
                end

                f_save_proc_outputs(Y, cuts_data,  params);

            else
                fprintf('%s already exists, moving on...\n', params.save_fname)
            end

        end
    end
end

fprintf('All done\n')