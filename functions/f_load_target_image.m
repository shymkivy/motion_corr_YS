function params = f_load_target_image(params)
params_moco = params.params_moco;

% maybe broken let me know
moco_image_target = cell(params.num_planes,1);

fpath1 = [params.save_dir_cuts '\' params_moco.im_target_fname];
if ~isempty(params_moco.im_target_fname)
    if iscell(params_moco.im_target_fname) % this should be list of string names
        if numel(params_moco.im_target_fname) == num_planes
            for n_pl = 1:num_planes
                moco_image_target{n_pl} = [];
                if numel(params_moco.im_target_fname{n_pl})
                    loaded1 = 0;
                    [~, ~, ext1] = fileparts(params_moco.im_target_fname{n_pl});
                    if sum(strcmpi(ext1, {'.tif', '.tiff'}))
                        if exist(fpath1, 'file')
                            moco_image_target{n_pl} = imread(fpath1, ext1(2:end));
                            loaded1 = 1;
                        end
                    end
                    if ~loaded1
                        warning('Moco image for pl%d unreadable, not using it: %s', n_pl, params_moco.im_target_fname);
                    end
                end
            end
        else
            error('im_target_fname input should be cell containing number of planes used, with corresponding average image inputs, or empty strings');
        end
    else % either string or mat file
        [~, ~, ext1] = fileparts(params_moco.im_target_fname);
        
        if sum(strcmpi(ext1, {'.tif', '.tiff'}))
            if num_planes == 1
                moco_image_target = imread(fpath1, ext1(2:end));
            else
                error('im_target_fname string tif inputs only available for one plane data, use cell instead');
            end
        else
            cuts_target_fname = [save_dir_cuts '\' params_moco.im_target_fname '_h5cutsdata.mat'];
            if exist(cuts_target_fname, 'file')
                moco_init_load = load(cuts_target_fname);
                for n_pl = 1:num_planes
                    if ~isempty(moco_init_load.cuts_data{n_pl}.image_target)
                        moco_image_target{n_pl} = moco_init_load.cuts_data{n_pl}.image_target;
                    else
                        warning("specified im_target_fname cuts mat file has no target image");
                    end
                end
            else
                warning("specified im_target_fname cuts mat file doesn't exist, using no target");
            end
        end
    end
end
params.params_moco.image_target = moco_image_target;

end