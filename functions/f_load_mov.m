function [Y, params] = f_load_mov(params)

if ~isfield(params, 'num_planes'); params.num_planes = 1; end
if ~isfield(params, 'use_prairie_mpl_tags'); params.use_prairie_mpl_tags = 1; end
if ~isfield(params, 'mpl_tags'); params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; end % multiplane data tags in prairie
if ~isfield(params, 'prairie_chan_tag'); params.prairie_chan_tag = 'Ch2'; end
if ~isfield(params, 'h5_movie_tag'); params.h5_movie_tag = '/mov'; end
if ~isfield(params, 'load_tif_format'); params.load_tif_format = ''; end

num_planes = params.num_planes;

[~, ~, ext1] = fileparts(params.load_fname);

load_path = sprintf('%s\\%s', params.load_dir, params.load_fname);
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
    Y_full = bigread3(load_path, 1, [], params.load_tif_format);
elseif load_type == 3
    Y_full = h5read(load_path, params.h5_movie_tag);
end

if num_planes > 1
    if ~use_mpl_tags
        last_time = size(Y_full,3);
        params.ave_trace_full = squeeze(mean(mean(Y_full, 1),2));
        figure; plot(params.ave_trace_full);
        title('Full ave trace');
        for n_pl = 1:num_planes
            ind_mpl = n_pl:num_planes:last_time;
            Y{n_pl} = Y_full(:,:,ind_mpl);
        end
    end
else
    Y{1} = Y_full;
end
clear Y_full;

end