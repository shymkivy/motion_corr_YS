function [Y, params] = f_load_mov_old(params)

if ~isfield(params, 'num_planes'); params.num_planes = 1; end
if ~isfield(params, 'use_prairie_mpl_tags'); params.use_prairie_mpl_tags = 1; end
if ~isfield(params, 'prairie_mpl_tags'); params.prairie_mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; end % multiplane data tags in prairie
if ~isfield(params, 'prairie_chan_tag'); params.prairie_chan_tag = 'Ch2'; end
if ~isfield(params, 'h5_movie_tag'); params.h5_movie_tag = '/mov'; end
if ~isfield(params, 'load_tif_format'); params.load_tif_format = ''; end
if ~isfield(params, 'prairie_multipage_tif'); params.prairie_multipage_tif = true; end

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

    % multipage format they started stacking multiple frames in one file
    % multipage tseries, they stack time frames in mulpile files
    % multipage z-stack, they stack planes in same file and each volume is different file
    % old tseries each frame is separate file
    % old zseries each plane has name tag Ch2_000001 ...
    % old zseries with no tags is just sequentially saved planes and volumes by time
    if (~params.prairie_multipage_tif) && (num_planes > 1) && (params.use_prairie_mpl_tags)
        tags1 = params.prairie_mpl_tags;
    else
        tags1 = params.prairie_chan_tag;
    end
    Y_full = f_collect_prairie_tiffstacks2(load_path, tags1);
    
    if ~params.prairie_multipage_tif
        if (num_planes > 1) && (~params.use_prairie_mpl_tags)
            % old zseries way of where planes are registerred sequentally
            % (pix x pix) x planes x frames
            Y2 = cat(3,Y_full{:});
            siz1 = size(Y2);
            last_time = siz1(3);
            for n_pl = 1:num_planes
                ind_mpl = n_pl:num_planes:last_time;
                Y{n_pl,1} = Y2(:,:,ind_mpl);
            end
            clear Y2;
        else
            % (pix x pix) x (planes x frames)
            for n_pl = 1:num_planes
                Y{n_pl,1} = cat(3,Y_full{n_pl,:});
            end
        end
    else
        if num_planes > 1
            % (pix x pix x planes) x frames
            Y2 = cat(4,Y_full{:}); 
            siz1 = size(Y2);
            for n_pl = 1:siz1(3)
                Y{n_pl,1} = reshape(Y2(:,:,n_pl,:), siz1(1), siz1(2), siz1(4));
            end
            clear Y2;
        else
            % (pix x pix x frames) x files
            Y{1,1} = cat(3,Y_full{:});
        end
        
    end
elseif load_type == 2
    Y_full = bigread4(load_path);
elseif load_type == 3
    Y_full = h5read(load_path, params.h5_movie_tag);
end

if load_type > 1
    if num_planes > 1
        siz1 = size(Y_full);
        last_time = prod(siz1(3:end));
        Y2 = reshape(Y_full, siz1(1), siz1(2), last_time);
        figure; plot(params.ave_trace_full);
        title('Full ave trace');
        params.ave_trace_full = squeeze(mean(mean(Y2, 1),2));
        for n_pl = 1:num_planes
            ind_mpl = n_pl:num_planes:last_time;
            Y{n_pl} = Y2(:,:,ind_mpl);
        end
        clear Y2;
    end
end

clear Y_full;

end