function Y = f_load_mov(fpath, params)

if ~exist('params', 'var'); params = struct(); end
if ~isfield(params, 'h5_movie_tag'); params.h5_movie_tag = '/mov'; end
if ~isfield(params, 'manually_split_planes'); params.manually_split_planes = 0; end   % 0, or number of planes to split

disp(params.save_fname);
[~, ~, ext1] = fileparts(fpath);

if ~numel(ext1)
    if exist(fpath, 'dir') % is a directory
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

if load_type == 1   % prairie directory
    Y = f_load_prairie(fpath);
elseif load_type == 2
    Y{1,1} = bigread4(fpath);
elseif load_type == 3
    Y{1,1} = h5read(fpath, params.h5_movie_tag);
end

if params.manually_split_planes > 1 % code for manually splitting into planes
    num_planes = params.manually_split_planes;
    siz1 = size(Y{1,1});
    Y2 = Y{1,1};
    Y = cell(num_planes,1);
    for n_pl = 1:num_planes
        ind_mpl = n_pl:num_planes:siz1(3);
        Y{n_pl} = Y2(:,:,ind_mpl);
    end
    clear Y2;
end

end