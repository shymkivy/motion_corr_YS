function Y_full = f_collect_prairie_tiffstacks(load_path, tags)
% load the PrairieView output tif images and combine them into single
% movie stack h5 or tif. h5 is much faster

if ~exist(load_path, 'dir'); Error(['Pipeline Error: Directory ' load_path ' does not exist']); end
if ~exist('tags', 'var'); tags = '_'; end % else; tag2 = ['*' tag1]; 

if iscell(tags)
    num_tag_im = zeros(numel(tags),1);
    for n_tag = 1:numel(tags)
        dir_names = dir([load_path, '\*' tags{n_tag} '*.tif']);
        num_tag_im(n_tag) = numel(dir_names);
    end
    tags2 = tags(logical(num_tag_im));
else
    tags2 = {tags};
end

dir_names = dir([load_path, '\*' tags2{1} '*.tif']);
if ~issorted({dir_names.name})
    warning('file names order is not sorted');
end
info = imfinfo([load_path, '\',  dir_names(1).name], 'tif');
num_tags = numel(tags2);
num_file = numel(dir_names);
num_frames = numel(info);
if ~num_frames; Error('Pipeline Error: No tiff files in specified folters'); end

Y_full = cell(num_tags, num_file);
for n_tag = 1:numel(tags2)
    dir_names = dir([load_path, '\*' tags2{n_tag} '*.tif']);
    if ~issorted({dir_names.name})
        warning('file names order is not sorted, doing sorting, check results');
        [~, idx1] = sort({dir_names.name});
        dir_names = dir_names(idx1);
    end
    fpath = [load_path, '\',  dir_names(1).name];
    tremp_frame = imread(fpath, 1);

    dim = size(tremp_frame);
    if ~isempty(dir_names)
        num_file = numel(dir_names);
        info = imfinfo([load_path, '\',  dir_names(1).name], 'tif');
        num_frames = numel(info);
        info = imfinfo([load_path, '\',  dir_names(num_file).name], 'tif');
        num_frames_last = numel(info);
        frames_all = num_frames*(num_file-1) + num_frames_last;
        hh = waitbar(0, sprintf('Loading Prairie tiffs; %s', strrep(tags2{n_tag},'_',' ')));
        for n_fl = 1:(num_file-1)
            num_frames_pre = num_frames*(n_fl-1);
            Y = zeros([dim, num_frames], 'uint16');
            for n_fr = 1:num_frames
                Y(:,:,n_fr) = imread([load_path, '\', dir_names(n_fl).name], n_fr);
                waitbar((num_frames_pre+n_fr)/frames_all, hh);
            end
            Y_full{n_tag,n_fl} = Y;
        end
        n_fl = num_file;
        num_frames_pre = num_frames*(n_fl-1);
        Y = zeros([dim, num_frames_last], 'uint16');
        for n_fr = 1:num_frames_last
            Y(:,:,n_fr) = imread([load_path, '\', dir_names(n_fl).name], n_fr);
            waitbar((num_frames_pre+n_fr)/frames_all, hh);
        end
        Y_full{n_tag,n_fl} = Y;
        close(hh)
    end
end

end