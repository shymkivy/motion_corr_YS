function Y = f_load_prairie(load_path)
% new prairie load where it reads info and loads appropriately

if ~exist(load_path, 'dir'); Error(['Pipeline Error: Directory ' load_path ' does not exist']); end

dir_names = dir([load_path, '\*.tif']);
if ~issorted({dir_names.name})
    warning('file names order is not sorted');
end
% reads info only once and assumes all files in dir are saved in same fmt
info = imfinfo([load_path, '\',  dir_names(1).name], 'tif');
info0 = if_get_info(info);

% get data either from companion ome file or the first image in dir
companion = dir([load_path, '\*.ome']);
if numel(companion)
    fp = fopen([load_path, '\' , companion.name], 'r' );
    text1 = fread(fp,[1 inf],'*char');
    fclose(fp);
else
    text1  = info(1).ImageDescription;
end
idx_sta = regexp(text1,'<TiffData');
idx_end = regexp(text1, '</TiffData>');

num_frames_all = numel(idx_sta);
firstC = zeros(num_frames_all,1);
firstT = zeros(num_frames_all,1);
firstZ = zeros(num_frames_all,1);
fnames = cell(num_frames_all,1);

for n_fl = 1:num_frames_all
    temp_txt = text1(idx_sta(n_fl):idx_end(n_fl));
    [~,idxe1] = regexp(temp_txt,'FirstC="', 'once');
    [idxs2,~] = regexp(temp_txt(idxe1+1:end),'"', 'once');
    firstC(n_fl) = str2double(temp_txt(idxe1+1:idxe1+idxs2-1));
    
    [~,idxe1] = regexp(temp_txt,'FirstT="', 'once');
    [idxs2,~] = regexp(temp_txt(idxe1+1:end),'"', 'once');
    firstT(n_fl) = str2double(temp_txt(idxe1+1:idxe1+idxs2-1));
    
    [~,idxe1] = regexp(temp_txt,'FirstZ="', 'once');
    [idxs2,~] = regexp(temp_txt(idxe1+1:end),'"', 'once');
    firstZ(n_fl) = str2double(temp_txt(idxe1+1:idxe1+idxs2-1));

    [~,idxe1] = regexp(temp_txt,'FileName="', 'once');
    [idxs2,~] = regexp(temp_txt(idxe1+1:end),'"', 'once');
    fnames{n_fl} = temp_txt(idxe1+1:idxe1+idxs2-1);
end
[fnames_uq,fstart,~] = unique(fnames);

frame_counts = [diff(fstart); num_frames_all-fstart(end)+1];

num_file = numel(fnames_uq);

tic;
info1 = info0;
imData = zeros(info0.height, info0.width, num_frames_all, info0.form);
hh = waitbar(0, 'Loading Prairie tiffs');
for n_fl = 1:num_file
    info1.num_frames = frame_counts(n_fl);
    imData(:,:,fstart(n_fl):fstart(n_fl)+frame_counts(n_fl)-1) = if_read_frame([load_path, '\',  fnames_uq{n_fl}], info1);
    waitbar(n_fl/num_file, hh);
end
close(hh)
toc

planes_uq = unique(firstZ);
num_planes = numel(planes_uq);

Y = cell(num_planes,1);
if num_planes > 1
    for n_pl = 1:num_planes
        idx_pl = planes_uq(n_pl) == firstZ;
        Y{n_pl,1} = imData(:,:,idx_pl);
    end
else
    Y{1,1} = imData;
end

end

function info2 = if_get_info(info)
    
    
    bd=info(1).BitDepth;
    if (bd==64)
        form='uint64';
    elseif(bd==32)
        form='uint32';
    elseif (bd==16)
        form='uint16';
    elseif (bd==8)
        form='uint8';
    end
    
    he=info(1).ByteOrder;
    if strcmp(he,'big-endian')
        machinefmt = 'ieee-be';
    else
        machinefmt = 'ieee-le';
    end
    
    info2 = struct();
    info2.form = form;
    info2.machinefmt = machinefmt;
    info2.num_frames = numel(info);
    info2.height = info(1).Height;
    info2.width = info(1).Width;
    info2.num_strips = numel(info(1).StripOffsets);
    info2.StripOffsets = zeros(info2.num_frames, info2.num_strips);
    for n_fr = 1:info2.num_frames
        info2.StripOffsets(n_fr,:) = info(n_fr).StripOffsets;
    end
    info2.RowsPerStrip = info(1).RowsPerStrip;
    info2.strip_starts = (0:(info2.num_strips-1))*info2.RowsPerStrip+1;
end

function imData = if_read_frame(fpath, info2)

    imData = zeros(info2.height, info2.width, info2.num_frames, info2.form);
    fp = fopen(fpath , 'rb');
    for n_fr = 1:info2.num_frames
        for n_str = 1:info2.num_strips
            fseek(fp, info2.StripOffsets(n_fr, n_str),'bof');
            tmp1 = fread(fp, [info2.width info2.RowsPerStrip], info2.form, 0, info2.machinefmt)';
            imData(info2.strip_starts(n_str):info2.strip_starts(n_str)+info2.RowsPerStrip-1,:,n_fr) = cast(tmp1, info2.form);
        end
    end
    fclose(fp);

end
