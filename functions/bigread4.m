function imData=bigread4(path_to_file,sframe,num2read, form)
%reads tiff files in Matlab bigger than 4GB, allows reading from sframe to sframe+num2read-1 frames of the tiff - in other words, you can read page 200-300 without rading in from page 1.
%based on a partial solution posted on Matlab Central (http://www.mathworks.com/matlabcentral/answers/108021-matlab-only-opens-first-frame-of-multi-page-tiff-stack)
%Darcy Peterka 2014, v1.0
%Darcy Peterka 2014, v1.1 
%Darcy Peterka 2016, v1.2(bugs to dp2403@columbia.edu)
%Eftychios Pnevmatikakis 2016, v1.3 (added hdf5 support)
%Yuriy 2018, v3 changed so you could read big stacks (>4gb) again
%Yuriy 2025, v4 made compatible with multiple strips per frame, alse sped up by removing the cell in the end
%Program checks for bit depth, whether int or float, and byte order.  Assumes uncompressed, non-negative (i.e. unsigned) data.
%
% Usage:  my_data=bigread('path_to_data_file, start frame, num to read);
% "my_data" will be your [M,N,frames] array.
%Will do mild error checking on the inputs - last two inputs are optional -
%if nargin == 2, then assumes second is number of frames to read, and
%starts at frame 1

[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif')
    
    %get image info
    info = imfinfo(path_to_file);
    
    blah=size(info);
    
    % image stacks over ~4GB only have info for first string, so we need to use description to get num of frames and caculate offsets 
    big_stack = 0;
    if (prod(blah) >1) && (blah(1) == 1)
        he=info.ImageDescription;
        numFramesStr = regexp(he, 'images=(\d*)', 'tokens');
        numFrames = str2double(numFramesStr{1}{1});
        
        big_stack = numFrames > 1;
    else
        numFrames = blah(1);
    end
    
    num_tot_frames=numFrames;

    %should add more error checking for args... very ugly code below.  works
    %for me after midnight though...
    if nargin<2
        sframe = 1;
    end
    if ~exist('num2read', 'var') || isempty(num2read)
        num2read=numFrames-sframe+1;
    else
        numFrames = min(numFrames, numFrames-sframe+1);
    end
    if sframe<=0
        sframe=1;
    end
    if num2read<1
        num2read=1;
    end
    if sframe>num_tot_frames
        sframe=num_tot_frames;
        num2read=1;
        disp('starting frame has to be less than number of total frames...');
    end
    if (num2read+sframe<= num_tot_frames+1)
        lastframe=num2read+sframe-1;
    else
        num2read=numFrames-sframe+1;
        lastframe=num_tot_frames;
        disp('Hmmm...just reading from starting frame until the end');
    end

    bd=info(1).BitDepth;
    he=info(1).ByteOrder;
    bo=strcmp(he,'big-endian');
    if ~exist('form', 'var') || isempty(form)
        if (bd==64)
            form='uint64';
        elseif(bd==32)
            form='uint32';
        elseif (bd==16)
            form='uint16';
        elseif (bd==8)
            form='uint8';
        end
    end


    % Use low-level File I/O to read the file
    fp = fopen(path_to_file , 'rb');
    % The StripOffsets field provides the offset to the first strip. Based on
    % the INFO for this file, each image consists of 1 strip.

    he=info(1).StripOffsets;
    num_strips = numel(he);
    RowsPerStrip = info(1).RowsPerStrip;
    
    %finds the offset of each strip in the movie.  Image does not have to have
    %uniform strips, but needs uniform bytes per strip/row.
    ofds=zeros(numFrames,num_strips);
    if ~big_stack
        for ii=1:numFrames
            ofds(ii,:)=info(ii).StripOffsets;
        end
    else
        % if we have big stack, calculate offsets manually.
        ofds(1) = info(1).StripOffsets;
        strip_shift = info(1).StripByteCounts;
        for ii=2:numFrames
            ofds(ii)= ofds(ii-1) + strip_shift;
        end
    end
    
    sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(num2read+sframe-1),' of ',num2str(num_tot_frames), ' total frames'];
    disp(sframemsg)
    pause(.2)
    %go to start of first strip
    fseek(fp, ofds(1,1), 'bof');
    %framenum=numFrames;
    framenum=num2read;

    he_w=info(1).Width;
    he_h=info(1).Height;
    
    % specify format
    if (bo)
        machinefmt = 'ieee-be';
    else
        machinefmt = 'ieee-le';
    end
    
    if strcmpi(form,'uint64')
        machinefmt = [machinefmt, '.l64'];
    end
    
    imData = zeros(he_h, he_w, framenum, form);
    strip_starts = (0:(num_strips-1))*RowsPerStrip+1;
    for cnt = sframe:lastframe
        for n_strip = 1:num_strips
            %cnt;
            fseek(fp,ofds(cnt, n_strip),'bof');
            tmp1 = fread(fp, [he_w RowsPerStrip], form, 0, machinefmt)';
            imData(strip_starts(n_strip):strip_starts(n_strip)+RowsPerStrip-1,:,cnt-sframe+1) = cast(tmp1,form);
        end
        %imData{cnt-sframe+1}=frame1;
    end
    
    %imData=cell2mat(imData);
    %imData=reshape(imData,[he_h,he_w,framenum]);
    fclose(fp);
    disp('Finished reading images');
    
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5')
    info = hdf5info(path_to_file);
    dims = info.GroupHierarchy.Datasets.Dims;
    if nargin < 2
        sframe = 1;
    end
    if nargin < 3
        num2read = dims(end)-sframe+1;
    end
    num2read = min(num2read,dims(end)-sframe+1);
    imData = h5read(path_to_file,'/mov',[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]);
else
    error('Unknown file extension. Only .tiff and .hdf5 files are currently supported');
end