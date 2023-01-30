function f_save_mov_YS(data, file, h5tag)

[filepath,~,ext] = fileparts(file);

if isempty(filepath)
    filepath = pwd;
end
    
if ~exist(filepath, 'dir')
    mkdir(filepath)
end


if strcmpi(ext,'.bin')
    ops.RegFile = fullfile(file);
    fid = fopen(ops.RegFile, 'w');
    dwrite = int16(data);
    fwrite(fid, dwrite, 'int16');
    fclose(fid);
end

if strcmpi(ext,'.h5') || strcmpi(ext,'.hdf5')
    if ~exist('h5tag', 'var')
        h5tag = '/Y';
    end

    if exist(file, 'file')
        delete(file)
        disp(['Replacing previous file ' file]);
    end
    h5create(file, h5tag, size(data), 'Datatype', 'uint16');
    h5write(file,h5tag,uint16(data));
end



end