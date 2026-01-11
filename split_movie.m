
pwd2 = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath([pwd2 '\functions']));

mov_dir = 'D:\VR\data_proc\L\movies\';
%mov_dir = 'G:\data\Auditory\caiman_data_missmatch\movies';

fname = {'L_10_21_25_im1.h5',...
         %'M10_im12_A1_ammn2_5_31_20.h5',...
         %'M10_im13_A1_freq_grating1_5_31_20.h5',...
         %'M10_im14_A1_freq_grating2_5_31_20.h5',...
         };

max_size = 20000;  % 30 seems close to max that imagej can open
num_save = 1;

for n_fl = 1:numel(fname)
    [~, fname2, ext1] = fileparts(fname{n_fl});

    Y = h5read([mov_dir '\' fname2 ext1],'/mov');

    [d1, d2, T] = size(Y);
    
    num_blocks = ceil(T/max_size);
    start1 = 1;
    for n_bl = 1:min([num_save, num_blocks])
        end1 = min((start1 + max_size - 1), T);
        f_save_mov_YS(Y(:,:,start1:end1), sprintf('%s\\%s_pt%d%s', mov_dir, fname2, n_bl, ext1), '/mov');
        start1 = end1 + 1;
    end
end

