function nr_corr = f_mc_compute_nonrigid_motion(Y, params)

[d1, d2, T] = size(Y);

block_size = params.nonrigid_block_size;
block_overlap = params.nonrigid_block_overlap;
target_frame = params.target_frame;

ops.regLambda = params.reg_lambda;
ops.maxregshift = block_size; 

num_split_m = round((d1 - 2*block_overlap)/(block_size - 2*block_overlap));
num_split_n = round((d1 - 2*block_overlap)/(block_size - 2*block_overlap));

cent_coords_m = round(linspace(1+block_size/2, d1 - block_size/2, num_split_m));
cent_coords_n = round(linspace(1+block_size/2, d2 - block_size/2, num_split_n));

num_blocks = num_split_m*num_split_n;

coords_mn = zeros(num_blocks, 2);
for nn = 1:num_split_n
    for mm = 1:num_split_m
        n_bl = (nn-1)*num_split_n + mm;
        coords_mn(n_bl, :) = [cent_coords_m(mm) cent_coords_n(nn)];
    end
end

figure; hold on;
im1 = imagesc(target_frame); 
axis equal tight;
for n_bl = 1:num_blocks
    plot(coords_mn(n_bl, 2), coords_mn(n_bl, 1), 'ko')
    rectangle('Position',[coords_mn(n_bl, 2)-block_size/2 coords_mn(n_bl, 1)-block_size/2 block_size block_size])
end
im1.Parent.YDir = 'reverse';

dsall_m = zeros(num_blocks, T);
dsall_n = zeros(num_blocks, T);
corr_all_mn = zeros(num_blocks, T);
tic;
fprintf('Block #/%d', num_blocks);
for n_bl = 1:num_blocks
    fprintf('..%d', n_bl);

    m_low = max(round((coords_mn(n_bl, 1) - block_size/2)),1);
    m_high = min(round((coords_mn(n_bl, 1) + block_size/2)),d1);

    n_low = max(round((coords_mn(n_bl, 2)  - block_size/2)),1);
    n_high = min(round((coords_mn(n_bl, 2) + block_size/2)),d2);

    idx_m = m_low:m_high;
    idx_n = n_low:n_high;

    ref_block = target_frame(idx_m, idx_n);
    Y_cut = Y(idx_m, idx_n, :);
    
    if 0
        figure; imagesc(ref_block);
        figure; imagesc(Y_cut(:,:,1));
    end
    [dsall_mn, ~, ops1] = f_suite2p_reg_compute_cpu(Y_cut, ref_block, ops);
    dsall_m(n_bl,:) = dsall_mn(:,1);
    dsall_n(n_bl,:) = dsall_mn(:,1);
    corr_all_mn(n_bl,:) = ops1{1}.CorrFrame;
end
fprintf('\nDone; compute duration=%.1fsec \n', toc);


nr_corr.dsall_m = dsall_m;
nr_corr.dsall_n = dsall_n;
nr_corr.corr_all_mn = corr_all_mn;
nr_corr.cent_coords_n = cent_coords_n2;
nr_corr.cent_coords_m = cent_coords_m2;
nr_corr.params = params;

if 0
    figure; hold on;
    for n_bl = 1:num_blocks
        plot(dsall_m(n_bl,:) - (n_bl-1)*5)
    end
    title('dsall m')

    figure; hold on;
    for n_bl = 1:num_blocks
        plot(dsall_n(n_bl,:) - (n_bl-1)*5)
    end
    title('dsall n')
    
    thresh1 = max(corr_all_mn(:));
    figure; hold on; axis tight
    for n_bl = 1:num_blocks
        plot(corr_all_mn(n_bl,:) - (n_bl-1)*thresh1*1.1);
        plot(zeros(T,1) - (n_bl-1)*thresh1*1.1, '--r');
        plot(zeros(T,1) + thresh1 - (n_bl-1)*thresh1*1.1, '--g');
    end
    title('corr')
end

end