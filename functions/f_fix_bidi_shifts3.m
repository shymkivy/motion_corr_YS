function [Y_bidi, bidi_out] = f_fix_bidi_shifts3(Y, params)
% Y should be whole movie (d1 * d2 * T)
if ~exist('params', 'var'); params = struct(); end
if ~isfield(params, 'fix_range'); params.fix_range = -20:20; end
if ~isfield(params, 'smooth_std'); params.smooth_std = [0 0.5 2]; end   % smoothing in x, y, z
if ~isfield(params, 'num_iterations'); params.num_iterations = 1; end
if ~isfield(params, 'plot_stuff'); params.plot_stuff = 0; end
if ~isfield(params, 'title_tag'); params.title_tag = 0; end
if ~isfield(params, 'do_interp'); params.do_interp = 0; end
if ~isfield(params, 'normalize'); params.normalize = 1; end

if ~isfield(params, 'laser_open_frac')
    % this is the fraction of full sine period of res galvo that microscope uses
    % prairie version 5
    params.laser_open_frac = .65;%.797; % measured 51us up 13us down (says 63 total)28.8; % 14.4
end

smooth_std = params.smooth_std;
fix_range = params.fix_range;
num_iterations = params.num_iterations;
plot_stuff = params.plot_stuff;
title_tag = params.title_tag;
normalize = params.normalize;
interp_density = 1;

[d1, d2, T] = size(Y);
type1 = class(Y);

Y2 = single(Y);
fprintf('Bidi fix; ');
%% undo interpolations that correct for res galvo movement
% for small shifts is is not needed 
if params.do_interp
    fprintf('Interpolating\n');
    Y2 = f_bidi_res_galvo_interp(Y2, params.laser_open_frac, interp_density, true);
else
    fprintf('Not interpolating\n');
end

%% correct bidi
best_shifts = zeros(T,num_iterations);
best_corr_vals = zeros(T, num_iterations);
best_zero_vals = zeros(T, num_iterations);
num_range = numel(fix_range);

idx_1 = 1:2:d1;
idx_2 = 2:2:d1;

for n_rep = 1:num_iterations
    fprintf('Bidi fix iter %d; ', n_rep)
    
    Y_odd = Y2(idx_1,:,:);
    Y_even = Y2(idx_2,:,:);
    
    %f_save_mov_YS(Y_odd(:,:,1:min(25000, size(Y_odd,3))), 'odd.h5', '/mov');
    %f_save_mov_YS(Y_even(:,:,1:min(25000, size(Y_even,3))), 'even.h5', '/mov');
    
    if sum(smooth_std>0)
        fprintf('smoothing: ');
        tic;
        Y_odd = f_smooth_movie(Y_odd, smooth_std);
        Y_even = f_smooth_movie(Y_even, smooth_std);
        fprintf('smooth duration=%.1fsec; ', toc);
    end
    
    tic;
    fprintf('bidi estimation: ');
    for n_fr = 1:T
        frame_odd = double(Y_odd(:,:,n_fr));
        frame_even = double(Y_even(:,:,n_fr));
        if normalize
            frame_odd = frame_odd/norm(frame_odd, 'fro');
            frame_even = frame_even/norm(frame_even, 'fro');
        end

        corr_vals = zeros(num_range,1);
        max_shift = d2 - max(abs(fix_range));
        for n_sh = 1:num_range
            sh1 = fix_range(n_sh);
            
            start_even = max([1-sh1, 1]);
            end_even = start_even + max_shift-1;
            frame_even_sh = frame_even(:, start_even:end_even);
            start_odd = max([1+sh1, 1]);
            end_odd = start_odd + max_shift-1;
            frame_odd2 = frame_odd(:,start_odd:end_odd);

            if normalize
                frame_even_sh = frame_even_sh/norm(frame_even_sh, 'fro');
                frame_odd2 = frame_odd2/norm(frame_odd2, 'fro');
            end

            corr_val1 = frame_odd2.*frame_even_sh;
            corr_vals(n_sh) = mean(mean(corr_val1));
            if sh1 == 0
                best_zero_vals(n_fr, n_rep) = mean(mean(corr_val1));
            end
        end
        [best_corr_vals(n_fr, n_rep), idx] = max(corr_vals);
        best_shifts(n_fr, n_rep) = fix_range(idx);
        
    end
    fprintf('compute duration=%.1fsec; ', toc);
    
    % apply to the raw
    tic;
    best_shifts_all = sum(best_shifts,2);
    Y2 = f_bidi_apply_shift(Y2, best_shifts_all);
    fprintf('apply durration=%.1fsec\n', toc);
    %f_save_mov_YS(Y_bidi, ['C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\movies\test_bidi_iter' num2str(n_rep) '.h5'], '/mov');
    %f_save_mov_YS(Y_bidi, ['F:\AC_data\test_bidi_iter' num2str(n_rep) '.h5'], '/mov');
    
%     k = 14;
%     Y_odd_frame_sh = circshift(Y_odd_frame_ip, round(k), 2);
%     Y_even_frame_sh = circshift(Y_even_frame_ip, -round(k), 2);
%     
%     rgb_im2 = zeros(128,256,3);
%     rgb_im2(:,:,1) = Y_odd_frame_sh;
%     rgb_im2(:,:,2) = Y_even_frame_sh; 
%     rgb_im2 = rgb_im2/max(rgb_im2(:))*1.5;
%     
%     figure; imagesc(rgb_im2)
    
end

% undo  interp
if params.do_interp
    fprintf('Deinterpolating\n');
    Y2 = f_bidi_res_galvo_interp(Y2, params.laser_open_frac, interp_density, false);
end
Y_bidi = f_set_dtype(Y2, type1);
bidi_out.params = params;
bidi_out.best_shifts = best_shifts;
bidi_out.best_corr_vals = best_corr_vals;
bidi_out.best_zero_vals = best_zero_vals;

if plot_stuff
    
    fr_idx = 2000:3000;
    figure; 
    imagesc(mean(Y2(:,:,fr_idx),3));
    title(sprintf('example ave; laser frac = %.2f', params.laser_open_frac));

    rgb_im = zeros(round(d1/2),d2,3);
    rgb_im(:,:,1) = mean(Y2(idx_1,:,fr_idx),3);
    rgb_im(:,:,2) = mean(Y2(idx_2,:,fr_idx),3);
    rgb_im = rgb_im/max(rgb_im(:))*1.5;
    figure; 
    imagesc(rgb_im)
    title(sprintf('example odd even lines; laser frac = %.2f', params.laser_open_frac));
    
    colors1 = parula(num_iterations);
    figure;
    ax1 = subplot(2,1,1); hold on; axis tight;
    for n_rep = 1:num_iterations
        plot(best_shifts(:,n_rep), 'color', colors1(n_rep,:));
    end
    plot(sum(best_shifts,2), 'k');
    title(['best shift ' title_tag], 'interpreter', 'none');
    ax2 = subplot(2,1,2); hold on; axis tight;
    for n_rep = 1:num_iterations
        plot(best_corr_vals(:,n_rep)-best_zero_vals(:,n_rep), 'color', colors1(n_rep,:));
    end
    %plot(sum(best_corr_vals-best_zero_vals,2), 'k', 'Linewidth', 1);
    title(['best corr -zero ' title_tag], 'interpreter', 'none');
    linkaxes([ax1 ax2], 'x');
end
    

end