function [Y, cuts_data, params] = f_compute_bidi_wrap(Y, cuts_data, params)
%Y_pre_corr = Y;
% Y = Y_pre_corr;

params_bidi = params.params_bidi;
if ~isfield(cuts_data{1}, 'bidi_out')
    % compute
    for n_pl = 1:params.num_planes
        fprintf('Bidi cuts for %s %s; plane %d of %d\n', params.save_fname, cuts_data{n_pl}.title_tag, n_pl, params.num_planes);
        params_bidi.title_tag = cuts_data{n_pl}.title_tag;
        [~, cuts_data{n_pl}.bidi_out] = f_fix_bidi_shifts3(Y{n_pl}, params_bidi);
    end
    save(params.cuts_fname, 'params', 'cuts_data');
end

bidi_shifts_all = cell(params.num_planes,1);
for n_pl = 1:params.num_planes
    bidi_shifts_all{n_pl} = sum(cuts_data{n_pl}.bidi_out.best_shifts,2);
end
bidi_shifts_all = cat(2,bidi_shifts_all{:});
% only use top 3
if isfield(params_bidi, 'use_planes')
    mean_tag = num2str(params_bidi.use_planes(1):min([params_bidi.use_planes(2) params.num_planes]));
    mean_bidi_shifts = round(mean(bidi_shifts_all(:,params_bidi.use_planes(1):min([params_bidi.use_planes(2) params.num_planes])),2));
else
    mean_tag = 'all';
    mean_bidi_shifts = round(mean(bidi_shifts_all(:,1:min([3 params.num_planes])),2));
end

colors1 = parula(params.num_planes);
figure; hold on;
for n_pl = 1:params.num_planes
    plot(bidi_shifts_all(:,n_pl), 'color', colors1(n_pl, :))
end
plot(mean_bidi_shifts, 'k');
title(['computed bidi shifts all planes; black-average pl ' mean_tag])

% apply
for n_pl = 1:params.num_planes
    tic;
    T = size(Y{n_pl},3);
    type1 = class(Y{n_pl});
    if cuts_data{n_pl}.bidi_out.params.do_interp
        num_blocks = ceil(T/params.block_size);
        start1 = 1;
        fprintf('Applying bidi in blocks; block #/%d:', num_blocks)
        for n_bl = 1:num_blocks
            fprintf('..%d', n_bl);
            end1 = min((start1 + params.block_size-1), T);
            Y_temp = single(Y{n_pl}(:,:,start1:end1));
            if cuts_data{n_pl}.bidi_out.do_interp
                Y_temp = f_bidi_res_galvo_interp(Y_temp, cuts_data{n_pl}.bidi_out.params.laser_open_frac, 1, true);
            end
            Y_temp = f_bidi_apply_shift(Y_temp, bidi_shifts_all(start1:end1,:));
            if cuts_data{n_pl}.bidi_out.do_interp
                Y_temp = f_bidi_res_galvo_interp(Y_temp, cuts_data{n_pl}.bidi_out.params.laser_open_frac, 1, false);
            end
            Y{n_pl}(:,:,start1:end1) = f_set_dtype(Y_temp, type1);
            start1 = end1 + 1;
        end
    else
        Y{n_pl} = f_bidi_apply_shift(Y{n_pl}, bidi_shifts_all);
    end
    fprintf('\nDone with bidi apply; compute time = %.1f\n', toc);
end

params.proc_steps = [params.proc_steps '_bidi'];
if params.save_all_steps
    f_save_mov_wrap(Y, params, [cuts_data{n_pl}.title_tag params.proc_steps]);
end

end