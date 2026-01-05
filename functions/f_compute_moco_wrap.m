function [Y, cuts_data, params] = f_compute_moco_wrap(Y, cuts_data, params)

%Y_pre_corr = Y;
%Y = Y_pre_corr;

fprintf('Moco rigid method = %d\n', params.moco_rigid_method);
if params.do_nonrigid
    fprintf('Moco nonrigind method = %d\n', params.moco_nonrigid_method); 
end

params_moco = params.params_moco;
params = f_load_target_image(params);

% correct movie to itself
for n_pl = 1:params.num_planes
    if ~isfield(cuts_data{n_pl}, 'dsall') || params.overwrite_moco_rigid
        fprintf('%s %s\n', params.save_fname, cuts_data{n_pl}.title_tag);
        [~, mc_out] = f_mc_rigid(Y{n_pl}, params.params_moco);
        cuts_data{n_pl}.dsall = mc_out.dsall;
        cuts_data{n_pl}.corr_all = mc_out.corr_all;
        cuts_data{n_pl}.corr_all_z = mc_out.corr_all_z;
    end
end

save(params.cuts_fname, 'params', 'cuts_data');

% add all coorection and use median actoss planes
[~, dsall1_use] = f_mc_dsall_proc(cuts_data, params_moco.medfilt);

f_mc_plot_cuts_data(cuts_data, params.save_fname);

%Y2 = Y_pre_moco;
for n_pl = 1:params.num_planes
    mov_type = class(Y{n_pl});
    Y{n_pl} = f_suite2p_reg_apply(Y{n_pl}, dsall1_use);
    Y{n_pl} = f_set_dtype(Y{n_pl}, mov_type);
    %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
end

% apply global offset to input target image
for n_pl = 1:params.num_planes
    Y_temp = single(Y{n_pl});
    cuts_data{n_pl}.image_target = mean(Y_temp,3);
    cuts_data{n_pl}.image_target_std = std(Y_temp,0,3);
end

ds_base_all = zeros(params.num_planes, 2);
for n_pl = 1:params.num_planes
    cuts_data{n_pl}.ds_base = [0 0];
    if ~isempty(params.params_moco.image_target{n_pl})
        cuts_data{n_pl}.image_target_external = params.params_moco.image_target{n_pl};
        cuts_data{n_pl}.ds_base = f_suite2p_reg_compute(cuts_data{n_pl}.image_target, cuts_data{n_pl}.image_target_external, params_moco.reg_lambda_base);
    end
    ds_base_all(n_pl, :) = cuts_data{n_pl}.ds_base;
end
        
save(params.cuts_fname, 'params', 'cuts_data');

if ~isempty(params_moco.im_target_fname)
    figure; plot(ds_base_all); title('correction to external input database');
end

dsall1_use2 = ones(size(dsall1_use));
for n_pl = 1:params.num_planes
    mov_type = class(Y{n_pl});
    if sum(ds_base_all(n_pl, :))
        Y{n_pl} = f_suite2p_reg_apply(Y{n_pl}, dsall1_use2.*ds_base_all(n_pl, :));
        Y{n_pl} = f_set_dtype(Y{n_pl}, mov_type);
        %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
    end
end
 
if params.moco_zero_edge
    for n_pl = 1:params.num_planes
        Y{n_pl} = f_mc_zero_edges(Y{n_pl}, dsall1_use, ds_base_all(n_pl, :));
    end
end

proc_steps = [params.proc_steps '_moco'];
if params.save_all_steps
    f_save_mov_wrap(Y, params, [cuts_data{n_pl}.title_tag params.proc_steps]);
end

if params.do_nonrigid
    %Y_pre_corr = Y;
    %Y = Y_pre_corr;
    
    for n_pl = 1:params.num_planes
        if ~isfield(cuts_data{n_pl}, 'nr_corr_data') || params.overwrite_moco_nonrigid
            [~, mc_out2] = f_mc_nonrigid(Y{n_pl}, params.params_moco);
            cuts_data{n_pl}.nr_corr_data = mc_out2.nr_corr;
        end
    end

    f_mc_plot_nrcuts_data(cuts_data, params.save_fname)
    
    for n_pl = 1:params.num_planes
        Y{n_pl} = f_mc_apply_nonrigid_corr(Y{n_pl}, cuts_data{n_pl}.nr_corr_data);
    end
    params.proc_steps = [proc_steps '_nonrigid'];
    if params.save_all_steps
        f_save_mov_wrap(Y, params, [cuts_data{n_pl}.title_tag params.proc_steps]);
    end
end

end