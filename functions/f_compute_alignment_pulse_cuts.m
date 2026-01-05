function [Y, cuts_data, params] = f_compute_alignment_pulse_cuts(Y, params)

fprintf('Pulse corp method = %d\n', params.align_pulse_crop_method);

% load cuts data
params.num_planes = size(Y,1);
params.cuts_fname = [params.save_dir_cuts '\' params.save_fname '_h5cutsdata.mat'];
if exist(params.cuts_fname, 'file')
    load_data = load(params.cuts_fname);
    cuts_data = load_data.cuts_data;
else
    cuts_data = cell(params.num_planes,1);
end

if ~isfield(cuts_data{1}, 'vid_cuts_trace')
    T_all = zeros(params.num_planes,1);
    for n_pl = 1:params.num_planes
        T_all(n_pl) = size(Y{n_pl},3);
    end
    vid_cuts_trace_all = true(max(T_all),1);
    for n_pl = 1:params.num_planes
        cuts_data{n_pl} = params;
        if params.num_planes>1
            cuts_data{n_pl}.title_tag = sprintf('_mpl%d_pl%d', params.num_planes, n_pl);
        else
            cuts_data{n_pl}.title_tag = '';
        end
        T1 = size(Y{n_pl},3);
        cuts_data{n_pl}.ave_trace = reshape(mean(mean(Y{n_pl}, 1),2), T1, 1);
        cuts_data{n_pl} = f_compute_align_cuts(cuts_data{n_pl});
        vid_cuts_trace_all(1:T_all(n_pl)) = vid_cuts_trace_all(1:T_all(n_pl)).*cuts_data{n_pl}.vid_cuts_trace;
    end
    for n_pl = 1:params.num_planes
        cuts_data{n_pl}.vid_cuts_trace = vid_cuts_trace_all(1:T_all(n_pl));
    end
    save(params.cuts_fname, 'params', 'cuts_data');
end

% apply cuts
params.proc_steps = [params.proc_steps, '_cut'];
for n_pl = 1:params.num_planes
    Y{n_pl}(:,:,~cuts_data{n_pl}.vid_cuts_trace) = [];
end

if params.save_all_steps
    f_save_mov_wrap(Y, params, [cuts_data{n_pl}.title_tag params.proc_steps]);
end

end