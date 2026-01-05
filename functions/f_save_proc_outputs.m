function f_save_proc_outputs(Y, cuts_data,  params)

for n_pl = 1:params.num_planes
    params.save_mov_path = [params.save_dir_movie '\' params.save_fname cuts_data{n_pl}.title_tag '.h5'];
    params.cuts_data = cuts_data{n_pl};
    
    f_save_mov_YS(Y{n_pl}, params.save_mov_path, '/mov');
    
    if params.save_indiv_h5info
        save([params.save_dir '\' params.save_fname cuts_data{n_pl}.title_tag '_h5cutsinfo.mat'], 'params');
    end
    
    if params.save_in_parts
        num_frames = size(Y{1},3);
        part_size = params.save_in_parts_size;
        intervals = 1:part_size:num_frames;
        for n_int = 1:numel(intervals)
            f_save_mov_YS(Y{n_pl}(:,:,intervals(n_int):min(intervals(n_int)+part_size-1, num_frames)), [params.save_dir_movie '\' params.save_fname cuts_data{n_pl}.title_tag '_pt' num2str(n_int) '.h5'], '/mov');
        end
    end

    tmp_fig = figure; imagesc(squeeze(mean(Y{n_pl},3)));
    title([params.save_fname ' Ave prjection ' cuts_data{n_pl}.title_tag], 'Interpreter', 'none');
    axis tight equal;
    saveas(tmp_fig,[params.save_dir_movie '\ave_proj\' params.save_fname cuts_data{n_pl}.title_tag '_ave_proj.fig']);
end

end

