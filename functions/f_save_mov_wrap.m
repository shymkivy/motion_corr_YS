function f_save_mov_wrap(Y, params, title_tag)

for n_pl = 1:params.num_planes
    f_save_mov_YS(Y{n_pl}(:,:,1:min(params.save_all_steps_frames, size(Y{n_pl},3))), [params.save_dir_movie '\' params.save_fname title_tag '.h5'], '/mov'); % 
end

end