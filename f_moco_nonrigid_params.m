function params_moco = f_moco_nonrigid_params(params_moco, moco_nonrigid_method)

% list of nonrigid methods                                     
if moco_nonrigid_method == 1
    params_moco.nonrigid_smooth_std = [0.5 0.5 6];
    params_moco.nonrigid_reg_lambda = [.5];
    params_moco.nonrigid_block_size = 60;
    params_moco.nonrigid_block_overlap = 40;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3];
    
elseif moco_nonrigid_method == 2
    params_moco.nonrigid_smooth_std = [0.5 0.5 1];
    params_moco.nonrigid_reg_lambda = [.5];
    params_moco.nonrigid_block_size = 60;
    params_moco.nonrigid_block_overlap = 10;
    params_moco.nonrigid_block_smooth = [0.5 0.5 1];

elseif moco_nonrigid_method == 3
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [.5];
    params_moco.nonrigid_block_size = 40;
    params_moco.nonrigid_block_overlap = 30;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]

elseif moco_nonrigid_method == 4
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [.5];
    params_moco.nonrigid_block_size = 30;
    params_moco.nonrigid_block_overlap = 15;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]
end

end