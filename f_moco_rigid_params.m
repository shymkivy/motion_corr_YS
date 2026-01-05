function params_moco = f_moco_rigid_params(params_moco, moco_rigid_method)

% descriptions of different moco methods
if moco_rigid_method == 1 % regular multiplane
    params_moco.num_iterations = 2;
    
    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 1;...
                              0.5 0.5 0.5];

    params_moco.reg_lambda = [1, 2];

elseif moco_rigid_method == 12 % regular multiplane
    params_moco.num_iterations = 4;
    
    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 1;...
                              0.5 0.5 0.5];

    params_moco.reg_lambda = [1, 2];

elseif moco_rigid_method == 2 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [0, 2];
elseif moco_rigid_method == 21 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [1, 2];
                          
elseif moco_rigid_method == 22 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 2;...
                              0.5 0.5 2];
                          
    params_moco.reg_lambda = [1, 2];
                          
elseif moco_rigid_method == 23 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [1, 2, 2, 4];
elseif moco_rigid_method == 24 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [1, 2];
                          
elseif moco_rigid_method == 25 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 6; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 2;...
                              0.5 0.5 1;...
                              0.5 0.5 0;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [.1];

elseif moco_rigid_method == 26 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 3; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 2;...
                              0.5 0.5 1;... % was 3 for missmatch
                              0.5 0.5 0;...
                              0.5 0.5 0;...
                              0.5 0.5 0;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [0.01, .1];

elseif moco_rigid_method == 27 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 6; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [1, 2, 2, 4];

elseif moco_rigid_method == 3 % multiplane super noisy; dream/chrmine
    
    params_moco.num_iterations = 4; % 

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1, 2];

elseif moco_rigid_method == 31 % 3 extended (echo 10hz datasets)
    
    params_moco.num_iterations = 6; % 

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 2;
                              0.5 0.5 1;
                              0.5 0.5 1];
                          
    params_moco.reg_lambda = [1, 2];
                          
elseif moco_rigid_method == 32 % even more super noisy; dream/chrmine
    
    params_moco.num_iterations = 3; % 

    params_moco.smooth_std = [1 1 12;...
                              1 1 7;... 
                              1 1 5;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1, 2];
elseif moco_rigid_method == 0 % no iterations
    params_moco.num_iterations = 1;
    params_moco.smooth_std = [0.5 0.5 0.5];
    params_moco.reg_lambda = [2 .5];
end


end