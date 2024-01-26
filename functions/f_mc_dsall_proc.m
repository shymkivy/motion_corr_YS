function [dsall1, dsall1_all] = f_mc_dsall_proc(cuts_data, do_medfilt)

num_planes = numel(cuts_data);

dsall1 = cell(num_planes, 1);
for n_pl = 1:num_planes
    dsall1{n_pl} = sum(cat(3,cuts_data{n_pl}.dsall{:}),3);
end
dsall1_all = median(cat(3,dsall1{:}),3);

if do_medfilt
    dsall1_all = medfilt1(dsall1_all, 3);
end


end