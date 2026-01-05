function Y = f_bidi_res_galvo_interp(Y, laser_open_frac, interp_density, undo_interp)
% for large bidi fixes habe to first undo the interpolation, then redo after
% to undo it set undo to true
% redo it set undo to fales
if ~exist('undo_interp', 'var'); undo_interp = true; end

[d1, d2, T] = size(Y);

deg_per_fov = 180 * laser_open_frac;
y0 = 1:d1;
%z0 = 1:T;

deg0 = linspace(-deg_per_fov/2, deg_per_fov/2, d2/interp_density);
x0 = sin(deg0/360*2*pi);
x0n = x0 - min(x0);
x0n = x0n/max(x0n)*(d2-1)+1;

x_coords = 1:interp_density:d2;

if undo_interp
    %Y = interp3(x_coords, y0', z0, Y, x0n, y0', z0, 'linear');
    for n_f = 1:T
        Y(:,:,n_f) = interp2(x_coords, y0', Y(:,:,n_f), x0n, y0', 'linear');
    end
else
    %Y = interp3(x0n, y0', z0, Y, x_coords, y0', z0, 'linear');
    for n_f = 1:T
        Y(:,:,n_f) = interp2(x0n, y0', Y(:,:,n_f), x_coords, y0', 'linear');
    end
end

end