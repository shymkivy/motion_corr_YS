% computes registration offset of data from mean image ops.mimg
% subpixel registration is kriging interpolation
% inputs:
%%% data = [ly, lx, nFrames]
%%% ops = options
%%% removeMean - if 1, then subtract mean offsets and realign to mean offset
%%% (for pre-registration alignment)
% outputs:
%%% offsets = [nFrames, 2]
%%% corr    = [nFrames]
%%% regdata = [ly, lx, nFrames] (only computed if removeMean is not empty!)

% slightly modified from suite2p

function [dv,corr, corr_z] = regoffKriging_YS2(data, ops, removeMean)

refImg = ops.mimg;
subpixel = getOr(ops, {'subPixel' 'SubPixel'}, 10); % subpixel factor
useGPU = getOr(ops, 'useGPU', false);
phaseCorrelation = getOr(ops, {'phaseCorrelation', 'PhaseCorrelation'}, true);
% maximum shift allowed
maxregshift = getOr(ops, 'maxregshift', round(.1*max(ops.Ly,ops.Lx))); 
% slope on taper mask preapplied to image. was 2, then 1.2
maskSlope   = getOr(ops, 'maskSlope', 2); 
% SD pixels of gaussian smoothing applied to correlation map (MOM likes .6)
smoothSigma = 1.15;
reg_lambda = getOr(ops, 'regLambda', 0); % some regularization param ys

% if subpixel is still inf, threshold it for new method
if ops.kriging
    subpixel = min(10, subpixel);
end

[ly, lx, nFrames] = size(data);

% Taper mask
[ys, xs] = ndgrid(1:ly, 1:lx);
ys = abs(ys - mean(ys(:)));
xs = abs(xs - mean(xs(:)));
mY      = max(ys(:)) - 4;
mX      = max(xs(:)) - 4;
maskMul = single(1./(1 + exp((ys - mY)/maskSlope)) ./(1 + exp((xs - mX)/maskSlope)));
maskOffset = mean(refImg(:))*(1 - maskMul);

% Smoothing filter in frequency domain
hgx = exp(-(((0:lx-1) - fix(lx/2))/smoothSigma).^2);
hgy = exp(-(((0:ly-1) - fix(ly/2))/smoothSigma).^2);
hg = hgy'*hgx;
fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));

% fft of reference image 
eps0          = single(1e-10);
cfRefImg = conj(fftn(refImg));
if phaseCorrelation
    absRef   = abs(cfRefImg);
    cfRefImg = cfRefImg./(eps0 + absRef) .* fhg;
end

if useGPU
    batchSize = getBatchSize(ly*lx);
    eps0      = gpuArray(eps0);
    cfRefImg    = gpuArray(cfRefImg);
    maskMul = gpuArray(maskMul);
    maskOffset = gpuArray(maskOffset);    
else
    batchSize = 1000;
end

% allow max shifts +/- lcorr
lpad   = 3;
lcorr  = min(maxregshift, floor(min(ly,lx)/2)-lpad);

% only need a small kernel +/- lpad for smoothing
[x1,x2] = ndgrid([-lpad:lpad]);
xt = [x1(:) x2(:)]';
if useGPU
    xt = gpuArray(single(xt));
end

if ops.kriging
    % compute kernels for regression (smoothing factor)
    sigL     = .85; % kernel width in pixels
    Kx = kernelD(xt,xt,sigL*[1;1]);
    linds = [-lpad:1/subpixel:lpad];
    [x1,x2] = ndgrid(linds);
    xg = [x1(:) x2(:)]';
    if useGPU
        xg = gpuArray(single(xg));
    end
    Kg = kernelD(xg,xt,sigL*[1;1]);
    Kmat = Kg/Kx;
end

% loop over batches
dv = zeros(nFrames, 2);
corr = zeros(nFrames, 1);
corr_mean = zeros(nFrames, 1);
corr_std = zeros(nFrames, 1);

% adding penalty for longer displacement YS
%meshgrid()
pen_mask_full = 1./(sqrt((xs*reg_lambda).^2 + (ys*reg_lambda).^2)+1);
%[y2, x2] = ndgrid(-lpad:1:lpad);
%pen_mask_ccmat = 1./(sqrt((x2*reg_lambda(2)).^2 + (y2*reg_lambda(2)).^2)+1);
%figure; imagesc(pen_mask);
%figure; plot(pen_mask(:,27));

nBatches = ceil(nFrames/batchSize);
for bi = 1:nBatches
    fi = ((bi - 1)*batchSize + 1):min(bi*batchSize, nFrames);
    
    if useGPU
        batchData = gpuArray(single(data(:,:,fi)));
    else
        batchData = single(data(:,:,fi));
    end
    
    corrMap = fft2(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, batchData)));
    
    %keyboard;
    if phaseCorrelation
        corrMap = bsxfun(@times, corrMap./(eps0 + abs(corrMap)), cfRefImg);
    else
        corrMap = bsxfun(@times, corrMap, cfRefImg);
    end
    
    % compute correlation matrix
    corrClip = real(ifft2(corrMap));
    corrClip = fftshift(fftshift(corrClip, 1), 2);

    %% subpixel registration
    % allow only +/- lcorr shifts
    idxx = floor(lx/2)+1+[-lcorr:lcorr];
    idxy = floor(ly/2)+1+[-lcorr:lcorr];
    cc0 = corrClip(idxy, idxx,:);
    %.*pen_mask;      % added penalty mask for long movements YS
    pen_mask2 = pen_mask_full(idxy, idxx);
    
    cc0_reg = cc0.*pen_mask2;

    if subpixel > 1
        % kriging subpixel

        [cmax,ii]   = max(reshape(cc0_reg, [], numel(fi)),[],1);
        
        [iy, ix] = ind2sub((2*lcorr+1) * [1 1], ii);
        
        dl       = single(-lpad:1:lpad);
        if useGPU
            dl   = gpuArray(dl);
            ccmat = gpuArray.zeros(numel(dl), numel(dl), numel(fi), 'single');
        else
            ccmat = zeros(numel(dl), numel(dl), numel(fi), 'single');
        end
        mxpt        = [iy(:)+floor(ly/2) ix(:)+floor(lx/2)] - lcorr;
        for j = 1:numel(fi)
            % matrix +/- lpad surrounding max point
            ccmat(:, :, j) = corrClip(mxpt(j,1)+dl, mxpt(j,2)+dl, j);
        end
        %
        %ccmat = ccmat.*pen_mask_ccmat; % YS added
        ccmat = reshape(ccmat,[], numel(fi));
        if ops.kriging
            % regress onto subsampled grid
            ccb         = Kmat * ccmat;
            
            % find max of grid
            [cx,ix]     = max(ccb, [], 1);
            [ix11,ix21] = ind2sub(numel(linds)*[1 1],ix);
            mdpt        = floor(numel(linds)/2)+1;
            dv0         = bsxfun(@minus, ([ix11' ix21'] - mdpt)/subpixel + mxpt, ...
                [floor(ly/2) floor(lx/2)]) - 1;
        else
            yshift      = xt(1,:) * ccmat;
            xshift      = xt(2,:) * ccmat;
            dv0         = bsxfun(@rdivide, [yshift' xshift'], sum(ccmat, 1)') + mxpt;
            dv0         = bsxfun(@minus, dv0, [floor(ly/2) floor(lx/2)] + 1);
            if isfinite(subpixel)
                dv0 = round(dv0 * subpixel) / subpixel;
            end
            cx     = max(ccmat, [], 1);
        end
        dv(fi,:) = gather_try(dv0);
        corr(fi)  = gather_try(cx);
    % otherwise just take peak of matrix
    else
        
        [cmax,iy]  = max(cc0_reg,[],1);
        [cx, ix]   = max(cmax,[],2);

        iy = reshape(iy(sub2ind([size(iy,2) size(iy,3)], ix(:), (1:size(iy,3))')),...
            1, 1, []);
        
        dv0 = [iy(:)-lcorr ix(:)-lcorr]-1;
        dv(fi,:)  = gather_try(dv0);
        corr(fi) = gather_try(cx(:));
    end
    
    corr_mean(fi) = mean(mean(cc0,2),1);
    corr_std(fi) = std(reshape(cc0, numel(idxy)*numel(idxx), []), [], 1);
    
    if 0

        figure; imagesc(mean(cc0(:,:,1:3),3))
        
        dv00 = gather_try(dv0);

        n_fr = 3;
        figure; hold on; axis equal tight;
        for n_fr = 1:numel(fi)
            imagesc(cc0(:,:,n_fr))
            plot(dv00(n_fr, 2)+1+lcorr, dv00(n_fr, 1)+1+lcorr, 'ro')
            pause(1/30)
            drawnow;
        end


        figure;  hold on; axis equal tight;
        imagesc(cc0_reg(:,:,n_fr))
        plot(dv00(n_fr, 2)+1+lcorr, dv00(n_fr, 1)+1+lcorr, 'ro')

        mov1 = cc0;
        mov1 = mov1 - min(mov1(:));
        mov1 = mov1/max(mov1(:));
        mov1 = mov1 * 2^12;

        f_save_mov_YS(mov1, 'F:\AC_data\caiman_data_missmatch\movies\M10_im1_test1.h5');

    end

end

corr_z = (corr - mean(corr_mean))/mean(corr_std);

if 0
    figure; plot(corr_z)
    figure; plot(corr)
    figure; plot(dv)

end

if nargin > 2 && removeMean
    dv = bsxfun(@minus, dv, mean(dv,1));
end


