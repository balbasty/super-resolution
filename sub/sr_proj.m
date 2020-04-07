function out = sr_proj(mode, in, odim, xmat, ymat, opt, slinfo)
% Apply projection matrix (or its adjoint)
%
% FORMAT sr_proj('A',   in, xdim, xmat, ymat, opt, slinfo) [FORWARD]
% FORMAT sr_proj('At',  in, ydim, xmat, ymat, opt, slinfo) [ADJOINT]
% FORMAT sr_proj('AtA', in, xdim, xmat, ymat, opt, slinfo) [ADJxFOR]

spm_diffeo('boundary', 1);

switch lower(opt.mode(1))
    case 'd' % Denoise
        idim = [size(in) 1];
        idim = idim(1:3);
        if all(all(xmat == ymat)) && all(odim == idim)
            % Projection = identity
            out = in;
        else
            % Need to use pull/push
            switch lower(mode)
                case 'a'
                    out = proj_A_denoise(in, odim, xmat, ymat);
                case 'at'
                    out = proj_At_denoise(in, odim, xmat, ymat);
                case'ata'
                    out = proj_AtA_denoise(in, odim, xmat, ymat);
            end
        end
    case 's' % Super-resolution
        switch lower(mode)
            case 'a'
                info = proj_info_sr(odim, xmat, size(in), ymat, slinfo);
                out = proj_A_sr(in, info);
            case 'at'
                info = proj_info_sr(size(in), xmat, odim, ymat, slinfo);
                out = proj_At_sr(in, info);
            case'ata'
                info = proj_info_sr(odim, xmat, size(in), ymat, slinfo);
                out = proj_AtA_sr(in, info);
        end
end

% === DENOISE =============================================================
function out = proj_A_denoise(in, xdim, xmat, ymat)
w   = warps_affine(xdim, ymat\xmat);
out = spm_diffeo('pull', in, w);
% -------------------------------------------------------------------------
function out = proj_At_denoise(in, ydim, xmat, ymat)
xdim = [size(in) 1];
xdim = xdim(1:3);
w    = warps_affine(xdim, ymat\xmat);
out  = spm_diffeo('push', in, w, ydim);
% -------------------------------------------------------------------------
function out = proj_AtA_denoise(in, xdim, xmat, ymat)
ydim = [size(in) 1];
ydim = ydim(1:3);
w    = warps_affine(xdim, ymat\xmat);
out  = spm_diffeo('pull', in, w);
out  = spm_diffeo('push', out, w, ydim);
% -------------------------------------------------------------------------

% === SUPER-RES ===========================================================
function info = proj_info_sr(xdim, xmat, ydim, ymat, slinfo)
info = struct;
% Voxel sizes
info.dimx = [xdim 1];
info.dimx = info.dimx(1:3);
info.dimy = [ydim 1];
info.dimy = info.dimy(1:3);
info.matx = xmat;
info.maty = ymat;
info.vsx  = sqrt(sum(xmat(1:3,1:3).^2));
info.vsy  = sqrt(sum(ymat(1:3,1:3).^2));
% Template voxel size when projected in observed space (rotation only)
info.vsxy  = diag([info.vsx 1])*(xmat\ymat);
info.vsxy  = sqrt(sum(info.vsxy(1:3,1:3).^2));
info.matxy = info.matx/diag([info.vsx 1])*diag([info.vsxy 1]);
info.dimxy = ceil(xdim.*info.vsx./info.vsxy);
% Smoothing kernel
info.fwhm  = info.vsx ./ info.vsxy;
info.fwhm  = (1 - slinfo.gap) .* info.fwhm;    % If gap: smooth less
info.sd    = info.fwhm/(2*sqrt(2*log(2)));     % FWHM to st. dev.
info.shape = slinfo.profile;                   % Pulse shape
info.acc   = slinfo.accumulate;                % Accumulate or Average
% Offset to add to projected template space
info.off  = -ceil(info.sd*6);
info.mato = eye(4);
info.mato(1:3,4) = info.off;
% Apply offset
info.dimxy = info.dimxy + 2*abs(info.off);
info.matxy = info.matxy*info.mato;
% -------------------------------------------------------------------------
function x = proj_A_sr(x, info)
w = warps_affine(info.dimxy, info.maty\info.matxy);
x = spm_diffeo('pullc', x, w);                    % Template to intermediate
% spm_smooth(x, x, info.fwhm);                     % Smooth intermediate
x = smooth(x, info.fwhm, info.shape, info.acc);  % Smooth intermediate
w = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('pull', x, w);                    % Intermediate to observed
% -------------------------------------------------------------------------
function x = proj_At_sr(x, info)
w = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('push', x, w, info.dimxy);        % Observed to intermediate
% spm_smooth(x, x, info.fwhm);                     % Smooth intermediate
x = smooth(x, info.fwhm, info.shape, info.acc);  % Smooth intermediate
w = warps_affine(info.dimxy, info.maty\info.matxy);
x = spm_diffeo('pushc', x, w, info.dimy);         % Intermediate to template
% -------------------------------------------------------------------------
function x = proj_AtA_sr(x, info)
w1 = warps_affine(info.dimxy, info.maty\info.matxy);
w2 = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('pullc', x, w1);                  % Template to intermediate
% spm_smooth(x, x, info.fwhm);                    % Smooth intermediate
x = smooth(x, info.fwhm, info.shape, info.acc); % Smooth intermediate
x = spm_diffeo('pull', x, w2);                  % Intermediate to observed
x = spm_diffeo('push', x, w2, info.dimxy);      % Observed to intermediate
% spm_smooth(x, x, info.fwhm);                    % Smooth intermediate
x = smooth(x, info.fwhm, info.shape, info.acc); % Smooth intermediate
x = spm_diffeo('pushc', x, w1, info.dimy);       % Intermediate to template

% === HELPERS =============================================================
function psi = warps_affine(lat, mat)
id  = warps_identity(lat);
psi = reshape(reshape(id,[prod(lat) 3])*mat(1:3,1:3)' + mat(1:3,4)',[lat 3]);
if lat(3) == 1, psi(:,:,:,3) = 1; end
% -------------------------------------------------------------------------
function id = warps_identity(d)
id = zeros([d(:)' 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
% -------------------------------------------------------------------------
function x = smooth(x, width, shape, mode)
kern = cell(1,3);
for d=1:3
    kern{d} = smoothkern(width(d), shape(d), mode, 1);
end
i  = (length(kern{1}) - 1)/2;
j  = (length(kern{2}) - 1)/2;
k  = (length(kern{3}) - 1)/2;
spm_conv_vol(x, x, kern{1}, kern{2}, kern{3}, -[i,j,k]);
% -------------------------------------------------------------------------
function ker = smoothkern(w, shape, acc, basis)
% width - FWHM: ratio between high-res and low-res voxel sizes
% shape - Shape of the slice profile (0=rect, 1=tri, 2=gauss)
% acc   - Average (f) or accumulate (t) intensities
% basis - Basis used to encode the high-res image (0=nn, 1=tri)
%
% The convolution kernel is obtained by convolving analytically the
% slice profile with the image basis function. See spm_smoothkern.
%
% `acc` decides if the integral of the kernel should be 1 or `w`
% (it might not be optimal for the Gaussian profile, as the non-normalised
% integral is w*sqrt(pi/(2*log(2))). I need to check how Gaussian slice 
% profiles really integrate.).
% A (better?) solution could be to optimise the intensity scaling
% parameter.
if w == 0, ker = 1; return; end
switch shape
    case 0
        % RECTANGLE
        if basis
            lim = ceil((w+2)/2);
            x   = -lim:lim;
            neg_low  = min(max(x-w/2,-1), 0);
            neg_high = max(min(x+w/2, 0),-1);
            pos_low  = min(max(x-w/2, 0), 1);
            pos_high = max(min(x+w/2, 1), 0);
            ker = integrate_poly(neg_low, neg_high, 1,  1) ...
                + integrate_poly(pos_low, pos_high, 1, -1);
        else
            lim = ceil((w+1)/2);
            x   = -lim:lim;
            ker = max(min(x+0.5,w/2) - max(x-0.5,-w/2), 0);
        end
        if ~acc, ker = ker/w; end
    case 1
        % TRIANGLE
        if basis
            lim = ceil((2*w+2)/2);
            x = -lim:lim;
            neg_neg_low  = min(max(x,  -1), 0);
            neg_neg_high = max(min(x+w, 0),-1);
            neg_pos_low  = min(max(x,   0), 1);
            neg_pos_high = max(min(x+w, 1), 0);
            pos_neg_low  = min(max(x-w,-1), 0);
            pos_neg_high = max(min(x,   0),-1);
            pos_pos_low  = min(max(x-w, 0), 1);
            pos_pos_high = max(min(x,   1), 0);
            ker = integrate_poly(neg_neg_low, neg_neg_high, 1+x/w,  1+x/w-1/w, -1/w) ...
                + integrate_poly(neg_pos_low, neg_pos_high, 1+x/w, -1-x/w-1/w,  1/w) ...
                + integrate_poly(pos_neg_low, pos_neg_high, 1-x/w,  1-x/w+1/w,  1/w) ...
                + integrate_poly(pos_pos_low, pos_pos_high, 1-x/w, -1+x/w+1/w, -1/w);
        else
            lim = ceil((2*w+1)/2);
            x   = -lim:lim;
            neg_low  = min(max(x-0.5,-w), 0);
            neg_high = max(min(x+0.5, 0),-w);
            pos_low  = min(max(x-0.5, 0), w);
            pos_high = max(min(x+0.5, w), 0);
            ker = integrate_poly(neg_low, neg_high, 1,  1/w) ...
                + integrate_poly(pos_low, pos_high, 1, -1/w);
        end
        if ~acc, ker = ker/w; end
    case {2 'g'}
        % GAUSSIAN
        lim = ceil((4*w+1+basis)/2);
        ker = spm_smoothkern(w, -lim:lim, basis);
        if acc, ker = ker*w; end
    case 's'
        % SINC
        if basis
            error('Sinc*Tri not implemented yet.');
        else
            lim = ceil((4*w+1)/2);
            x = -lim:lim;
            ker = sinint(pi*(x/w+0.5)) - sinint(pi*(x/w-0.5));
        end
    otherwise
        error('Unknown shape parameter: should be 0, 1  or 2');
end
% -------------------------------------------------------------------------
function k = integrate_poly(l, h, varargin)
% FORMAT k = integrate_poly(l, h, a, b, c, ...)
% Integrate polynomial a+b*x+c*x^2+... on [l,h]
k  = 0;
hh = h;
ll = l;
for i=1:numel(varargin)
    if varargin{i} ~= 0
        k = k + (varargin{i}/i).*(hh-ll);
    end
    hh = hh.*h;
    ll = ll.*l;
end