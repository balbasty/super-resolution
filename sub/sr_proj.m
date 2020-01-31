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
info.vsx = sqrt(sum(xmat(1:3,1:3).^2));
info.vsy = sqrt(sum(ymat(1:3,1:3).^2));
% Template voxel size when projected in observed space (rotation only)
info.vsxy  = diag([info.vsx 1])*(xmat\ymat);
info.vsxy  = sqrt(sum(info.vsxy(1:3,1:3).^2));
info.matxy = info.matx/diag([info.vsx 1])*diag([info.vsxy 1]);
info.dimxy = ceil(xdim.*info.vsx./info.vsxy);
% Smoothing kernel
info.fwhm  = info.vsx ./ info.vsxy;
info.fwhm(info.fwhm <= 1) = 0;
info.fwhm(slinfo.profile == 0) = 0;         % Rect profiles: no smooting
info.fwhm = (1 - slinfo.gap) .* info.fwhm;  % If gap: smooth less
info.sd = info.fwhm/(2*sqrt(2*log(2)));     % FWHM to st. dev.
% Offset to add to projected template space
info.off  = -ceil(info.sd*4);
info.mato = eye(4);
info.mato(1:3,4) = info.off;
% Apply offset
info.dimxy = info.dimxy + 2*abs(info.off);
info.matxy = info.matxy*info.mato;
% -------------------------------------------------------------------------
function x = proj_A_sr(x, info)
w = warps_affine(info.dimxy, info.maty\info.matxy);
x = spm_diffeo('pull', x, w);               % Template to intermediate
spm_smooth(x, x, info.fwhm);                % Smooth intermediate
w = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('pull', x, w);               % Intermediate to observed
% -------------------------------------------------------------------------
function x = proj_At_sr(x, info)
w = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('push', x, w, info.dimxy);   % Observed to intermediate
spm_smooth(x, x, info.fwhm);                % Smooth intermediate
w = warps_affine(info.dimxy, info.maty\info.matxy);
x = spm_diffeo('push', x, w, info.dimy);    % Intermediate to template
% -------------------------------------------------------------------------
function x = proj_AtA_sr(x, info)
w1 = warps_affine(info.dimxy, info.maty\info.matxy);
w2 = warps_affine(info.dimx, info.matxy\info.matx);
x = spm_diffeo('pull', x, w1);              % Template to intermediate
spm_smooth(x, x, info.fwhm);                % Smooth intermediate
x = spm_diffeo('pull', x, w2);              % Intermediate to observed
x = spm_diffeo('push', x, w2, info.dimxy);  % Observed to intermediate
spm_smooth(x, x, info.fwhm);                % Smooth intermediate
x = spm_diffeo('push', x, w1, info.dimy);   % Intermediate to template

% === HELPERS =============================================================
function psi = warps_affine(lat, mat)
id  = warps_identity(lat);
psi = reshape(reshape(id,[prod(lat) 3])*mat(1:3,1:3)' + mat(1:3,4)',[lat 3]);
if lat(3) == 1, psi(:,:,:,3) = 1; end
% -------------------------------------------------------------------------
function id = warps_identity(d)
id = zeros([d(:)' 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
