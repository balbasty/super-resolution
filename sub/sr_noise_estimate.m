function [noise,mu_val,info] = sr_noise_estimate(Scans,K)
% Estimate average noise from a series of images
% FORMAT [noise,mu_val] = sr_noise_estimate(Scans)
% Scans  - nifti structures or filenames of images
% K      - Number of Rician mixture components
% noise  - standard deviation estimate
% mu_val - expectation of more intense Rician
% info   - This struct can be used for plotting the fit as:
%              plot(info.x(:),info.p,'--',info.x(:), ...
%                   info.h/sum(info.h)/info.md,'b.', ...
%                   info.x(:),info.sp,'r');
% _________________________________________________________________________
%  Copyright (C) 2012-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% Yael Balbastre: make it work for any input type, not just nifti.
% $Id: spm_noise_estimate.m 7599 2019-05-30 13:50:41Z mikael $

if nargin < 2, K = 2; end
getmu = nargout >= 2;

% - Compute number of 3D scans + convert to array/file_array
if ~iscell(Scans), Scans = {Scans}; end
N = 0;
for i=1:numel(Scans)
    if ischar(Scans{i})
        Scans{i} = nifti(Scans{i});
    end
    if isa(Scans{i}, 'nifti')
        Scans{i} = Scans{i}.dat;
    end
    dim = [size(Scans{i}) 1];
    Scans{i} = reshape(Scans{i}, dim(1), dim(2), dim(3), prod(dim(4:end)));
    N = N + size(Scans{i}, 4);
end

% - Estimate noise
noise  = zeros(N,1);
mu_val = zeros(N,1);
info   = struct('x',[],'h',[],'p',[],'sp',[],'md',[]);
n      = 0;
for i=1:numel(Scans)
    if isa(Scans{i}, 'file_array')
        inttype = spm_type(Scans{i}.dtype(1:(end-3)),'intt');
        slope   = Scans{i}.scl_slope;
    else
        inttype = isinteger(Scans{i});
        slope   = 1;
    end
    for k=1:size(Scans{i}, 4)
        n = n + 1;
        f = Scans{i}(:,:,:,k);
        [noise(n), mu_val(n), info(n)] = estimate_noise_1(f, K, inttype, slope, getmu);
    end
end

end

% Helper: estimate noise from one 3D scan
function [noise, mu_val, info] = estimate_noise_1(f, K, inttype, slope, getmu)

if inttype
    f(f==max(f(:))) = 0;
    x      = double(0:slope:max(f(:)));
    [h,x]  = hist(f(f~=0),x);
else
    x      = (0:1023)*(max(f(:))/1023);
    f(f==max(f(:))) = 0;
    [h,x]  = hist(f(f~=0 & isfinite(f)),x);
end
[mg,nu,sd,info] = spm_rice_mixture(double(h(:)),double(x(:)),K);
noise           = min(sd);

mu_val = NaN;
if getmu && nargout>=2
    x          = -nu.^2./(2*sd.^2);
    msk        = x>-20;
    Laguerre   = exp(x(msk)/2).*((1-x(msk)).*besseli(0,-x(msk)/2)-x(msk).*besseli(1,-x(msk)/2));
    Ey         = zeros(size(sd));
    Ey( msk)   = sqrt(pi*sd(msk).^2/2).*Laguerre;
    Ey(~msk)   = nu(~msk);
    mu_val     = max(Ey);
end
    
end