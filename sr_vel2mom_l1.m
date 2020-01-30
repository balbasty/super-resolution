function [y,Dy] = sr_vel2mom_l1(y, prec, vs, w)
% Compute gradient of L1 regularisation operator.
%
%   In a reweighted least squares scheme, the L1 regulariser
%   [tr(sqrt(f'K'L^2'Kf))] is replaced by a reweighted L2 regulariser 
%   [tr(f'K'LWLKf)]. This function return the gradient of the reweighted 
%   regulariser [K'LWLKf] and the forward operator [LKf]
%
% FORMAT [KWKf,Kf] = sr_vel2mom_l1(f, prec, vs, w)
% f    - {nx ny nz}            - Input image (single channel)
% prec - {1}                   - Regularisation precision for each channel
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% vs   - {3}                   - Voxel size
% KWKf - {nx ny nz}            - Gradient of regulariser
% Kf   - {nx ny nz 3 2}        - Forward differential operator 

if nargout > 1
    Dy  = sqrt(prec) * sr_imgrad(y, vs);
end

spm_field('boundary', 1);
y = spm_field('vel2mom1', y, w, [vs 1], prec);
