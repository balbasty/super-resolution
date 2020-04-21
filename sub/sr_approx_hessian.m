function H = sr_approx_hessian(ydim, xdim, xmat, ymat, opt, slice)
% FORMAT H = sr_approx_hessian(ydim, xdim, xmat, ymat, opt, slice)

% -------------------------------------------------------------------------
% Compute approximate Hessian
H = sr_proj('AtA', ones(ydim, 'single'), xdim, xmat, ymat, opt, slice);