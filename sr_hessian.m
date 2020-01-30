function H = sr_hessian(c, in0, out, opt)
% FORMAT H = sr_hessian(c, in, out, [opt])
% in  - Input data structure
% out - Model data structure 
% opt - Structure of parameters with [optional] fields:
%       . subsample - Subsampling distance in mm [Inf=no]
%       . verbose   - Verbosity level [0]

% -------------------------------------------------------------------------
% Options
if nargin < 4, opt = struct; end
if ~isfield(opt, 'verbose'),   opt.verbose   = 0;     end

% -------------------------------------------------------------------------
% Read model info
ydim = out.dim;  % Model dimensions
ymat = out.mat;  % Model orientation matrix
if opt.log
    y0   = exp(single(out.dat(:,:,:,c)));
end
H    = zeros(ydim, 'single');

% -------------------------------------------------------------------------
% Iterate over repeats
for r=1:numel(in0)
    if opt.verbose > 0, fprintf('.'); end

    in    = in0{r};
    if isfield(in, 'slice'), slice = in.slice; else, slice = []; end
    xdim  = in.dim;      % Observed lattice
    xmat  = in.mat;      % Observed orientation matrix
    lam   = in.lam;      % Noise precision

    % ---------------------------------------------------------------------
    % Compute gradient and Hessian in observed space
    H1 = sr_proj('AtA', ones(ydim, 'single'), xdim, xmat, ymat, opt, slice);
    if opt.log, H1 = H1 .* y0.^2; end
    H  = H + lam * H1; clear H1
end