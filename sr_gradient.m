function [llx,g,H] = sr_gradient(c, in0, out, opt)
% FORMAT [llx,g,H] = sr_gradient(in, out, [opt])
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
ydim = out.dim;   % Model dimensions
ymat = out.mat;  % Model orientation matrix
y0   = single(out.dat(:,:,:,c));
if opt.log, y0 = exp(y0); end
llx  = 0;
if nargout > 1
    g = zeros(ydim, 'single');
    if nargout > 2
        H = zeros(ydim, 'single');
    end
end

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
    % Compute residuals
    x   = single(in.dat());
    y   = sr_proj('A', y0, xdim, xmat, ymat, opt, slice); % Project to low-res
    msk = isfinite(x) & x > 0 & isfinite(y);
    y(~msk) = 0;
    x(~msk) = 0;
    clear msk;
    res = y - x;
    clear x y

    % ---------------------------------------------------------------------
    % Compute log-likelihood
    llx = llx + 0.5 * lam * sum(res(:).^2, 'double');

    % ---------------------------------------------------------------------
    % Compute gradient and Hessian in observed space
    if nargout > 1
        g1 = sr_proj('At', res, ydim, xmat, ymat, opt, slice); clear res
        if opt.log, g1 = g1 .* y0; end
        g  = g + lam * g1; clear g1
        if nargout > 2
            H1 = sr_proj('AtA', ones(ydim, 'single'), xdim, xmat, ymat, opt, slice);
            if opt.log, H1 = H1 .* y0.^2; end
            H  = H + lam * H1; clear H1
        end
    end
end