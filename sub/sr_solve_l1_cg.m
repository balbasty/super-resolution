function dy = sr_solve_l1_cg(H, g, w, lam, vs, opt, dy)
% Conjugate gradient solver for L1 spatial regularisation.
%
% FORMAT d = sr_solve_l1_cg(H, g, w, lam, vs, [opt], [dy0])
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% lam  - {nf}                  - Regularisation value for each channel (1/b^2)
% vs   - {3}                   - Voxel size
% opt  - {struct}              - Structure of options with fields:
%      . nbiter                - Number of CG iterations    [10]
%      . tolerance             - Early stopping tolerance   [1E-3]
%      . verbose               - Print progress             [true]
%      . precond               - Use preconditioner         [false]
% d0   - {nx ny nz nf}         - Step: initial guess
% d    - {nx ny nz nf}         - Step: d = H\g

if nargin < 6, opt = struct; end
if ~isfield(opt, 'nbiter'),    opt.nbiter    = 10;    end
if ~isfield(opt, 'tolerance'), opt.tolerance = 1E-3;  end
if ~isfield(opt, 'verbose'),   opt.verbose   = true;  end
if ~isfield(opt, 'precond'),   opt.precond   = false;  end

% Neumann boundary conditon
spm_field('boundary', 1);
fmg    = [2 2];
lam    = lam(:)';
wbnd   = double(max(w(:)));
lambnd = lam * wbnd;

if all(w(:)==1)
    % We can use FMG
    for k=1:size(g,4)
        dy(:,:,:,k) = spm_field(H(:,:,:,k), g(:,:,:,k), [vs 0 1 0 fmg], lambnd(k));
    end
    return
elseif nargin < 7 || isempty(dy)
    dy = zeros(size(g), 'single');
end

% Prior term
function y = prior(x)
    y = spm_field('vel2mom1', x, w, [vs   1  ], lam);
end

% Solve inversion using conjugate-gradient
iM = @(b) spm_field(H, b, [vs 0 1 0 fmg], lambnd);
HH = @(b) spm_field('Atimesp', H, b) + prior(b);
if opt.precond, precond = iM;
else,           precond = []; end
dy = sr_cg(HH, g, dy, precond, opt.nbiter, opt.tolerance, opt.verbose);

end