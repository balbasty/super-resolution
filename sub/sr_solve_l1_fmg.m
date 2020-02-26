function dy = sr_solve_l1_fmg(H, g, w, lam, vs, opt)
% Relaxation solver for L1 spatial regularisation.
%
% FORMAT d = sr_solve_l1_fmg(H, g, w, lam, vs, [opt])
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% lam  - {nf}                  - Regularisation for each channel (1/b^2)
% vs   - {3}                   - Voxel size
% opt  - {struct}              - Structure of options with fields:
%      - nbcycle               - Number of full multigrid cycles [2]
%      . nbiter                - Number of relaxation iterations [2]
% d    - {nx ny nz nf}         - Step: d = H\g


if nargin < 6, opt = struct; end
if ~isfield(opt, 'nbcycle'),   opt.nbcycle   = 2;    end
if ~isfield(opt, 'nbiter'),    opt.nbiter    = 2;    end

% Neumann boundary conditon
spm_field('boundary', 1);
fmg = [opt.nbcycle opt.nbiter];
lam = lam(:)';

% We use a majoriser of the true Hessian
wbnd   = double(max(w(:)));
lambnd = lam * wbnd;
dy     = zeros(size(g), 'single');
for k=1:size(g,4)
    dy(:,:,:,k) = spm_field(H(:,:,:,k), g(:,:,:,k), [vs 0 1 0 fmg], lambnd(k));
end

