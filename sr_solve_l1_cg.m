function dy = sr_solve_l1_cg(H, g, w, lam, vs, opt)
% Conjugate gradient solver for L1 spatial regularisation.
%
% FORMAT d = sr_solve_l1_cg(H, g, w, lam, vs, [opt])
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% lam  - {nf}                  - Regularisation value for each channel (1/b^2)
% vs   - {3}                   - Voxel size
% opt  - {struct}              - Structure of options with fields:
%      . nbiter                - Number of CG iterations [10]
%      . tolerance             - Early stopping tolerance [1E-3]
%      . verbose               - Print progress [true]
%      . precond               - Use preconditioner [true]
% d    - {nx ny nz nf}         - Step: d = H\g

if nargin < 6
    opt = struct;
end
if ~isfield(opt, 'nbiter'),    opt.nbiter    = 10;    end
if ~isfield(opt, 'tolerance'), opt.tolerance = 1E-3;  end
if ~isfield(opt, 'verbose'),   opt.verbose   = true;  end
if ~isfield(opt, 'precond'),   opt.precond   = true;  end

% Neumann boundary conditon
spm_field('boundary', 1);
fmg = [2 2];
lam = lam(:)';

% Initial guess using a majoriser of the true Hessian
wbnd   = double(max(w(:)));
lambnd = lam * wbnd;
dy     = spm_field(H, g, [vs 0 1 0 fmg], lambnd);

if all(w(:)==1)
    % No need for CG
    return
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
    
%     % Matrix-vector product for a diagonal/symmetric Matrix
%     function y = matvec(A,x)
%         y = zeros(size(x), 'like', x);
%         K = size(x,4);
%         ind = symIndices(size(x,4));
%         for n1=1:K
%             ind1 = ind(n1,:);
%             y(:,:,:,n1) = y(:,:,:,n1) + A(:,:,:,ind1(n1)) .* x(:,:,:,n1);
%             for n2=(n1+1):K
%                 ind2 = ind1(n2);
%                 Ak = A(:,:,:,ind2);
%                 y(:,:,:,n1) = y(:,:,:,n1) + Ak .* x(:,:,:,n2);
%                 y(:,:,:,n2) = y(:,:,:,n2) + Ak .* x(:,:,:,n1);
%             end
%         end
%     end
