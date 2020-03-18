function [dy,info] = sr_solve_l1_relax(H, g, w, lam, vs, opt, dy)
% Relaxation solver for L1 spatial regularisation.
%
% FORMAT d = sr_solve_l1_relax(H, g, w, lam, vs, [opt], [d0])
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% lam  - {nf}                  - Regularisation for each channel (1/b^2)
% vs   - {3}                   - Voxel size
% opt  - {struct}              - Structure of options with fields:
%      . nbiter                - Number of CG iterations    [10]
%      . tolerance             - Early stopping tolerance   [1E-3]
%      . verbose               - Print progress             [true]
% d0   - {nx ny nz nf}         - Step: initial guess
% d    - {nx ny nz nf}         - Step: d = H\g


% Neumann boundary conditon
spm_field('boundary', 1);
fmg  = [2 2];
lam  = lam(:)';
wbnd = double(max(w(:)));

if all(w(:)==1)
    % We can use FMG
    lambnd = lam * wbnd;
    for k=1:size(g,4)
        dy(:,:,:,k) = spm_field(H(:,:,:,k), g(:,:,:,k), [vs 0 1 0 fmg], lambnd(k));
    end
    return
elseif nargin < 7 || isempty(dy)
    dy = zeros(size(g), 'single');
end

% Smoothing term: use diagonal of membrane kernel
scl = spm_diffeo('kernel', [3 3 3], [vs 0 1 0 0 0]);
scl = double(abs(scl(1,1,1)));
scl = scl*wbnd;

% Define forward and inverse operators
function Ax = A(x,ind)
% Forward operator: A = E + F
    spm_field('boundary', 1);
    Ax = spm_field('vel2mom1', x, w, [vs   1  ], lam);
    if nargin < 2
        Ax = Ax + H.*x;
    else
        Ax   = reshape(Ax, [], size(Ax,4));
        Ax   = Ax(ind,:);
        Ax   = reshape(Ax, [], 1, 1, size(Ax,2));
        Hsub = reshape(H, [], size(H,4));
        Hsub = Hsub(ind,:);
        Hsub = reshape(Hsub, [], 1, 1, size(Hsub,2));
        xsub = reshape(x, [], size(x,4));
        xsub = xsub(ind,:);
        xsub = reshape(xsub, [], 1, 1, size(xsub,2));
        Ax   = Ax + Hsub.*xsub;
        Ax   = reshape(Ax, [], size(Ax,4));
    end
end
function x = iE(x,ind)
% Forward operator: inv(E)
    spm_field('boundary', 1);
    if nargin < 2
        for kk=1:size(x,4)
            x(:,:,:,kk) = spm_field(H(:,:,:,kk) , x(:,:,:,kk) , [1 1 1 scl 0 0], lam(kk));
        end
    else
        Hsub = reshape(H, [], size(H,4));
        Hsub = Hsub(ind,:);
        Hsub = reshape(Hsub, [], 1, 1, size(Hsub,2));
        x = reshape(x, [], size(x,4));
        x = x(ind,:);
        x = reshape(x, [], 1, 1, size(x,2));
        for kk=1:size(x,4)
            x(:,:,:,kk) = spm_field(Hsub(:,:,:,kk) , x(:,:,:,kk) , [1 1 1 scl 0 0], lam(kk));
        end
        x = reshape(x, [], size(x,4));
    end

end

% Relax
[dy,info] = sr_relax3(@A, g, @iE, dy, opt);

end

