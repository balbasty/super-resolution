function u = sr_uncertainty_l1(H, lam, vs, w)
% Diagonal uncertainty of reconstruction maps with MTV regularisation
%
% FORMAT u = sr_uncertainty_l1(H, prec, vs, w)
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% prec - {nf}                  - Regularisation for each channel (1/b^2)
% vs   - {3}                   - Voxel size
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% u    - {nx ny nz nf}         - Posterior uncertainty (= variance)

vs  = vs(:)';
lam = lam(:)';

spm_field('boundary', 1);
spm_diffeo('boundary', 1);
u   = zeros(size(H), 'single');
for c=1:size(H,4)
    u(:,:,:,c) = spm_field('diaginv1', H(:,:,:,c), w, [vs 1], lam(c));
end
