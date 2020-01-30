function dy = sr_solve_l2(H, g, prec, vs)
% Solver for L2 spatial regularisation.
%
% FORMAT d = solve(H, g, prec, vs)
% H    - {nx ny nz nf} - Field of Hessian matrices (diagonal)
% g    - {nx ny nz nf} - Field of gradients
% prec - {nf}          - Regularisation precision for each channel
% vs   - {3}           - Voxel size
% d    - {nx ny nz nf} - Step: d = H\g

fmg = [2 2];
spm_field('boundary', 1);
dy = zeros(size(g), 'single');
for k=1:size(g,4)
    dy(:,:,:,k) = spm_field(H(:,:,:,k), g(:,:,:,k), [vs 0 1 0 fmg], prec(k));
end