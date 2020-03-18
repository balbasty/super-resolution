function dy = sr_solve_l0(H,g)
% Solver for L0 spatial regularisation.
%
% FORMAT d = sr_solve_l0(H,g)
% H - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g - {nx ny nz nf}         - Field of gradients
% d - {nx ny nz nf}         - Step: d = H\g

spm_field('boundary', 1);

dy = zeros(size(g), 'single');
for k=1:size(g,4)
    dy(:,:,:,k) = spm_field(H(:,:,:,k), g(:,:,:,k));
end
