function D = sr_imdiv(X,vs,type)  
% Computes the divergence of an image (with voxel size)
% FORMAT div = utils.imdiv(img,(vs),(type))
% img   - Image in "gradient space" (Nx * Nz * Nz * ... * Ndim * Ntype)
% vs    - Voxel size [1]
% type  - Finite difference type '+'/'-'/['+-']/'-+'
% div   - Divergence (Nx * Nz * Nz * ...)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% -------------------------------------------------------------------------
% Default parameters
if nargin < 2, vs = ones(1,'like',X); end
if nargin < 3, type = '+-'; end

dim = size(X);
if numel(type) > 1
    dimt = numel(dim);
    dim = dim(1:end-1);
else
    dimt = numel(dim) + 1;
end
dimd = numel(dim);
dim  = dim(1:end-1);
vs   = padarray(vs(:)', [0 max(0,3-numel(dim))], 'replicate', 'post');

% -------------------------------------------------------------------------
% Sum finite differences
D = zeros(dim,'like', X);
Z = sqrt(numel(type));
for i=1:numel(size(D))
    for t=1:numel(type)
        Xit = select_slices(X, [dimd dimt], {i t}) ./ (vs(i) * Z);
        switch type(t)
            case '+'
                D = D + cat(i,      -utils.select_slices(Xit, i, 1), ...
                               -diff(utils.select_slices(Xit, i, 1:(dim(i)-1)), 1, i), ...
                                     utils.select_slices(Xit, i, dim(i)-1));
            case '-'
                D = D + cat(i,      -utils.select_slices(Xit, i, 2), ...
                               -diff(utils.select_slices(Xit, i, 2:dim(i)), 1, i), ...
                                     utils.select_slices(Xit, i, dim(i)));
        end
    end
end