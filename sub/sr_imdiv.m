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
vs   = sr_padarray(vs(:)', [0 max(0,3-numel(dim))], 'replicate', 'post');

% -------------------------------------------------------------------------
% Sum finite differences
D = zeros(dim,'like', X);
Z = sqrt(numel(type));
for i=1:numel(size(D))
    for t=1:numel(type)
        Xit = select_slices(X, [dimd dimt], {i t}) ./ (vs(i) * Z);
        switch type(t)
            case '+'
                D = D + cat(i,      -select_slices(Xit, i, 1), ...
                               -diff(select_slices(Xit, i, 1:(dim(i)-1)), 1, i), ...
                                     select_slices(Xit, i, dim(i)-1));
            case '-'
                D = D + cat(i,      -select_slices(Xit, i, 2), ...
                               -diff(select_slices(Xit, i, 2:dim(i)), 1, i), ...
                                     select_slices(Xit, i, dim(i)));
        end
    end
end

function S = select_slices(X, dim, ind)
% Select a (multidimensional) slice from a (multidimensional) array.
%
% FORMAT slice = utils.select_slice(array, dim, ind)
% array - {array}  input array
% dim   - {vector} dimension(s) along which to select a slice
% ind   - {[cell of] vector[s]} indices to select in each dimension
% slice - {array}  slice
%
% EXAMPLES
% >> % Create 5-dimensional array
% >> X = rand(10,10,10,10,10);
% >>
% >> % Select the 2nd slice along the 3rd dimension
% >> % Equivalent to: S = X(:,:,2,:,:)
% >> S = utils.select_slice(X, 3, 2);
% >>
% >> % Select the 1st and 2nd slice along the 1st dimension and the 5th
% >> % along the 4th dimension
% >> % Equivalent to: S = X([1 2],:,:,5,:)
% >> S = utils.select_slice(X, [1 4], {[1 2], 5});
sizeX    = size(X);
sub      = struct;
sub.type = '()';
sub.subs = repmat({':'}, [1 numel(sizeX)]);
if iscell(ind)
    assert(numel(dim) == numel(ind));
    for d=1:numel(dim)
        sub.subs{dim(d)} = ind{d};
    end
else
    sub.subs{dim} = ind;
end
S = subsref(X,sub);