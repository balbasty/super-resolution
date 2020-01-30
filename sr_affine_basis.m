function B = sr_affine_basis(code)
% Generate a basis for the Lie algebra of affine matrices.
%
%   Affine transformation matrices are encode using their Lie algebra.
%   An affine matrix can then be generated by matrix-exponentiating its
%   projection on that basis, i.e.: M = expm(sum_k ak * Bk)
%
% FORMAT  = sr_affine_basis(group)
% group - Name of the group that must be encoded (D=0|2|3):
%         . T(D):  Translation
%         . SO(D): Special orthogonal (= translation + rotation)
%         . SE(D): Special Euclidean (= translation + rotation + scale)
% B     - 4x4xK Affine basis

% Author: John Ashburner

% This should probably be re-done to come up with a more systematic way 
% of defining the groups.
g     = regexpi(code,'(?<code>\w*)\((?<dim>\d?)\)','names');
if ~isfield(g, 'dim'), g.dim = '3'; end
g.dim = str2double(g.dim);
if (g.dim ~=0 && g.dim~=2 && g.dim~=3)
    error('Can not use size');
end
switch g.dim
case 0
    B = zeros(4,4,0);
case 2
    switch g.code
    case 'T' 
        B        = zeros(4,4,2);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
    case 'SO'
        B        = zeros(4,4,1);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
    case 'SE' 
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(1,2,3) =  1;
        B(2,1,3) = -1;
    otherwise
        error('Unknown group.');
    end
case 3
    switch g.code
    case 'T' 
        B        = zeros(4,4,3);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
    case 'SO' 
        B        = zeros(4,4,3);
        B(1,2,1) =  1;
        B(2,1,1) = -1;
        B(1,3,2) =  1;
        B(3,1,2) = -1;
        B(2,3,3) =  1;
        B(3,2,3) = -1;
    case 'SE' 
        B        = zeros(4,4,6);
        B(1,4,1) =  1;
        B(2,4,2) =  1;
        B(3,4,3) =  1;
        B(1,2,4) =  1;
        B(2,1,4) = -1;
        B(1,3,5) =  1;
        B(3,1,5) = -1;
        B(2,3,6) =  1;
        B(3,2,6) = -1;
    otherwise
        error('Unknown group.');
    end
end