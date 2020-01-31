function [lat,mat,vs] = sr_mean_space(datmat, datlat, vs, crop)
% Compute the (mean) model space from individual spaces.
%
%   This procedure relies on finding a barycentre of the individual
%   orientation matrices, as described in:
%       Pennec & Arsigny, "Exponential Barycenters of the Canonical  
%       Cartan Connection and Invariant Means on Lie Groups."
%       Matrix Information Geometry, Springer, pp.123-168, 2012.
%       https://dx.doi.org/10.1007/978-3-642-30232-9_7
%       https://hal.inria.fr/hal-00699361
%   Beforehand, additional 90 degree rotations and flips are removed from
%   the individual orientation matrices.
%   Afterwards, remaining shears are removed from the mean orientation
%   matrix.
%
% FORMAT [lat,mat,vs] = sr_mean_space(datmat, datlat, vs, crop)
% datmat - 4x4xN - Subjects orientation matrices 
% datlat - 3xN   - Subjects lattice dimensions
% vs     - 3     - Model space voxel size [NaN = mean from subjects]
% crop   - 3     - Percentage to crop along each dimension [0]

% Author: John Ashburner

if nargin < 4, crop = 0;   end
if nargin < 3, vs   = NaN; end

crop = sr_padarray(crop(:), [3-numel(crop) 0], 'replicate', 'post');
vs   = sr_padarray(vs(:),   [3-numel(vs)   0], 'replicate', 'post');

% Rigid-body matrices computed from exp(p(1)*B(:,:,1)+p(2)+B(:,:,2)...)
%--------------------------------------------------------------------------
is2d = datlat(3,1) == 1;
if is2d, B = sr_affine_basis('SE(2)');
else,    B = sr_affine_basis('SE(3)'); end

% Find combination of 90 degree rotations and flips that brings all
% the matrices closest to axial
%--------------------------------------------------------------------------
datmat0 = datmat;
pmatrix = [1, 2, 3;
           2, 1, 3;
           3, 1, 2;
           3, 2, 1;
           1, 3, 2;
           2, 3, 1];
for n=1:size(datmat,3)
    vs1   = sqrt(sum(datmat(1:3,1:3,n).^2));
    R     = datmat(:,:,n)/diag([vs1 1]);
    R     = R(1:3,1:3);
    minss = Inf;
    minR  = eye(3);
    for i=1:6 % Permute (= 'rotate + flip') axes
        R1 = zeros(3);
        R1(pmatrix(i,1),1)=1;
        R1(pmatrix(i,2),2)=1;
        R1(pmatrix(i,3),3)=1;
        for j=0:7 % Mirror (= 'flip') axes
            F  = diag([bitand(j,1)*2-1, bitand(j,2)-1, bitand(j,4)/2-1]);
            R2 = F*R1;
            ss = sum(sum((R/R2-eye(3)).^2));
            if ss<minss
                minss = ss;
                minR  = R2;
            end
        end
    end
    rdim = abs(minR*datlat(:,n));
    R2   = inv(minR);
    minR = [R2 R2*((sum(R2,1)'-1)/2.*(rdim+1)); 0 0 0 1];
    datmat(:,:,n) = datmat(:,:,n)*minR;
end

% Average of these matrices
%--------------------------------------------------------------------------
mat = spm_meanm(datmat);

% If average involves shears, then find the closest matrix that does not
% require them
%--------------------------------------------------------------------------
p = spm_imatrix(mat);
if sum(p(10:12).^2)>1e-8

    % Zooms computed from exp(p(7)*B2(:,:,1)+p(8)*B2(:,:,2)+p(9)*B2(:,:,3))
    %----------------------------------------------------------------------
    B2        = zeros(4,4,3);
    B2(1,1,1) = 1;
    B2(2,2,2) = 1;
    B2(3,3,3) = 1;

    p      = zeros(9,1); % Parameters
    for it=1:10000
        [R,dR] = spm_dexpm(p(1:6),B);  % Rotations + Translations
        [Z,dZ] = spm_dexpm(p(7:9),B2); % Zooms

        M  = R*Z; % Voxel-to-world estimate
        dM = zeros(4,4,6);
        for n=1:6, dM(:,:,n)   = dR(:,:,n)*Z; end
        for n=1:3, dM(:,:,n+6) = R*dZ(:,:,n); end
        dM = reshape(dM,[16,9]);

        d   = M(:)-mat(:);   % Difference
        gr  = dM'*d;         % Gradient
        Hes = dM'*dM;        % Hessian
        p   = p - Hes\gr;    % Gauss-Newton update
        if sum(gr.^2)<1e-8, break; end
    end
    mat = M;
end

% Set required voxel size
%--------------------------------------------------------------------------
vs0 = vs;
vs  = sqrt(sum(mat(1:3,1:3).^2));
vs0(~isfinite(vs0)) = vs(~isfinite(vs0));
mat = mat * diag([vs0(:)./vs(:); 1]);
vs  = sqrt(sum(mat(1:3,1:3).^2));

% Ensure that the FoV covers all images, with a few voxels to spare
%--------------------------------------------------------------------------
mn    =  Inf*ones(3,1);
mx    = -Inf*ones(3,1);
for n=1:size(datmat0,3)
    dm      = datlat(:,n);
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1];
    M   = mat\datmat0(:,:,n);
    vs1 = M(1:3,:)*corners;
    mx  = max(mx,max(vs1,[],2));
    mn  = min(mn,min(vs1,[],2));
end
mx    = ceil(mx);
mn    = floor(mn);
o     = -crop.*(mx - mn);
if is2d, o(1:3) = 0; end
lat = (mx - mn + (2*o + 1))';
mat = mat * [eye(3) mn - (o + 1); 0 0 0 1];
end