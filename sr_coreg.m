function M = coreg(fref,varargin)
% Co-register a set of moving volumes with a fixed reference.
%
%   This function wraps `spm_coreg` and adds a few functionalities. It 
%   can take several moving images as input, and can register in a 
%   coarse-to-fine fashion by using less and less smoothing of the joint
%   histogram. Like `spm_coreg`, inputs must be nifti files (in-memory
%   arrays are not handled). However, orientation matrices can be
%   overriden.
%
% FORMAT T = utils.coreg(ref, mov1, mov2, ..., [fwhm])
% ref  - Filename or {Filename Matrix} of reference image
% mov  - Filename or {Filename Matrix} of moving image
% fwhm - Vector of decreasing FWHM for histogram smoothing [21 14 7]
% T    - 4x4xN transformation matrices.
%        Use M <- T\M to bring `ref` and `mov` in alignement.
%
% Note that `spm_vol` structures can be given as input instead of
% filenames.

% -------------------------------------------------------------------------
% FWHM
fwhms = [];
if ~isempty(varargin) && isnumeric(varargin{end})
    fwhms = varargin{end};
    varargin = varargin(1:end-1);
end
if isempty(fwhms)
    fwhms = [21 14 7];
end
fwhms = fwhms(:)';

% -------------------------------------------------------------------------
% REFERENCE
if iscell(fref)
    matref = fref{2};
    fref   = fref{1};
else
    matref = [];
end
ref = spm_vol(fref);
if ~isempty(matref)
    ref.mat = matref;
end

% -------------------------------------------------------------------------
% MOVING
M = repmat(eye(4), [1 1 numel(varargin)]);
for v=1:numel(varargin)
    fmov = varargin{v};
    if iscell(fmov)
        matmov = fmov{2};
        fmov   = fmov{1};
    else
        matmov = [];
    end
    mov = spm_vol(fmov);
    if isempty(matmov)
        mat = mov.mat;
    else
        mat = matmov;
    end
    for fwhm=fwhms
        mov.mat = M(:,:,v)\mat;
        q = spm_coreg(ref,mov,struct('fwhm',[fwhm fwhm]));
        M(:,:,v) = spm_matrix(q(:)') * M(:,:,v);
    end
end
    