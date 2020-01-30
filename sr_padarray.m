function Y = sr_padarray(X, padsize, method, direction)
% FORMAT X = sr_padarray(X, padsize, [method], [direction])
% X         - An array
% padsize   - Padding size along each dimension of the array (>= 0)
% method    - 'circular', 'replicate', 'symmetric' or a value [0]
% direction - 'pre'/'post'/['both']
%
% Note that:
% . 'circular'  corresponds to the boundary condition of an FFT
% . 'symmetric' corresponds to the boundary condition of a DCT-II
%
% If padsize < 0, it is set to 0 instead.


% Possible extensions: 
% . method and direction could have a different value per dimension
% . other methods: antisymmetric (DST/Dirichlet)
%                  symmetry wrt voxel center rather than voxel side (DCT-I)

dim = size(X);
padsize = reshape(padsize, 1, []);

% Select appropriate indexing method
if ~isscalar(method)
    switch lower(method)
        case 'circular'
            fillidx = @circidx;
        case 'replicate'
            fillidx = @repidx;
        case 'symmetric'
            fillidx = @symidx;
        otherwise
            error('Unknoen method %s', method);
    end
end

% Output size a bigger for 'both' than for 'pre'/'post'
if strcmpi(direction, 'both')
    padfactor = 2;
else
    padfactor = 1;
end

% Compute output dimensions
dim     = [dim ones(1,max(numel(padsize)-numel(dim),0))];
padsize = [padsize ones(1,zeros(numel(dim)-numel(padsize),0))];
padsize(padsize < 0) = 0;
padsize(~isfinite(padsize)) = 0;
newdim  = dim + padsize * padfactor;

S = struct('type', '()', 'subs', {repmat({':'}, [1 numel(newdim)])});

if isscalar(method)
    % Allocate output volume
    if iscell(X)
        Y = cell(newdim);
        if isscalar(method)
            [Y{:}] = deal(method);
        end
    else
        if ~ischar(X)
            Y = method * ones(newdim, 'like', X);
        else
            Y = repmat(method, newdim);
        end
    end
    % Copy input volume
    switch lower(direction)
        case 'post'
            for d=1:numel(dim)
                S.subs{d} = 1:dim(d);
            end
        otherwise
            for d=1:numel(newdim)
                S.subs{d} = padsize(d)+(1:dim(d));
            end
    end
    Y = subsasgn(Y,S,X);
else
    % Sample input volume
    for d=1:numel(newdim)
        switch lower(direction)
            case 'post'
                S.subs{d} = fillidx(1:newdim(d), dim(d));
            otherwise
                S.subs{d} = fillidx((1:newdim(d)) - padsize(d), dim(d));
        end
    end
    Y = subsref(X,S);
end

% -------------------------------------------------------------------------
% FILLING FUNCTIONS

function i = circidx(i,n)

i          = i - 1;
nonneg     = (i >= 0);
i(nonneg)  = mod(i(nonneg), n);
i(~nonneg) = mod(n + mod(i(~nonneg), n), n);
i          = i + 1;

function i = repidx(i,n)

pre        = (i <= 0);
post       = (i > n);
i(pre)     = 1;
i(post)    = n;

function i = symidx(i,n)

i          = i - 1;
n2         = n*2;
pre        = (i < 0);
i(pre)     = n2 - 1 - mod(-i(pre)-1, n2);
i(~pre)    = mod(i(~pre), n2);
post       = (i >= n);
i(post)    = n2 - i(post) - 1;
i          = i + 1;
