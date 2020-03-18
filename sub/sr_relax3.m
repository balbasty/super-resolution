function [x,info] = sr_relax3(A,b,iE,x,opt)
% Relaxation solver (A = Block-diagonal + Band-diagonal).
%
%   Relax can be used to solve (large) linear systems of the form A*x = b, 
%   where A = E + F. E should be efficiently invertible (block-diagonal), 
%   and F should possess a local spatial structure (band-diagonal). 
%   Because of this local support, updates can be made inplace, following a 
%   checkerboard scheme. Updates take the form x = x + (E+s*I)\(b-A*x),
%   where `s` is a user-defined smoothing term.
%
% FORMAT [x,[info]] = sr_relax3(A,b,[iE],[x0],[opt])
%
% INPUT
% -----
% A       - Function handle for the linear operator `A*x = (E+F)*x`
% b       - Target array
% iE      - Function handle for the inverse operator `(E+s*I)\x`, where `s`
%           is a (user-defined) smoothing term.
% x0      - Initial guess [0]
% opt     - Structure of options with (optional) fields:
%           . nbiter    - Maximum number of iterations [32]
%           . tolerance - Tolerance for early stopping [1E-3]
%           . band      - Local (3d) support of the band matrix:
%                         ['checker'] or vector of bandwidth.
%           . verbose   - Verbosity level [0]
%           . sumtype   - Accumulator type 'native'/['double']
%
% OUTPUT
% ------
% x       - Solution to the linear system
% info    - Structure with fields:
%           . nbiter - Effective number of iterations
%           . rr     - Root mean squared residuals (per iteration)
%           . time   - Time (per iteration)
%
% NOTA BENE
% ---------
% . Function handles should be of the form f(x,[ind]), where `x` is the  
%   input vector and `ind` is a submask of voxels for which the problem 
%   must be solved.
%
% . Here, the spatial structure is assumed to be 3D, so 'vectors' should 
%   be 4D-arrays, where the first 3 dimensions correspond to 'spatial'
%   dimensions, and the 4-th dimension is a 'feature' dimension.
%
% . Note that summing in single is twice as fast as summing in double, but 
%   the additional precision is often needed.
%
% . More information on this scheme in:
%       Ashburner, J., 2007. A fast diffeomorphic image registration 
%       algorithm. Neuroimage, 38(1), pp.95-113.

DEBUG = false;

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------
if nargin < 4 || isempty(x) || (isscalar(x) && isnan(x))
    x = zeros(size(b), 'like', b);
end
if nargin < 5
    opt = struct;
end
if ~isfield(opt, 'nbiter'),    opt.nbiter    = 10;        end
if ~isfield(opt, 'tolerance'), opt.tolerance = 1E-3;      end
if ~isfield(opt, 'band'),      opt.band      = 'checker'; end
if ~isfield(opt, 'verbose'),   opt.verbose   = true;      end
if ~isfield(opt, 'sumtype'),   opt.sumtype   = 'double';  end
if isnumeric(opt.band)
    opt.band = sr_pad(opt.band(:)', [0 3-numel(opt.band)], 'replicate', 'post');
end

if DEBUG
    sr_plot_interactive('Init', 'Relaxation', 'Residuals', 'Iteration');
end

% -------------------------------------------------------------------------
% Initial residuals
% -------------------------------------------------------------------------
dim = [size(b) 1];
dim = dim(1:3);
bb  = sqrt(sum(b(:).^2, opt.sumtype));
r   = b - A(x);
rr  = sqrt(sum(r(:).^2, opt.sumtype));
rr  = rr/bb;
if opt.verbose
    fprintf('Relax: ');
end

if DEBUG
    ee   = A(x)-2*b;                                % A-norm of the error
    ee   = sum(ee(:).*x(:), opt.sumtype);           % (up to a constant)
    sr_plot_interactive('Set', 0, rr);
end

% -------------------------------------------------------------------------
% Checkerboard
% -------------------------------------------------------------------------
checker = ischar(opt.band) && strcmpi(opt.band, 'checker');
if checker
    checkerboard = reshape(1:(2*ceil(dim(1)/2)), [], 1, 1);
    checkerboard = bsxfun(@plus, checkerboard, reshape(1:(2*ceil(dim(2)/2)), 1, [], 1));
    checkerboard = bsxfun(@plus, checkerboard, reshape(1:(2*ceil(dim(3)/2)), 1, 1, []));
    checkerboard = rem(checkerboard, 2);
    checkerboard = checkerboard(1:dim(1),1:dim(2),1:dim(3));
end

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------
b = reshape(b, [], size(b,4));
start = tic;
time  = [];
for it=1:opt.nbiter

    if checker
    % ---------------------------------------------------------------------
    % Red and black scheme

        for j=[0 1]
            sub = (checkerboard == j);

            if false % (FASTER BUT NOT EXACT)
            % Compute residuals
            r = reshape(r, [], size(r,4));
            r(sub,:) = b(sub,:) - A(x,sub);
            r = reshape(r, [dim size(r,2)]);
            end

            % Solve easy (smoothed) inverse problem
            x = reshape(x, [], size(x,4));
            x(sub,:) = x(sub,:) + iE(r,sub);
            x = reshape(x, [dim size(x,2)]);
            
            if true % (SLOWER BUT EXACT)
            % Compute residuals
            r = reshape(b, [dim size(b,2)]) - A(x);
            end
        end

    else
    % ---------------------------------------------------------------------
    % 'band' scheme

        for k0=1:opt.band(3)
            k = k0:opt.band(3):dim(3);
            for j0=1:opt.band(2)
                j = j0:opt.band(2):dim(2);
                for i0=1:opt.band(1)
                    i = i0:opt.band(1):dim(1);

                    sub = sub2ind(dim,i,j,k);

                    % Compute residuals
                    r = reshape(r, [], size(r,4));
                    r(sub,:) = b(sub,:) - A(x,sub);
                    r = reshape(r, [dim size(r,2)]);

                    % Solve easy (smoothed) inverse problem
                    x = reshape(x, [], size(x,4));
                    x(sub,:) = x(sub,:) + iE(r,sub);
                    x = reshape(x, [dim size(x,2)]);

                end
            end

        end

    end

    % Compute residuals
    rr0  = rr(end);
    rr1  = sqrt(sum(r(:).^2, opt.sumtype));
    rr1  = rr1/bb;
    rr   = [rr rr1];
    if DEBUG
        time = [time toc(start)];
        e1   = reshape(A(x), [], size(b,2))-2*b;
        e1   = sum(e1(:).*x(:), opt.sumtype);
        ee   = [ee e1];
        sr_plot_interactive('Set', it, rr1);
    end
    if rr1 < opt.tolerance
        if opt.verbose, fprintf('o'); end
        break
    elseif rr1 > rr0
        if opt.verbose, fprintf('x'); end
        break;
    elseif opt.verbose
        fprintf('.');
    end
end
if opt.verbose
    fprintf('\n');
end

sr_plot_interactive('Clear');
info.nbiter = it;
info.rr     = rr;
if DEBUG
    info.time   = time;
    info.ee     = ee;
end