function o = sr_opt_defaults(o)
% Set default options.
%
% FORMAT opt_out = sr_opt_defaults(opt_in)
% opt_in  - User-defined option structure [empty]
% opt_out - All user-defined options are kept untouched, the other are set 
%           to their default value.

% --- Default values
if nargin < 1, o = struct; end
o = setdefault(o, 'mode',             'denoise');        % Mode (denoise|superres)
o = setdefault(o, 'log',              false);            % Solve for log-images
o = setdefault(o, 'itermax',          20);               % Maximum number of iterations
o = setdefault(o, 'tolerance',        1E-4);             % Gain threshold for early stopping (per scale)
o = setdefault(o, 'out.folder',       '.');              % Output folder [same as input if empty]
o = setdefault(o, 'out.mem',          'map');            % map/load output volumes
o = setdefault(o, 'reg.mode',         1);                % Regularisation mode (0=None|1=L1|2=L2)
o = setdefault(o, 'reg.value',        1E2);              % Regularisation value (L2: 1/sig^2, L1: 1/b^2)
o = setdefault(o, 'reg.smo' ,         1E-3);             % RLS smoother
o = setdefault(o, 'vs',               NaN);              % Reconstruction voxel size (Nan=from input)
o = setdefault(o, 'fov',              0);                % Field of view (0=bounding box|n=index of input volume)
o = setdefault(o, 'coreg.do',         true);             % Co-register volumes first
o = setdefault(o, 'coreg.fwhm',       [21 14 7]);        % Co-register volumes first
o = setdefault(o, 'threads',          'automatic');      % Number of threads (automatic|integer)
o = setdefault(o, 'verbose',          1);                % Verbosity (0=quiet|[1]=print|2=plot)
o = setdefault(o, 'input.mat',        eye(4));           % Default orientation matrix
o = setdefault(o, 'slice.dir',        'thickest');       % Default slice direction (thickest|all)
o = setdefault(o, 'slice.gap',        1/3);              % Default slice gap
o = setdefault(o, 'armijo',           [2 1]);            % Default slice gap

% o = setdefault(o, 'solver.type',      {'cg' 'relax'});   % Solver type (fmg|relax|cg)
% o = setdefault(o, 'solver.nbiter',    [40 10]);          % Number of iterations of the linear solver
% o = setdefault(o, 'solver.tolerance', 1E-3);             % Solver gain threshold 
% o = setdefault(o, 'solver.verbose',   1);                % Solver verbosity 
% o = setdefault(o, 'solver.sumtype',   'double');         % Solver accumulator type (double|native)
% o = setdefault(o, 'solver.precond',   false);            % [CG only] use preconditioner 

% --- Reformat options
o.vs               = sr_padarray(o.vs(:)',               [0 3-numel(o.vs)],                                  'replicate', 'post');
% o.solver.nbiter    = sr_padarray(o.solver.nbiter(:)',    [0 numel(o.solver.type)-numel(o.solver.nbiter)],    'replicate', 'post');
% o.solver.tolerance = sr_padarray(o.solver.tolerance(:)', [0 numel(o.solver.type)-numel(o.solver.tolerance)], 'replicate', 'post');
% o.solver.verbose   = sr_padarray(o.solver.verbose(:)',   [0 numel(o.solver.type)-numel(o.solver.verbose)],   'replicate', 'post');
% o.solver.sumtype   = sr_padarray(o.solver.sumtype(:)',   [0 numel(o.solver.type)-numel(o.solver.sumtype)],   'replicate', 'post');
% o.solver.precond   = sr_padarray(o.solver.precond(:)',   [0 numel(o.solver.type)-numel(o.solver.precond)],   'replicate', 'post');

% --- Get solver handle
% switch lower(o.solver.type)
%     case 'relax'
%         o.solver.fun = @sr_solve_l1_relax;
%     case 'cg'
%         o.solver.fun = @sr_solve_l1_cg;
% end

% --- Ensure output folder exists
path = strsplit(o.out.folder, filesep);
if isempty(path{1}), path{1} = filesep; end
% TODO: ^ not sure about this, especially on windows...
for i=1:numel(path)
    folder = fullfile(path{1:i});
    if ~exist(folder, 'dir')
        st = mkdir(folder);
        if ~st
            error('Cannot create output folder %s', folder);
        end
    end
end

% -------------------------------------------------------------------------
function o = setdefault(o, field, value, exists)
% Set default value in option structure. The value is only set if none 
% existed previously.
%
% FORMAT opt = setdefault(opt, field, value)
% opt   - Structure
% field - Hierarchy if fields (cell, or '.' separated)
% value - Default value
%
% EXAMPLE
% >> opt.my.field1 = 1;
% >> opt = setdefault(opt, 'my.field1', 2);
% >> opt.my.field1
% ans = 1
% >> opt = setdefault(opt, 'my.field2', 3);
% >> opt.my.field2
% ans = 3
if nargin < 4, exists = false; end
if ~iscell(field), field = strsplit(field, '.'); end

if isempty(field)
    if ~exists, o = value; end
else
    exists = isfield(o, field{1});
    if ~exists, o.(field{1}) = []; end
    o.(field{1}) = setdefault(o.(field{1}), field(2:end), value, exists);
end