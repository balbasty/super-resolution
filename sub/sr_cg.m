function [x,info] = sr_cg(A,b,x,iM,nit,tol,verbose,sumtype)
% Conjugate gradient solver.
% CG can be used to solve (large) linear systems: A*x = b
%
% FORMAT [x,info] = sr_cg(A,b,[x0],[iM],[nit],[tol],[verbose],[sumtype])
% A       - Function handle for the linear operator A*x (left-hand side)
% b       - Target array (right-hand side)
% x0      - Initial guess [0]
% iM      - Function handle for the (inverse) preconditioning matrix [id]
% nit     - Maximum number of iterations [32]
% tol     - Tolerance for early stopping [1E-3]
% verbose - Verbosity level [0]
% sumtype - Accumulator type 'native'/['double']
% x       - A\b
% info    - Structure with fields:
%           . nbiter - Effective number of iterations
%           . rr     - Root mean squared residuals (per iteration)
%           . time   - Time (per iteration)
%
% Note that summing in single is twice as fast as summing in double, but
% the additional precision is often needed.
%
% Reminder: CG ensures a monotonic decrease of the *A-norm* of the errors
% ||Ax-b||_A = (Ax-b)'*A\(Ax-b)

% Adapted from Mikael's cg_im_solver.

DEBUG = false;

if nargin < 3 || isempty(x)
    x       = zeros(size(b),'single');
end
if nargin < 4 || isempty(iM)
    iM       = @(y) y;
end
if nargin < 5 || isempty(nit) || isnan(nit)
    nit     = 32;
end
if nargin < 6 || isempty(tol) || isnan(tol)
    tol     = 1e-3;
end
if nargin < 7 || isempty(verbose) || isnan(verbose)
    verbose = true;
end
if nargin < 8 || isempty(sumtype) || isnan(sumtype)
    sumtype = 'double';
end

sr_plot_interactive('Init', 'Conjugate Gradient', 'Residuals', 'Iteration');

% -------------------------------------------------------------------------
% Initialisation  

bb = sqrt(sum(b(:).^2, sumtype));               % Norm of b: sqrt(b'*b)
r  = b - A(x);                                  % Residual: b - A*x
z  = iM(r);                                     % Preconditioned residual

rr = sqrt(sum(r(:).^2, sumtype))/bb;            % Norm of r: sqrt(r'*r)
rz = sum(r(:).*z(:), sumtype);                  % Inner product of r and z
p     = z;                                      % Initial conjugate directions p
beta  = 0;                                      % Initial step size
  
if DEBUG
    ee   = A(x)-2*b;                            % A-norm of the error
    ee   = sum(ee(:).*x(:), sumtype);           % (up to a constant)
    
    sr_plot_interactive('Set', 0, rr);
    start = tic;
    time  = [];
end

if verbose, fprintf('Conjugate gradient: '); end

% -------------------------------------------------------------------------
% Run algorithm
for j=1:nit
    % ---------------------------------------------------------------------
    % Calculate conjugate directions P which defines the direction of descent
    p = z + beta*p;
    clear z

    % ---------------------------------------------------------------------
    % Finds the step size of the conj. gradient descent
    Ap    = A(p);
    alpha = rz / sum(p(:).*Ap(:), sumtype);
    
    % ---------------------------------------------------------------------
    % Perform conj. gradient descent, obtaining updated X and R, using the 
    % calculated P and alpha
    x = x + alpha * p; 
    r = r - alpha * Ap;
    clear Ap
    
    % ---------------------------------------------------------------------
    % Check convergence
    rr0  = rr(end);
    rr1  = sqrt(sum(r(:).^2, sumtype))/bb;
    rr   = [rr rr1];
    if DEBUG
        time = [time toc(start)];
        e1   = A(x)-2*b;
        e1   = sum(e1(:).*x(:), sumtype);
        ee   = [ee e1];
        sr_plot_interactive('Set', j, rr1);
    end
    if rr1 > rr0
        if verbose, fprintf('x'); end
    elseif rr1 < tol
        fprintf('o');
        break;
    elseif verbose
        fprintf('.');
    end
    
    % ---------------------------------------------------------------------
    % Update preconditioned residual  
    z = iM(r);
    
    % ---------------------------------------------------------------------
    % Finds the step size for updating P
    rz0  = rz;
    rz   = sum(r(:).*z(:), sumtype);
    beta = rz / rz0;
end
if verbose, fprintf('\n'); end

sr_plot_interactive('Clear');
info.nbiter = j;
info.rr     = rr;
if DEBUG
    info.time   = time;
    info.ee     = ee;
end