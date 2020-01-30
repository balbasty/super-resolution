function [x,nit,rr] = sr_cg(A,b,x,iM,nit,tol,verbose,sumtype)
% Conjugate gradient solver.
% CG can be used to solve (large) linear systems: A*x = b
%
% FORMAT [x,nit,nrm,tol] = sr_cg(A,b,[x0],[iM],[nit],[tol],[verbose])
% A       - Function handle for the linear operator A*x (left-hand side)
% b       - Target array (right-hand side)
% x0      - Initial guess [0]
% iM      - Function handle for the (inverse) preconditioning matrix [id]
% nit     - Maximum number of iterations [32]
% tol     - Tolerance for early stopping [1E-3]
% verbose - Verbosity level [0]
% sumtype - Accumulator type 'native'/['double']
%
% Note that summing in single is twice as fast as summing in double, but
% the additional precision is often needed.


% Adapted from Mikael's cg_im_solver.

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

% -------------------------------------------------------------------------
% Initialisation  
bb = sqrt(sum(b(:).^2, sumtype));               % Norm of b: sqrt(b'*b)
r  = b - A(x); clear b                          % Residual: b - A*x
z  = iM(r);                                     % Preconditioned residual

rr = sqrt(sum(r(:).^2, sumtype))/bb;            % Norm of r: sqrt(r'*r)
rz = sum(r(:).*z(:), sumtype);                  % Inner product of r and z
p     = z;                                      % Initial conjugate directions p
beta  = 0;                                      % Initial step size
  
if verbose, fprintf('Conjugate gradient: '); end

% -------------------------------------------------------------------------
% Run algorithm
for j=1:nit
    % ---------------------------------------------------------------------
    % Calculate conjugate directions P which defines the direction of descent
    p = z + beta*p;
    % clear z

    % ---------------------------------------------------------------------
    % Finds the step size of the conj. gradient descent
    Ap    = A(p);
    alpha = rz / sum(p(:).*Ap(:), sumtype);
    
    % ---------------------------------------------------------------------
    % Perform conj. gradient descent, obtaining updated X and R, using the 
    % calculated P and alpha
    x = x + alpha * p; 
    r = r - alpha * Ap;
    % clear Ap
    
    % ---------------------------------------------------------------------
    % Check convergence
    rr0 = rr;
    rr  = sqrt(sum(r(:).^2, sumtype))/bb;    
    if rr > rr0
        if verbose, fprintf('x'); end
        rr = rr0;
        x  = x - alpha * p;
        break
    elseif rr < tol
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

nit   = j;