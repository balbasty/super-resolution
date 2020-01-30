function [nspm,nmatlab] = sr_threads(n, which)
% Set/Get the number of OMP threads.
%
% FORMAT [nspm,nmatlab] = sr_threads(n,which)
% Set the number of threads used by Matlab and/or SPM to n.
% - n can be:
%   . positive number:  set the number of threads to that number
%   . -1:               use the machine's max number of processors (omp_num_procs())
%   . 'automatic':      use Matlab's automatic setting
%   . 'matlab':         use Matlab's current setting
%   . 'spm':            use SPM's current setting
% - which can be:
%   . 'spm':            only sets the setting for SPM
%   . 'matlab':         only sets the setting for Matlab
%   . 'both':           sets the setting for both [default]
%
% FORMAT [nspm,nmatlab] = smx_threads
% Return the actual number of threads used.

if nargin < 2, which = 'both'; end

if nargin > 0
    if ischar(n)
        switch lower(n)
            case 'automatic'
                n0 = maxNumCompThreads('automatic');
                n  = maxNumCompThreads(n0);
            case 'spm'
                n = sscanf(getenv('SPM_NUM_THREADS'),'%d');
            case 'matlab'
                n = maxNumCompThreads;
            otherwise
                error('Unknown option %s', n);
        end
    end
    if any(strcmpi(which, {'spm', 'both'}))
        setenv('SPM_NUM_THREADS',sprintf('%d',n));
        try, spm_diffeo; end
    end
    if any(strcmpi(which, {'matlab', 'both'}))
        maxNumCompThreads(n);
    end
end
nspm    = sscanf(getenv('SPM_NUM_THREADS'),'%d');
nmatlab = maxNumCompThreads;
