function root = tbx_cfg_superresolution

if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','super-resolution'));
    addpath(fullfile(spm('Dir'),'toolbox','Longitudinal'));
    addpath(fullfile(spm('Dir'),'toolbox','Shoot'));
end


% -------------------------------------------------------------------------
plane         = cfg_menu;
plane.tag     = 'plane';
plane.name    = 'In-plane profile';
plane.help    = {
    'Profile of the forward kernel in the in-plane directions.'
    }';
plane.labels = {
                'Rectangle'
                'Gaussian'
                }';
plane.values = {
                0
                2
                }';
plane.val    = {2};
% -------------------------------------------------------------------------
slice         = cfg_menu;
slice.tag     = 'slice';
slice.name    = 'Slice profile';
slice.help    = {
    'Profile of the forward kernel in the slice (= thickest) direction.'
    }';
slice.labels = {
                'Rectangle'
                'Gaussian'
                }';
slice.values = {
                0
                2
                }';
slice.val    = {2};
% ----------------------------------------------v---------------------------
gap         = cfg_entry;
gap.tag     = 'gap';
gap.name    = 'Slice gap';
gap.help    = {
    'Gap between slices (in percent).'
    ['In clinical 2D acquisitions, there is usually a gap between ' ...
     'slices (i.e., only the central part of the slice id excited).']
    }';
gap.strtype = 'r';
gap.num     = [1 1];
gap.val     = {1/3};
% -------------------------------------------------------------------------
acc         = cfg_menu;
acc.tag     = 'acc';
acc.name    = 'Integration mode';
acc.help    = {
    'Integration mode of the forward kernel:'
    '* Average: intensity is assumed to be independent of the voxel size,'
    '* Accumulate: intensity is assumed to scale linearly with voxel size.'
    }';
acc.labels = {
                'Average'
                'Accumulate'
                }';
acc.values = {
                0
                1
                }';
acc.val    = {0};
% -------------------------------------------------------------------------
function elem = profile
elem         = cfg_branch;
elem.tag     = 'profile';
elem.name    = 'Forward profile';
elem.val     = {plane slice gap acc};
end
% -------------------------------------------------------------------------
function elem = vol(mode)
switch mode
case 'd'
elem         = cfg_files;
elem.tag     = 'files';
elem.name    = 'Files';
elem.help    = {
    'Select scans from this channel for processing.'
    }';
elem.filter  = 'image';
elem.ufilter = '.*';
elem.num     = [1 Inf];
elem.preview = @(f) spm_check_registration(char(f));
case 's'
elem        = cfg_branch;
elem.tag    = 'vol';
elem.name   = 'Volume';
elem.val    = {vol('d') profile};
end
end
% -------------------------------------------------------------------------
function elem = channel(mode)
elem                = cfg_branch;
elem.tag            = 'channel';
elem.name           = 'Channel';
switch mode
case 'd'
elem.val            = {vol(mode)};
case 's'
elem.val            = {cfg_repeat};
elem.val{1}.tag     = 'vols';
elem.val{1}.name    = 'Volumes';
elem.val{1}.help    = {'Select individual input volumes.'};
elem.val{1}.values  = {vol(mode)};
elem.val{1}.val     = {vol(mode)};
elem.val{1}.num     = [1 Inf];
end
end
% -------------------------------------------------------------------------
function elem = data(mode)
switch mode
case 'd'
    verb = 'denoised';
case 's'
    verb = 'super-resolved';
end
elem         = cfg_repeat;
elem.tag     = 'channels';
elem.name    = 'Channels';
elem.values  = {channel(mode)};
elem.help    = {
    'Specify the number of different channels.'
    ['If you have scans of different contrasts that represent the ' ...
     'same subject/organ/object, they can be ' verb ' together. ' ...
     'The model assumes that edges are localised in similar locations ' ...
     'across contrasts.']
    }';
elem.val     = {channel(mode)};
elem.num     = [1 Inf];
end
% -------------------------------------------------------------------------
function elem = modelog
elem        = cfg_menu;
elem.tag    = 'log';
elem.name   = 'Log-encoding';
elem.help   = {
    'Encode the reconstructed images by their log.'
    ['This ensures their positivity, and makes the regularisation ' ...
     'scale-independant. However, the algorithm becomes a bit slower.']
    }';
elem.labels = {'No' 'Yes'};
elem.values = {0 1};
elem.val    = {0};
end
% -------------------------------------------------------------------------
function elem = solver(val)
elem                = cfg_entry;
switch val
case 'fmg'
elem.tag            = 'fmg';
elem.name           = 'FMG';
elem.help           = {
    'Parameters of the full-multi-grid solver:'
    '(1) Number of Full Multigrid cycles,'
    '(2) Number of relaxation iterations per cycle.'
    };
elem.strtype        = 'i';
elem.val            = {[2 2]};
elem.num            = [1 2];
case 'cg'
elem.tag            = 'cg';
elem.name           = 'CG';
elem.help           = {
    'Number of conjugate gradient iterations.'
    };
elem.strtype        = 'i';
elem.val            = {50};
elem.num            = [1 1];
case 'relax'
elem.tag            = 'relax';
elem.name           = 'Relax';
elem.help           = {
    'Number of relaxation iterations.'
    };
elem.strtype        = 'i';
elem.val            = {0};
elem.num            = [1 1];
end
end
% -------------------------------------------------------------------------
function elem = solvers(val)
elem                = cfg_branch;
elem.tag            = 'solver';
elem.name           = 'Solver';
switch val
case 1
elem.val            = {solver('fmg') solver('cg') solver('relax')};
case 2
elem.val            = {solver('fmg')};
end
end
% -------------------------------------------------------------------------
function elem = regtype(val)
switch val
case 0
elem                = cfg_const;
elem.tag            = 'none';
elem.name           = 'None';
elem.val            = {0};
elem.hidden         = true;
case 1
elem                = cfg_branch;
elem.tag            = 'tv';
elem.name           = 'Joint Total-Variation';
elem.val            = {cfg_entry cfg_entry solvers(1)};
elem.val{1}.tag     = 'val';
elem.val{1}.name    = 'Value';
elem.val{1}.strtype = 'r';
elem.val{1}.val     = {5E4};
elem.val{1}.num     = [1 1];
elem.val{1}.help    = {
    'Regularisation value (expected variance of the gradients).'
    }';
elem.val{2}.tag     = 'smooth';
elem.val{2}.name    = 'Smoother';
elem.val{2}.strtype = 'r';
elem.val{2}.val     = {1E-3};
elem.val{2}.num     = [1 1];
elem.val{2}.help    = {
    'RLS smoothing term (should be small).'
    }';
case 2
elem                = cfg_branch;
elem.tag            = 'tkh';
elem.name           = 'Tikhonov';
elem.val            = {cfg_entry solvers(2)};
elem.val{1}.tag     = 'val';
elem.val{1}.name    = 'Value';
elem.val{1}.strtype = 'r';
elem.val{1}.val     = {5E4};
elem.val{1}.num     = [1 1];
elem.val{1}.help    = {
    'Regularisation value (expected variance of the gradients).'
    }';
end
end
% -------------------------------------------------------------------------
function elem = reg
elem        = cfg_choice;
elem.tag    = 'reg';
elem.name   = 'Regularisation type';
elem.help   = {
    'Type of image prior:'
    '* None: maximum-likelihood reconstruction,'
    '* Joint Total-Variation: edge-preserving prior (\ell_{2,1})'
    '* Tikhonov: smoothing prior (\ell_2)'
    }';
elem.values = {regtype(0) regtype(1) regtype(2)};
elem.val    = elem.values(2);
end
% -------------------------------------------------------------------------
function elem = coregtype(val)
switch val
case 0
elem                = cfg_const;
elem.tag            = 'no';
elem.name           = 'No';
elem.val            = {0};
elem.hidden         = true;
case 1
elem                = cfg_branch;
elem.tag            = 'yes';
elem.name           = 'Yes';
elem.val            = {cfg_entry};
elem.val{1}.tag     = 'fwhm';
elem.val{1}.name    = 'FWHM';
elem.val{1}.strtype = 'r';
elem.val{1}.val     = {[21 14 7]};
elem.val{1}.num     = [1 Inf];
elem.val{1}.help    = {
    'List of full-width at half-maximum to use for co-registration.'
    'The FWHM specifies the extent of smoothing of the joint-histogram.'
    'Starting with more smoothing can help with difficult cases.'
    }';
end
end
% -------------------------------------------------------------------------
function elem = coreg
elem        = cfg_choice;
elem.tag    = 'coreg';
elem.name   = 'Co-registration';
elem.help   = {
    'Co-register input volumes.'
    }';
elem.values = {coregtype(0) coregtype(1)};
elem.val    = elem.values(2);
end
% -------------------------------------------------------------------------
function elem = itermax
elem            = cfg_entry;
elem.tag        = 'max';
elem.name       = 'Max';
elem.help       = {
    'Maximum number of Gauss-Newton / Reweighted Least Squares iterations.'
    };
elem.strtype    = 'i';
elem.val        = {10};
elem.num        = [1 1];
end
% -------------------------------------------------------------------------
function elem = tol
elem            = cfg_entry;
elem.tag        = 'tol';
elem.name       = 'Tolerance';
elem.help       = {
    'Tolerance for early stopping.'
    'The tolerance relates to the gain in joint log-likelihood.'
    };
elem.strtype    = 'r';
elem.val        = {1E-3};
elem.num        = [1 1];
end
% -------------------------------------------------------------------------
function elem = armijo
elem            = cfg_entry;
elem.tag        = 'armijo';
elem.name       = 'Armijo factors';
elem.help       = {
    'Gauss-Newton damping factors.'
    ['If the damping factor is larger than 1, the Gauss-Newton update ' ...
    'step is reduced by as much. This prevents the algorithm to '...
    'overshoot during the first few iterations.']
    };
elem.strtype    = 'r';
elem.val        = {[2 1]};
elem.num        = [1 Inf];
end
% -------------------------------------------------------------------------
function elem = iter
elem        = cfg_branch;
elem.tag    = 'iter';
elem.name   = 'Iterations';
elem.val    = {itermax tol armijo};
end
% -------------------------------------------------------------------------
function elem = verbose
elem        = cfg_menu;
elem.tag    = 'verbose';
elem.name   = 'Verbosity';
elem.labels = {
    'Quiet'
    'Print'
    'Plot'
    };
elem. values = {
    0
    1
    2
    };
elem.val = {1};
end
% -------------------------------------------------------------------------
function elem = mapout
elem         = cfg_menu;
elem.tag     = 'map';
elem.name    = 'Memory Map';
elem.help    = {
    'Memory-map output data.'
    'Memory-mapping requires more i/o operations but saves RAM.'
    };
elem.labels  = {'No' 'Yes'};
elem.values  = {0 1};
elem.val     = {1};
end
% -------------------------------------------------------------------------
function elem = prefix
elem            = cfg_entry;
elem.tag        = 'prefix';
elem.name       = 'Prefix';
elem.help       = {'Prefix for output filenames'};
elem.strtype    = 's';
elem.num        = [1 Inf];
elem.val        = {'sr_'};
end
% -------------------------------------------------------------------------
function elem = output
elemi         = cfg_const;
elemi.tag     = 'same';
elemi.name    = 'Same as input';
elemi.help    = {'Output images are written in the same folder as the input files.'};
elemi.val     = {''};
elemi.hidden  = false;
elemo         = cfg_files;
elemo.tag     = 'folder';
elemo.name    = 'Folder';
elemo.help    = {'Select output folder'};
elemo.filter  = 'dir';
elemo.ufilter = '.*';
elemo.num     = [1 1];
elemo.hidden  = false;
elem          = cfg_choice;
elem.tag      = 'output';
elem.name     = 'Output';
elem.values   = {elemi elemo};
elem.val      = {elemi};
end
% -------------------------------------------------------------------------
function elem = io
elem        = cfg_branch;
elem.tag    = 'io';
elem.name   = 'IO';
elem.help   = {'Input/Output options.'};
elem.val    = {output prefix mapout};
end
% -------------------------------------------------------------------------
function elem = vs
elem        = cfg_entry;
elem.tag    = 'vs';
elem.name   = 'Voxel Size';
elem.help   = {
    'Target voxel size.'
    };
elem.strtype = 'r';
elem.val     = {NaN};
elem.num     = [1 1];
end
% -------------------------------------------------------------------------
function elem = opt(mode)
elem        = cfg_branch;
elem.tag    = 'opt';
elem.name   = 'Options';
elem.val    = {};
if mode == 's'
    elem.val{1} = vs;
end
elem.val    = [elem.val {io modelog reg coreg iter verbose}];
end
% -------------------------------------------------------------------------
function elem = sr(mode)
elem         = cfg_branch;
switch mode
case 'd'
    elem.tag  = 'denoise';
    elem.name = 'Denoising';
    elem.help = {'Denoise a multi-channel datasets.'};
case 's'
    elem.tag  = 'superres';
    elem.name = 'Super-Resolution';
    elem.help = {'Super-resolve a multi-channel datasets.'};
end
elem.val     = {data(mode) opt(mode)};
end
% -------------------------------------------------------------------------
function elem = mode
elem         = cfg_choice;
elem.tag     = 'mode';
elem.name    = 'Mode';
elem.help    = {'Choose between denoising and super-Resolution.'};
elem.values  = {sr('d') sr('s')};
elem.val     = elem.values(2);
end
% -------------------------------------------------------------------------
root         = cfg_exbranch;
root.tag     = 'sr';
root.name    = 'Denoising / Super-Resolution';
root.val     = {mode};
root.help    = {'Denoise and/or super-resolve multi-channel datasets.'};
root.prog    = @sr_run;
root.vout    = @vout_create;


end

%==========================================================================
function out = sr_run(job)

% ----
% Mode
% ----

opt      = struct;
opt.mode = fieldnames(job.mode);
opt.mode = opt.mode{1};

Nc = numel(job.mode.(opt.mode).channel);
in = cell(1,Nc);

% -----------------------
% Input and slice profile
% -----------------------

if opt.mode(1) == 's'
% Super-resolution
opt.slice.thickest   = cell(1,Nc);
opt.slice.other      = cell(1,Nc);
opt.slice.gap        = cell(1,Nc);
opt.slice.accumulate = cell(1,Nc);
for c=1:Nc
    in{c} = {};
    for v=1:numel(job.mode.(opt.mode).channel(c).vol)
        files       = job.mode.(opt.mode).channel(c).vol(v).files;
        thickest    = job.mode.(opt.mode).channel(c).vol(v).profile.slice;
        other       = job.mode.(opt.mode).channel(c).vol(v).profile.plane;
        gap         = job.mode.(opt.mode).channel(c).vol(v).profile.gap;
        accumulate  = job.mode.(opt.mode).channel(c).vol(v).profile.acc;
        Nf          = size(files,1);
        in{c}       = cat(1, in{c}, cellstr(job.mode.(opt.mode).channel(c).vol.files));
        opt.slice.thickest{c}   = cat(1, opt.slice.thickest{c},   thickest*ones(Nf,1));
        opt.slice.other{c}      = cat(1, opt.slice.other{c},      other*ones(Nf,1));
        opt.slice.gap{c}        = cat(1, opt.slice.gap{c},        gap*ones(Nf,1));
        opt.slice.accumulate{c} = cat(1, opt.slice.accumulate{c}, accumulate*ones(Nf,1));
    end
end
opt.vs = job.mode.(opt.mode).opt.vs;
else
% Denoising
for c=1:Nc
    in{c} = cellstr(job.mode.(opt.mode).channel(c).files);
end
end

% -------
% Options
% -------

if isfield(job.mode.(opt.mode).opt.io.output, 'same')
    opt.out.folder = '';
else
    opt.out.folder = job.mode.(opt.mode).opt.io.output.folder{1};
end
opt.out.prefix = job.mode.(opt.mode).opt.io.prefix;
if job.mode.(opt.mode).opt.io.map
    opt.out.mem = 'map';
else
    opt.out.mem = 'load';
end
opt.log = job.mode.(opt.mode).opt.log;
regmode = fieldnames(job.mode.(opt.mode).opt.reg);
regmode = regmode{1};
switch regmode
    case 'none'
        opt.reg.mode = 0;
    case 'tv'
        opt.reg.mode = 1;
    case 'tkh'
        opt.reg.mode = 2;
end
opt.reg.value    = job.mode.(opt.mode).opt.reg.(regmode).val;
opt.reg.smo      = job.mode.(opt.mode).opt.reg.(regmode).smooth;
opt.solver.fmg   = job.mode.(opt.mode).opt.reg.(regmode).solver.fmg;
opt.solver.cg    = job.mode.(opt.mode).opt.reg.(regmode).solver.cg;
opt.solver.relax = job.mode.(opt.mode).opt.reg.(regmode).solver.relax;
opt.coreg.do     = isfield(job.mode.(opt.mode).opt.coreg, 'yes');
if opt.coreg.do
    opt.coreg.fwhm = job.mode.(opt.mode).opt.coreg.yes.fwhm;
end
opt.itermax   = job.mode.(opt.mode).opt.iter.max;
opt.tolerance = job.mode.(opt.mode).opt.iter.tol;
opt.armijo    = job.mode.(opt.mode).opt.iter.armijo;
opt.verbose   = job.mode.(opt.mode).opt.verbose;

% ---
% Run
% ---

out = sr_fit(in, opt);
out = out.fnames;

end

%==========================================================================
function cdep = vout_create(job)

mode = fieldnames(job.mode);
mode = mode{1};
for c=1:numel(job.mode.(mode).channel)
    cdep(c) = cfg_dep;
    if mode(1) == 'd'
        cdep(c).sname = sprintf('Denoised contrast (%d)', c);
    else
        cdep(c).sname = sprintf('Super-resolved contrast (%d)', c);
    end
    cdep(c).src_output = substruct('.','channel','()',{c},'.','recon','()',{':'});
    cdep(c).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
end