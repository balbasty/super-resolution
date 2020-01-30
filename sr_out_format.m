function out = sr_out_format(dim, mat, in, opt)

C       = numel(in);
out     = struct;
out.mat = mat;
out.dim = dim;
% -------------------------------------------------------------------------
% Allocate disk/memory
if strcmpi(opt.out.mem, 'map')
    ndig = num2str(ceil(log10(C+0.5)));
    
    % Dir/base
    dir  = cell(1,C);
    base = cell(1,C);
    for c=1:C
        dir{c} = opt.out.dir;
        if isa(in{c}{1}.dat, 'file_array')
            fname = {in{1}{1}.dat.fname};
            fname = fname{1};
            if isempty(opt.out.dir)
                dir{c} = fileparts(fname);
            end
            [~,base{c}] = fileparts(fname);
            base{c}     = ['_' base{c}];
        else
            if isempty(opt.out.dir)
                dir{c}  = '.';
            end
            base{c} = '';
        end
    end
    
    % Recon
    dats = cell(1,numel(in));
    for c=1:C
        nii     = nifti;
        nii.dat = file_array(fullfile(dir{c}, sprintf(['sr_' '%0' ndig 'd' base{c} '.nii'], c)), dim, 'float32');
        nii.mat = mat;
        create(nii);
        dats{c} = nii.dat;
    end
    out.dat = cat(4, dats{:});
    
    % Weights
    if opt.reg.mode == 1
        nii     = nifti;
        nii.dat = file_array(fullfile(dir{1}, 'rls.nii'), dim, 'float32');
        nii.mat = mat;
        create(nii);
        out.rls = nii.dat;
    end
    
    % Uncertainty
    if opt.reg.uncertainty == 0
        dats = cell(1,numel(in));
        for c=1:C
            nii     = nifti;
            nii.dat = file_array(fullfile(dir{c}, sprintf(['u_sr_' '%0' ndig 'd' base{c} '.nii'], c)), dim, 'float32');
            nii.mat = mat;
            create(nii);
            dats{c} = nii.dat;
        end
        out.uncty = cat(4, dats{:});
    end
    
else
    out.dat = zeros([dim C], 'single');
    if opt.reg.mode == 1, out.rls = zeros(dim, 'single'); end
    if opt.reg.uncertainty == 0, out.uncty = zeros([dim C], 'single') + 1E-3; end
end

% -------------------------------------------------------------------------
% Initialise
out.lam = sr_padarray(opt.reg.value(:)', [0 C-numel(opt.reg.value)], 'replicate', 'post');
for c=1:C
    % Recon
    meanmu = 0;
    for r=1:numel(in{c}), meanmu = meanmu + in{c}{r}.mu; end
    meanmu = meanmu/numel(in{c});
    if opt.log, meanmu = log(meanmu); end
    out.dat(:,:,:,c) = meanmu;
    if ~opt.log, out.lam(c) = out.lam(c)/abs(meanmu).^2; end
    
    % Uncertainty
    if opt.reg.uncertainty == 0
        out.uncty(:,:,:,c) = 1E-3;
    end
end
% Weights
if opt.reg.mode == 1
    out.rls(:) = 1;
end