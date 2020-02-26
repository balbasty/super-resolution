function sr_out_write(out, in, opt)

C = size(out.dat,4);
if strcmpi(opt.out.mem, 'load')
    ndig = num2str(ceil(log10(C+0.5)));
    
    for c=1:C
        dir = opt.out.folder;
        if isa(in{c}{1}.dat, 'file_array')
            fname = {in{1}{1}.dat.fname};
            fname = fname{1};
            if isempty(opt.out.folder)
                dir = fileparts(fname);
            end
            [~,base] = fileparts(fname);
            base     = ['_' base];
        else
            if isempty(opt.out.folder)
                dir  = '.';
            end
            base = '';
        end
        if c == 1, dir1 = dir; end
        % Recon
        nii      = nifti;
        nii.dat  = file_array(fullfile(dir, sprintf(['sr_' '%0' ndig 'd' base '.nii'], c)), out.dim, 'float32');
        nii.mat  = out.mat;
        nii.mat0 = out.mat;
        create(nii);
        nii.dat(:,:,:) = out.dat(:,:,:,c);
    end
    
    % Weights
    if opt.reg.mode == 1
        nii      = nifti;
        nii.dat  = file_array(fullfile(dir1, 'rls.nii'), out.dim, 'float32');
        nii.mat  = out.mat;
        nii.mat0 = out.mat;
        create(nii);
        out.rls  = nii.dat;
    end
end