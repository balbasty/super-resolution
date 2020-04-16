function out = sr_out_write(out, opt)

C = size(out.dat,4);
% Recon
if strcmpi(opt.out.mem, 'load')
    dats = cell(1,C);
    for c=1:C
        nii      = nifti;
        nii.dat  = file_array(out.fnames{c}, out.dim, 'float32');
        nii.mat  = out.mat;
        nii.mat0 = out.mat;
        create(nii);
        if opt.log
            nii.dat(:,:,:) = exp(out.dat(:,:,:,c));
        else
            nii.dat(:,:,:) = out.dat(:,:,:,c);
        end
        dats{c} = nii.dat;
    end
    out.dat = cat(4, dats{:});
elseif opt.log
    for c=1:C
        out.dat(:,:,:,c) = exp(out.dat(:,:,:,c));
    end
end

% Weights
if opt.reg.mode == 1 && strcmpi(opt.out.mem, 'load')
    nii      = nifti;
    nii.dat  = file_array(out.fnames{C+1}, out.dim, 'float32');
    nii.mat  = out.mat;
    nii.mat0 = out.mat;
    create(nii);
    nii.dat(:,:,:) = out.rls(:,:,:);
    out.rls  = nii.dat;
end