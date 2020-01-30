function in = sr_in_coregister(in,fwhm)
% Coregister input volumes together.
%
%   The first volume (`in{1}{1}`) is taken as the reference.
%   This function alters the field 'mat' of each volume (but does
%   not alter the headers of the files on disk).
%
% FORMAT in = sr_coregister(in,[fwhm])

subtmpdir = make_dir(tempdir);
try
    % --- REFERENCE (vol == 1)
    refname = make_file(in{1}{1}, subtmpdir, 'ref');
    ref     = spm_vol(refname);
    ref.mat = in{1}.mat;

    % --- MOVING (vol > 1)
    for c=1:numel(in)
        for r=1:numel(in{1})
            if c == 1 && r == 1, continue; end
            movname      = make_file(in{c}{r}, subtmpdir, 'mov');
            in{c}{r}.trf = sr_coreg(ref, {movname, in{c}{r}.mat}, fwhm);
            in{c}{r}.mat = in{c}{r}.trf\in{c}{r}.mat;
            in{c}{r}.trf = in{c}{r}.mat0/in{c}{r}.mat;
            delete(movname);
        end
    end
    delete(refname);
catch
    rmdir(subtmpdir, 's');
end
% =========================================================================
function dirname = make_dir(dir)
if ~exist(dir,'dir'), error('Temporary directory %s does not exist', dir); end
dirname = fullfile(dir, ['spm super-resolution ' datestr(now)]);
ok = mkdir(dirname);
if ~ok, error('Could not create temporary directory'); end
% =========================================================================
function fname = make_file(in, dir, prefix)
fname   = fillfile(dir, [prefix '.nii']);
nii     = nifti;
nii.dat = file_array(fname, in.dim, class2type(in.dat));
nii.mat = in.mat;
create(nii);
nii.dat(:,:,:) = in.dat(:,:,:);