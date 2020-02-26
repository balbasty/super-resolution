function f = sr_plot_progress(in,out,ll,opt,figname)

% ---------------------------------------------------------------------
% Get figure object
if nargin < 5 || isempty(figname)
    figname = '[sr] fit';
end
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);


ncol = size(out.dat,4) + 1;
if isfield(out, 'uncty')
    nrow = 3;
else
    nrow = 2;
end
i = 1;

z = ceil(size(out.dat,3)/2);

% -------------------------------------------------------------------------
% Observed
for c=1:numel(in)
    subplot(nrow,ncol,i);
    dmu = [size(out.dat) 1];
    dmu = dmu(1:3);
    y   = warps_affine(dmu, in{c}{1}.mat\out.mat);
    y   = y(:,:,z,:);
    img = spm_diffeo('bsplinc', single(in{c}{1}.dat()), [0 0 0  0 0 0]);
    img = spm_diffeo('bsplins', img, y, [0 0 0  0 0 0]);
    % img = in{c}{1}.dat(:,:,z);
    imagesc(img);
    colormap('gray');
    colorbar;
    axis off
    title('Observed')
    i = i+1;
end

% -------------------------------------------------------------------------
% Free energy
if ~isempty(ll)
    ax = subplot(nrow,ncol,i);
    cla(ax);
    hold on
    colors = [0 0 0; hsv(size(ll,1)-1)];
    lines  = {2 1 1 1 1};
    titles = {'Total Energy' 'Data' 'Bound' 'TV'};
    for k=1:size(ll,1)
        p = plot(ll(k,:));
        p.Color = colors(k,:);
        p.Marker = 's';
        p.MarkerEdgeColor = 'none';
        p.LineWidth = lines{k};
    end
    hold off
    title('Free Energy');
    legend(titles{:});
end
i = i + 1;

% -------------------------------------------------------------------------
% Recon
for c=1:size(out.dat,4)
    subplot(nrow,ncol,i);
    img = out.dat(:,:,z,c);
    if opt.log, img = exp(img); end
    imagesc(img);
    colormap('gray');
    colorbar;
    axis off
    title('Recon')
    i = i+1;
end

% -------------------------------------------------------------------------
% Weights
if isfield(out, 'rls')
    subplot(nrow,ncol,i);
    img = out.rls(:,:,z);
    imagesc(img);
    colormap('gray');
    colorbar;
    axis off
    title('Weights')
end
i = i+1;

% -------------------------------------------------------------------------
% Uncertainty
if isfield(out, 'uncty')
    for c=1:size(out.dat,4)
        subplot(nrow,ncol,i);
        img = out.uncty(:,:,z,c);
        imagesc(img);
        colormap('gray');
        colorbar;
        axis off
        title('Uncertainty')
        i = i+1;
    end
end

drawnow

% === HELPERS =============================================================
function psi = warps_affine(lat, mat)
id  = warps_identity(lat);
psi = reshape(reshape(id,[prod(lat) 3])*mat(1:3,1:3)' + mat(1:3,4)',[lat 3]);
if lat(3) == 1, psi(:,:,:,3) = 1; end
% -------------------------------------------------------------------------
function id = warps_identity(d)
id = zeros([d(:)' 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));