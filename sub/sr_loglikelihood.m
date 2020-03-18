function [llx,lly] = sr_loglikelihood(y, out, in, opt)

% -------------------------------------------------------------------------
% Save old value and set new value
% -------------------------------------------------------------------------
y0 = single(out.dat());
out.dat(:,:,:,:) = y;
clear y

vs  = sqrt(sum(out.mat(1:3,1:3).^2));
vol = prod(vs);
K   = size(out.dat,4);

% -------------------------------------------------------------------------
% Data term
% -------------------------------------------------------------------------
if opt.verbose > 0, fprintf('Log-likelihood: data '); end
llx = 0;
for c=1:numel(in)
    llx = llx + sr_gradient(c, in{c}, out, opt);
end
if opt.verbose > 0, fprintf('\n'); end

% -------------------------------------------------------------------------
% Membrane
% -------------------------------------------------------------------------
if opt.reg.mode > 0
    if opt.verbose > 0, fprintf('Log-likelihood: membrane '); end
    switch opt.reg.mode
        % -----------------------------------------------------------------
        % L1 regularisation
        case 1
            w    = single(out.rls());
            lly  = 0.5 * vol * sum(w(:), 'double');
            w    = 1./w;
            for k=1:K
                if opt.verbose > 0, fprintf('.'); end
                y          = single(out.dat(:,:,:,k));
                Ly         = sr_vel2mom_l1(y, out.lam(k), vs, w);
                lly        = lly + 0.5 * vol * sum(y(:).*Ly(:), 'double');
                clear y Ly
            end

        % -----------------------------------------------------------------
        % L2 regularisation
        case 2
            lly = 0;
            for k=1:K
                if opt.verbose > 0, fprintf('.'); end
                y          = single(out.dat(:,:,:,k));
                Ly         = sr_vel2mom_l2(y, out.lam(k), vs);
                lly        = lly + 0.5 * vol * sum(y(:).*Ly(:), 'double');
                clear y Ly
            end
    end
    if opt.verbose > 0, fprintf('\n'); end
end

out.dat(:,:,:,:) = y0;