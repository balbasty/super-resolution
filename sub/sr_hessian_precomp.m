function Ha = sr_hessian_precomp(in, odim, omat, opt)

Ha = cellfun(@(x) cell(1,numel(x)), in, 'UniformOutput', false);
if opt.precomp_H
    one = ones(odim, 'single')
    for c=1:numel(in) % Loop over contrasts
        for r=1:numel(in{c}) % Loop over repeats
            if isfield(in{c}{r}, 'slice'), slice = in{c}{r}.slice; 
            else,                          slice = []; end
            Ha{c}{r} = sr_proj('AtA', one, in{c}{r}.dim, in{c}{r}.mat, ...
                                      omat, opt, slice);  
        end
    end
end