function scl = sr_diag_kernel(vs)

scl = spm_diffeo('kernel', [3 3 3], [vs 0 1 0 0 0]);
scl = double(abs(scl(1,1,1)));