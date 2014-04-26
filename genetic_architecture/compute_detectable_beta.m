% Check what effect sizes beta can be detected (not finished)
function compute_detectable_beta()

num_snps = 100000; % number of independent snps (out of LD)
n_samples = 10000; % number of individuals
res = 0.0001; f_vec = res:res:0.5;

var_expl_vec = 2 .* f_vec .* (1-f_vec);
z_score_cutoff = norminv(1 - 0.05 ./ num_snps); % when do effects start being detectable



