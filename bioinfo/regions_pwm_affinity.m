% Compute affinity of a DNA sequence to a pwm.
% We check the local affinity at each position and do some weighted sum
% over them where the relative strentgh of the site is determined by the
% 'temperature' random variable beta. Higher beta means that most sites are
% about the same whereas lower beta means that only the very strong sites
% contribute.
% In addition, we also can incorparate the conservation information, such
% that a site 'match' is determined by both the affinity to the pwm and
% it's total conservation.
%
% Input:
% pwm - the motif's pwm
% packed_seqs - the sequences to search for (in a  packed form)
% seqs_lens - a vector of sequences lengths
% do_log - flag saying if to log-transform the pwms
% strand - which strand to look for (should be both)
% beta - the temperature. -1 means threshold some affinity, -2 means take best site
% loglike_mat - a conservation measure (optional).
% alpha - top fraction of conserved to consider (optional. Default is 0.05)
%
% Output:
% regions_affinity - a number representing the total affinity of the regions to the pwm
%
function regions_affinity = regions_pwm_affinity(pwm, packed_seqs, seqs_lens, do_log, strand, ...
    beta, loglike_mat, alpha, varargin)

Assign24MammalsGlobalConstants;

% first compute the bound loc-scores
[bound_loc_scores_mat bound_loc_scores_mat_rev_comp] =  calc_all_sites_scores_matlab(pwm, ...
    packed_seqs, seqs_lens, do_log, [], strand); % assume they're all equal and take the first. Run on both strands

if(nargin >= 7) % here we have the input loglike_mat
    if(~exist('alpha', 'var'))
        alpha = 0.05; % conservation cutoff (fraction). Default is 0.05
    end
    if(isempty(alpha))
        alpha = 0.05; 
    end
    temp_len = 0;
    for i=1:length(bound_loc_scores_mat) % loop over regions
        temp_len = temp_len + length(loglike_mat{i});
    end
    q = quantile(cell2vec(loglike_mat, [], COLUMN), 1-alpha); % compute the cut-off for conservation. One for all regions
    for i=1:length(bound_loc_scores_mat) % loop over regions
        bound_loc_scores_mat{i}(loglike_mat{i} <= q) = -99999999;  % remove all the non-conserved ones and set them to minimum
        bound_loc_scores_mat_rev_comp{i}(loglike_mat{i} <= q) = -99999999;  % remove all the non-conserved ones and set them to minimum
    end
end

for i=1:length(bound_loc_scores_mat) % combine both strands (this somehow affects results greatly !!!)
     bound_loc_scores_mat{i} = [bound_loc_scores_mat{i} bound_loc_scores_mat_rev_comp{i}];
end    


n = length(bound_loc_scores_mat); % we assume it is a cell-array
regions_affinity = zeros(1,n, 'single'); beta = single(beta);  % make singles to save memory
switch beta
    case -1 % in this special flag, we just count the # of binding sites above a certain threshold
        alpha=0.002; cel_pwm{1} = pwm;
        threshold = pwms_to_thresholds(cel_pwm, alpha);
        for i=1:length(bound_loc_scores_mat) % Just count # of hits above a threshold ! (both strands) 
            regions_affinity(i) = sum(bound_loc_scores_mat{i} > threshold);  
        end
    case -2 % in this special case we just take the maximum affinity (log-likelihood ratio score)
        for i=1:length(bound_loc_scores_mat) % Just take the best hit ! (look at both strands!) 
            regions_affinity(i) = max(bound_loc_scores_mat{i});  
        end
    otherwise % here we actually compute the sum of some energies of binding
        for i=1:length(bound_loc_scores_mat)
            regions_affinity(i) = sum(exp(beta .* bound_loc_scores_mat{i}));
        end
end





