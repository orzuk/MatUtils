% Transfer a joint prob. distribution to case-control distribution (Pr(Z=1)= frac. of cases, default 0.5)
%
% Input:
% p_vec - joint probability of one or two genotypes and phenotype
% n_cases - number of cases (optional)
% n_controls - number of controls (optional)
%
% Output:
% q_vec - adjusted porbability vector considering #cases and #controls
%
function q_vec = pop_prob_to_case_control_prob(p_vec, n_cases, n_controls)

q_vec = p_vec;
if(exist('n_cases', 'var')) % compute cases fraction
    prevalence = n_cases ./ (n_cases + n_controls);
else
    prevalence = 0.5; % default: n_cases = n_controls 
end
if(length(prevalence) == 1)
    prevalence = repmat(prevalence, size(p_vec, 1), 1);
else
    prevalence = prevalence(1:size(p_vec,1)); % take first (we don't support different cases&controls unless also different p_vec)
end

switch size(p_vec,2) % assume row vector
    case 4 % single locus. Assume: P(00), P(01), P(10), P(11) with (X_1, Z)
        q_vec(:,[1 3]) = p_vec(:,[1 3]) .* repmat((1-prevalence) ./ sum(p_vec(:,[1 3]),2), 1, 2);
        q_vec(:,[2 4]) = p_vec(:,[2 4]) .* repmat(prevalence ./ sum(p_vec(:,[2 4]),2), 1, 2);
    case 6 % single locus. Genotypes. Assume: P(00), P(01),  P(10), P(11)  P(20), P(21) with (X_1, Z)
        q_vec(:,[1 2 3]) = p_vec(:,[1 2 3]) .* repmat((1-prevalence) ./ sum(p_vec(:,[1 2 3]),2), 1, 3);
        q_vec(:,[4 5 6]) = p_vec(:,[4 5 6]) .* repmat(prevalence ./ sum(p_vec(:,[4 5 6]),2), 1, 3);
    case 8 % epistasis. Assume: P(000), P(001), P(010), P(011), P(100), P(101), P(110), P(111) with (X_1, X_2, Z)
        q_vec(:,1:2:7) = p_vec(:,1:2:7) .* repmat((1-prevalence) ./ sum(p_vec(:,1:2:7),2), 1, 4);
        q_vec(:,2:2:8) = p_vec(:,2:2:8) .* repmat(prevalence ./ sum(p_vec(:,2:2:8),2), 1, 4);
        
    otherwise
        if(ismember(size(p_vec, 1), [4 6 8])) % transpose
            q_vec = pop_prob_to_case_control_prob(p_vec', n_cases, n_controls)';
        end
end
