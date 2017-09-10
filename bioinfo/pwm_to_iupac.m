% Convert a pwm to the closest iupac character
% Input: 
% pwm - input pwms 
% metric - what distance to iupac to use
% 
% Output: 
% iupac_seq - sequence in iupac format 
%
function iupac_seq = pwm_to_iupac(pwm, metric, varargin)

AssignGeneralConstants();
AssignStatsConstants();
iupac_vec = 'ACGTMRWSYKBDHVN';
iupac_pwms = derich_correct(iupac_to_pwm(iupac_vec), epsilon, 4);
num_cols = size(pwm,2); % num_iupac = length(iupac_vec);

if(~exist('metric', 'var') || isempty(metric))
    metric = EUCLIDIAN; 
end
switch metric
    case KULLBACK_LEIBLER
        pwm = derich_correct(pwm, epsilon, 4);
        distance_mat = repmat(sum(iupac_pwms .* log(iupac_pwms)), num_cols, 1)' - ... % sum_i Q_i * log Q_i
            iupac_pwms' * log(pwm); % sum_i Q_i * log P_i
    case EUCLIDIAN
        distance_mat = two_vecs_to_distance( iupac_pwms', pwm');
end
[~, closest_inds] = min(distance_mat); % find closest points
iupac_seq = iupac_vec(closest_inds);
