% Test a model given by a 8x1 probability vec for epistasis interaction
%
% Input:
% p_vec - vectof of 8 probabilities (2 genotypes and a phenotype)
% test_stat - what null non-epistatic model to use (e.g. additive, multiplicative etc.)
% epsilon - maximum deviation from null model allowed
%
% Output:
% epistasis_flag - flag indicating if model shows epistasis
% null_deviation - by how much is model different than null (no epistasis)
%
function [epistasis_flag null_deviation] = test_model_for_epistasis(p_vec, test_stat, epsilon)

p_z_given_x1_x2 = p_vec(2:2:8) ./ (p_vec(1:2:7) + p_vec(2:2:8)); % Pr(Z=1|x1,x2)
p_x1_x2 = p_vec(1:2:7) + p_vec(2:2:8); % Pr(x1,x2). A 2x2 (or 4x1) table


test_stat = strdiff(test_stat, '-analytic'); 
test_stat = strdiff(test_stat, 'analytic'); 
switch test_stat
    case {'additive', 'chi-square'}
        null_deviation = sum(p_z_given_x1_x2([1 4])) - sum(p_z_given_x1_x2([2 3]))
    case 'multiplicative'
        null_deviation = prod(p_z_given_x1_x2([1 4])) - prod(p_z_given_x1_x2([2 3]))
    case {'logistic', 'logit'}
        null_deviation = ...
            prod(p_vec([2 3 5 8]) ./ p_x1_x2) - prod(p_vec([1 4 6 7]) ./ p_x1_x2)
    case 'probit' % What should be the null expectation here?
        null_deviation = 0; % TEMP! WRONG !!!! (we assume no epistasis for now)
    otherwise 
        null_deviation = 0; % TEMP! WRONG !!!! (we assume no epistasis for now) 
end

epistasis_flag = abs(null_deviation) > epsilon; % allow some tolerance
