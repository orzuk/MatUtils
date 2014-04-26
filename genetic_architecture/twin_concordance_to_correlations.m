% Compute the MZ and DZ correlation on the liability scale from counts of twin studies.
% The pairwise concordance is  C/(C+D), where C is the number of concordant
% pairs and D is the number of discordant pairs.
% probandwise_concordance is 2C/(2C+D)
% 
% Input: 
% mu -
% pairwise_concordance -
% probandwise_concordance - 
% 
% Output: 
% r_twin - correlation in trait's value of twins
% 
function r_twin = twin_concordance_to_correlations(mu, ...
    pairwise_concordance,  probandwise_concordance)

% r_MZ = pairwise_concordance; % Temp - wrong!

threshold = norminv(1-mu); z_score = normpdf(threshold); % Get Gaussian height at the incidence threshold
mean_deviate = z_score / mu; % approximate mean liability of affected individuals
twin_threshold = norminv(1-probandwise_concordance);

r_twin = threshold - twin_threshold * ...
    sqrt(1 - (threshold^2 - twin_threshold^2)*(1-threshold / mean_deviate));
r_twin = r_twin / (mean_deviate + twin_threshold^2*(mean_deviate-threshold));


