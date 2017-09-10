% Estimate heritability using family segregation analysis.
% For continuous traits, use Falconer's formula:
% h = 2*(r_MZ - r_DZ) -  this should give the broad-sense (?) heritability.
% The narrow sense (?) heritability is given by:
% h_add = r_sib - r_ad, wherew ad is adoptee (i.e. someone with no genetic relation)
% (this is not exact: dominance variance do not go in but epistasis between
% different loci does go into the calculation)
%
% Input:
% family_tree - structure of family
% family_trait_values - value of trait in each member (multiple copies)
% scale - ???
% family_model - how to compute heritability (default: ACE model)
%
% Output:
% h - (broad sense) heritability
% h_add - (narrow sense) additive heritability
% h_liability - (narrow sense) additive heritability on the liability scale
% thresh - threshold for disease in liability model
% mean_liability_given_affected - liability of an affected person
% lambda_mz  - mz twin relative risk
% lambda_dz - dz twin relative risk
% C - shared envirounment component
% E - unique envirounment component
%
function [h h_add h_liability thresh mean_liability_given_affected lambda_mz lambda_dz] = ...
    family_segregation_to_heritability(family_tree, family_trait_values, scale, family_model, varargin)

AssignGeneralConstants;
if(~exist('scale', 'var') || isempty(scale))
    scale = 'observed';
end
if(~exist('family_model', 'var') || isempty(family_model))
    family_model = 'ACE';
end


prevalence = mean(family_trait_values(end,:)); % estimate of preavalence
lambda_mz = mean(family_trait_values(end,:).*family_trait_values(end-1,:)) / prevalence^2; % monozigotic twin relative risk
lambda_dz = mean(family_trait_values(end,:).*family_trait_values(end-2,:)) / prevalence^2; % dizigotic twin relative risk

r_MZ = corr(vec2column(family_trait_values(end,:)), vec2column(family_trait_values(end-1,:))); % compute correlation between monozigotic twins
r_DZ = corr(vec2column(family_trait_values(end,:)), vec2column(family_trait_values(end-2,:))); % compute correlation between dizogotic twins (siblings)


% what is the scale we look at
switch lower(scale)
    case {'observed', 'binary'} % compute heritability on the observed binary scale
        h = mz_twin_risk_to_heritability(lambda_mz, prevalence); % use conversion formula
        h_add = h; % for now .. wrong
        h_add = 2*(r_MZ - r_DZ); % assume we've estimating it as a binary heritability
        h_liability = heritability_scale_change(h_add, 'liability', prevalence); % transfer to liability
    case 'gaussian' % here all (both input&output) is gaussian
        r_MZ = corr(family_trait_values(end-1,:)', family_trait_values(end,:)'); % monozygotic twin values
        r_DZ = corr(family_trait_values(end-2,:)', family_trait_values(end,:)'); % dizygotic twin values
        r_AD = corr(family_trait_values(1,:)', family_trait_values(end,:)'); % unrelated individuals values (adoptees)
        h = 2*(r_MZ - r_DZ); % broad(?) sense heritability (ACE model)
        h_add = r_DZ - r_AD; % narrow sense heritability. Dizygotic are like sibling
        C = r_MZ - h; % = 2 r_DZ – r_MZ; % common envirounment component (ACE model)
        E = 1 - r_MZ; % unique envirounment component (ACE model)
        
    case 'liability' % compute heritability on the hidden liability scale (input is binary)
        thresh = norminv(1-prevalence); % compute liability threshold
        mean_liability_given_affected = exp(-thresh.^2./2) ./ (sqrt(2*pi) * prevalence); % mean value of liability given that you're affected
        
        
        %        r_MZ = familial_risk_to_heritability(lambda_mz, h_scale, mu, k_R)
        r_MZ = familial_risk_to_heritability(lambda_mz, 'liability', prevalence, 1);
        r_DZ = 0.5*familial_risk_to_heritability(lambda_dz, 'liability', prevalence, 0.5);
        switch family_model % How to compute heritability depends on model
            case 'exact_narrow'
                min_rho = 0; max_rho = 1; % find correlation using a stupid but safe bisection method. Assumes positive correlation
                while (max_rho - min_rho > epsilon)
                    mid_rho = (min_rho + max_rho)/2 % bisect
                    sigma_mat = [1 mid_rho ; mid_rho 1];
                    computed_lambda_mz = mvncdf([-thresh -thresh], [0 0], sigma_mat) / prevalence^2; % compute the bivariate integral (take x < -threshold)
                    
                    if(computed_lambda_mz < lambda_mz)
                        min_rho = mid_rho;
                    else
                        max_rho = mid_rho;
                    end
                end
                h = mid_rho; h_add = mid_rho; h_liability = h_add;
                
            case {'ACE', 'ADE'} % here esitmate heritability from a family by taking DZ and MZ risk
                switch family_model
                    case 'ACE'
                        h_liability = 2*(r_MZ - r_DZ); % standard ACE formula
                    case 'ADE'
                        h_liability = (4*r_DZ - r_MZ)/3; % biased ADE formula
                end
                h = h_liability; h_add = h_liability;
        end
        
end






