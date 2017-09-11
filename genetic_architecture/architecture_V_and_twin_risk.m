% Compute variance and pentrance from other architecture parameters 
% 
% Input: 
% mu - trait mean (prevalence)
% v_environment - variance due to environmental factors
% 
% Output: 
% V - trait variance 
% mz_twin_risk - probability of twin diseased given identical twin diseased
% 
function [V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment) % complete parameters (works also for vectors)

V = mu.*(1-mu); % just standard
mz_twin_risk = 1 - v_environment ./ mu; % standard


