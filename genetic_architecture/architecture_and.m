% Compute statistics for 'AND' of two architectures
function [mu V v_environment mz_twin_risk] = ...
    architecture_and(mu1, mu2, V1, V2, v_environment1, v_environment2) % Take and of two architectures
mu = mu1 .* mu2;
v_environment = mu - (mu1-v_environment1).*(mu2-v_environment2);
[V mz_twin_risk] = architecture_V_and_twin_risk(mu, v_environment);
