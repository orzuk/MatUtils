% get distribution of allelic ages from array of allele frequencies 
%
% Input: simulation_struct - structure with simulations results 
% 
% Output: 
% allele_age_vec - vector of allelic ages 
% 
function allele_age_vec = allele_age_from_sim_struct(simulation_struct)

[m, n] = size(simulation_struct.q); allele_age_vec=n+zeros(m,1); 
for i=1:m
    t=find(simulation_struct.q(i,:) == 0, 1, 'last');
    if(~isempty(t))
        allele_age_vec(i) = n-t; % find(simulation_struct.q(i,:) == 0, 1, 'last');
    end
end

