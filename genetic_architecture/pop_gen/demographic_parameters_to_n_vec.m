% Convert a set of demographic parameters into a vector of population sizes
% Input: 
% D - a structure with demographic parameters 
% i - index stating which demographic parameters to take  
%
% Output: 
% N_vec - a vector of population size at each generation 
% 
function N_vec = demographic_parameters_to_n_vec(D, i)

if(~exist('i', 'var') || isempty(i)) % default is first 
    i=1;
end
if(~isfield(D, 'num_params_vec'))
    D.num_params_vec = ones(3,1);
end
i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
if(~isfield(D, 'generations')) 
    D.generations = zeros(1, D.num_stages); D.expan_rate = zeros(1, D.num_stages); D.init_pop_size = zeros(1, D.num_stages);
    for j=1:D.num_stages
        D.init_pop_size(j) = D.init_pop_size_vec{j}(i_vec(j));
        D.generations(j) = D.generations_vec{j}(i_vec(j+D.num_stages));
        D.expan_rate(j) = D.expan_rate_vec{j}(i_vec(j+2*D.num_stages));
    end
end

% Here we expand: gen_vec, expan_rate_vec, init_pop_size_vec  to one vector N_vec containing population size at each generation
num_generations = sum(D.generations); 

N_vec = zeros(num_generations, 1); 

ctr=0;
for j=1:length(D.generations) % loop on #generation
   if(D.init_pop_size(j) == -1)
       init_size = N_vec(ctr); 
   else
       init_size = D.init_pop_size(j); 
   end
   N_vec((ctr+1):(ctr+D.generations(j))) = ceil(init_size * D.expan_rate(j) .^ (0:D.generations(j)-1)); % get (rounded) population size 
   ctr = ctr + D.generations(j); 
end
