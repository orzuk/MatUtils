% Compute an indep. dist. with a given KL distance from uniform
% 
% Input: 
% nodes - number of (binary) random variables
% eps_vec - minimal error to tolerate
% 
% Output: 
% Q_dists - distribution with given KL distance
% 
function Q_dists = FindIndependentDistWithGivenKLFromUniform(nodes, eps_vec)

tol = 0.000000000001;
alpha_vec = zeros(1,length(eps_vec)); 

% Prepare things to return
Q_dists = zeros(length(eps_vec), nodes);
Q = rand(1,nodes); % Generate random  vector

while( log(2)*nodes +  sum(Q .* log(Q)) + sum((1-Q).*log(1-Q)) <= max(eps_vec))
    % Prepare first the direction at random
    Q = rand(1, nodes);
end

for i = 1:length(eps_vec)
    alpha_min = 0; 

    % Determine alpha_max
    alpha_max = 1; %%%min((1-P)./Q) - tol;
 
    Q_min = (1-alpha_min)*0.5+alpha_min * Q; Q_max = (1-alpha_max)*0.5+alpha_max * Q;
    
    cross_min =  log(2)*nodes +  sum(Q_min .* log(Q_min)) + sum((1-Q_min).*log(1-Q_min));
    cross_max =  log(2)*nodes +  sum(Q_max .* log(Q_max)) + sum((1-Q_max).*log(1-Q_max));

    while( (cross_max-cross_min) > tol)
        new_alpha = (alpha_min + alpha_max)/2;
        Q_new = (1-new_alpha)*0.5+new_alpha * Q;
        new_cross = log(2)*nodes +  sum(Q_new .* log(Q_new)) + sum((1-Q_new).*log(1-Q_new));
        
        if(new_cross < eps_vec(i))
            cross_min = new_cross;
            alpha_min = new_alpha;
        else
            cross_max = new_cross;
            alpha_max = new_alpha;
        end    
    end
    
    alpha_vec(i) = (alpha_min + alpha_max)/2;
    Q_dists(i,:) = (1-alpha_vec(i))*0.5 + alpha_vec(i) * Q;         
end
    
