% Compute a distribution with a given KL distance from P 
function Q_dists = FindDistWithGivenKL(P, eps_vec)

tol = 0.000000000001;

alpha_vec = zeros(1,length(eps_vec)); 

% Prepare things to return
Q_dists = zeros(length(P), length(eps_vec));

Q=P;

while(cross_entropy(Q, P) <= max(eps_vec))
    % Prepare first the direction at random
    Q = rand(length(P),1);
    Q = Q ./ sum(Q); % Normalize        
end

for i = 1:length(eps_vec)
    alpha_min = 0; 

    % Determine alpha_max
    alpha_max = 1; %%%min((1-P)./Q) - tol;

    cross_min = cross_entropy((1-alpha_min)*P+alpha_min * Q, P);
    cross_max = cross_entropy((1-alpha_max)*P+alpha_max * Q, P);

    while( (cross_max-cross_min) > tol)
        new_alpha = (alpha_min + alpha_max)/2;
        new_cross = cross_entropy((1-new_alpha)*P+new_alpha * Q, P);
        
        if(new_cross < eps_vec(i))
            cross_min = new_cross;
            alpha_min = new_alpha;
        else
            cross_max = new_cross;
            alpha_max = new_alpha;
        end    
    end
    
    alpha_vec(i) = (alpha_min + alpha_max)/2;
    Q_dists(:,i) = (1-alpha_vec(i))*P + alpha_vec(i) * Q;         
end
    
