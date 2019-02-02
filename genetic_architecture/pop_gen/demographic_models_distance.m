% Compute how similar are two demographic models. Allow stretching of time (?)
% Input:
% N_vec1 - population size at each time for first model
% N_vec2 - population size at each time for second model
%
% Output:
% d - average difference in (log) population size between two models
%
function d = demographic_models_distance(N_vec1, N_vec2)

t1 = length(N_vec1); t2 = length(N_vec2);

if(t1>t2)
    d = demographic_models_distance(N_vec2, N_vec1);
else % here we know that t2>t1
    stretch_time=0;
    if(stretch_time) % compre all times stretched
        N_vec1 = interp1((1:t1).*t2./t1, N_vec1, (1:t2)', 'linear', N_vec1(1));
        d = mean((log(N_vec1)-log(N_vec2)).^2);
    else % compare only recent times
        d = mean((log(N_vec1)-log(N_vec2((t2-t1+1):end))).^2);
    end
end
