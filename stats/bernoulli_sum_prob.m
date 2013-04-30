% Compute the probability distribution of sum of N independent bernoulli random
% variables, with probabilities f_vec(i) of being one. 
% The complexity depends on the number of different probabilities. If
% they're all different the complexity is exponential in N 
function bernoulli_sum_probs  = bernoulli_sum_prob(f_vec)

N = length(f_vec); % get number of variables 
bernoulli_sum_probs = zeros(N+1,1); bernoulli_sum_probs(1) = 1; % Start with all mass at zero 

[f f_counts] = unique_with_counts(f_vec); 
if(length(f) == 1) % all probabilities are the same 
    bernoulli_sum_probs = binopdf(0:N, N, f);
else % here we've got multiple probabilities 
    for i=1:length(f) % loop over different probabilities (if they're all different it can be exponential ...) 
        tmp_bernoulli_sum_probs = binopdf(0:f_counts(i), f_counts(i), f(i));
        bernoulli_sum_probs = conv(bernoulli_sum_probs, tmp_bernoulli_sum_probs); 
        bernoulli_sum_probs = bernoulli_sum_probs(1:N+1); % cut tails 
    end
    if(isrowvector(f_vec))
        bernoulli_sum_probs = vec2row(bernoulli_sum_probs);
    else
        bernoulli_sum_probs = vec2column(bernoulli_sum_probs);
    end
end
