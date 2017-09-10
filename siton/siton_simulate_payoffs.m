% Compute recurrence formula for multi-step desicition game
%
% Input:
% num_turns - number of turns in the game
% payoff_dist - string representing the payoff distribution at each turn
% turn_cost - additional cost for each turn
% strategy - string representing the strategy chosen by the player
%
% Output:
% alpha_vec - vector of thresholds used by the strategy
% x_vec - vector of conditional expected payoffs from the strategy
%
function [alpha_vec x_vec] = siton_simulate_payoffs(num_turns, payoff_dist, turn_cost, strategy)

alpha_vec = zeros(num_turns, 1);
x_vec = zeros(num_turns, 1);


switch strategy
    case {'adaptive', 'optimal', 'smart'}
        alpha_vec(1) = icdf(payoff_dist, 0.00000000001, 0,1);
        x_vec(1) = diststat(payoff_dist, 1000); % Assume w.l.o.g. payoffs are ormalized to have mean one
        
        for i=2:num_turns
            turn_is = i
            B = cdf(payoff_dist, x_vec(i-1), 0, 1); % get general cumulative probability function
            if(i == 312)
                xxx = 234
            end
            A = cond_mean(payoff_dist, x_vec(i-1), 10000) * (1-B);
            x_vec(i) = A+B*x_vec(i-1);
        end
        
    case {'constant', 'constant-threshold', 'stupid'}
        
        
        % Here find alpha opt: the hard part
        max_alpha = icdf(payoff_dist, 0.9999, 0, 1); % take above median
        min_alpha = icdf(payoff_dist, 0.0001, 0, 1); % take below median
        mu = diststat(payoff_dist, 1000); 
        [alpha_opt opt_val] = fminbnd(@(x) stupid_bird(x,payoff_dist,mu,num_turns), min_alpha, max_alpha);
        [~, x_vec]  = stupid_bird(alpha_opt, payoff_dist, mu, num_turns)
        
%         while
%             alpha_vec(:) = opt_alpha; % assume we've found the best alpha
%             B = cdf(payoff_dist, opt_alpha); % get general cumulative probability function
%             A = cond_mean(payoff_dist, alpha_opt) * (1-B);
%             x_vec = A .* (1-B.^(1:num_turns)) ./ (1 - B); % set payoff
%         end
        
    case {'top_k', 'top', 'secretary'}
        
        
        
        
end

function [x x_vec]  = stupid_bird(alpha, payoff_dist, mu, num_turns)
alpha_vec = zeros(num_turns, 1);
alpha_vec(:) = alpha; % assume we've found the best alpha
B = cdf(payoff_dist, alpha, 0, 1); % get general cumulative probability function
A = cond_mean(payoff_dist, alpha) * (1-B);
x_vec = A .* (1-B.^(0:(num_turns-1))) ./ (1 - B) + mu*B.^(0:(num_turns-1)); % set payoff
x = -x_vec(end); % payoff when N turns are left. Take minus sign since we do minimizer
