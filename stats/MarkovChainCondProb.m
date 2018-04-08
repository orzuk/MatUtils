% Compute conditional probability to be in state i given that we didn't visit any state in A
% Input: 
% p0 - initial condition
% M - transition matrix 
% A - forbidden set
% method_str - 'numeric' exact computation, 'simulation' 
% iters - number of iterations to simulate 
% 
% Output: 
% p - vector of probabilities at time T with p_i = Pr(X_T=i | X_1,...,X_T \in A^c)
%
function p = MarkovChainCondProb(p0, M, A, T, method_str, iters)

n = length(M); % get number of states 
Ac = setdiff(1:n, A); % get complementary states 
p_to_Ac = sum(M(:,Ac), 2); % probability to get to A from each state 

P_Ac = zeros(n, T+1); % Prob(x_t = i , x_1, ..., x_t \in Ac)
P_Ac(:,1) = p0; % set init conditions 

if(~exist('iters', 'var') || isempty(iters))
    iters=10000;
end

switch method_str
    case 'numeric'        
        for t=1:T % recursion
            for j=Ac % loop on set of allowed states
                for i=Ac
                    P_Ac(j, t+1) = P_Ac(j, t+1) + P_Ac(i, t) * M(i,j); %                p_new(j) = p_old(i) * ??; 
                end
            end
        end
        P_Ac_sum = sum(P_Ac); 
        P_Ac_cond = P_Ac ./ repmat(P_Ac_sum, n, 1); 
        p = P_Ac_cond(:,T+1)'; % return last probabilities 
    case 'simulation'
        mc = dtmc(M);
        x = simulate(mc, T, 'X0', iters*p0); % simulate vector 
%        x = x(:,find(1-max(ismember(x, A)))); % remove chains passing through forbidden state
        p = hist(x(end,find(1-max(ismember(x, A)))), 1:n); p = p ./ sum(p); 
        
    case 'sim2' % simulate and reject on the fly 
        w = ones(iters, 1); % set weights vector 
        x = weighted_rand(p0, iters); x_new=zeros(iters, 1); % get initial samples 
        for t=1:T % loop on times 
            for i=Ac
               I = find(x == i);  
               x_new(I) = Ac(weighted_rand(M(i,Ac)./p_to_Ac(i), length(I))); 
               w(I) = w(I) * p_to_Ac(i); 
            end
            x = x_new; % update x 
        end
        p = weighted_hist(x, w, 1:n); p = p ./ sum(p); 
end



