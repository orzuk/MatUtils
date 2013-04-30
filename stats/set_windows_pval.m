function pvals_sorted = set_windows_pval(vals, N, min_dist, FDR_thresh)

% Assume n variables are uniformly distributed on [1, N]
% The function calculates the p-value of the 'windows', i.e. the distances
% between any two values, to be smaller than they really are.
% This can be useful for finding 'significant' 'dense' windows
% min_dist is the minimal 'distance' to start check. If 1, we check couples
% (neighbours). If 2, we check only couples for which there is at least one
% between them etc etc. It is recommended to take min_dist >= 3, since for
% smaller values we are not likely to get a significant pvalue anyway, and
% thus we only introduce more 'noise' to the FDR.
% FDR_thresh is the FDR threshold to be assigned to the pvalues after we
% have found them. If you don't want to use FDR, but get all the values,
% put FDR_thresh = -1
% Output a 3*r matrix, where r pvalues were found significant.
% The format is : [ left_side, right_side, pvalue ]


epsilon = 0.00000000000001; % Just for matlab f***ing warnings ...

% First set number of variables
n = length(vals);


if(min_dist >= n)
    pvals_sorted = zeros(1,3);
    pvals_sorted(1,1) = min(vals); pvals_sorted(1,2) = max(vals); pvals_sorted(1,3) = 0.999; 
else
    
    % Determine maximal window size
    %%%max_window = vals(n) - vals(1);
    
    % f is the local density function of the order statistics.
    % f(k, x) is the probability : Pr(x_(k) = x), where x_(k) is the k-th 
    % order statistic. (k=1 the minimum and k=n the maximum)
    f_order = zeros(n, N);
    
    
    p_vec0 = (1.0/N) * [1:N] - epsilon;
    p_vec1 = (1.0/N) * [0:N-1] + epsilon;
    n_vec = n*ones(1, N);
    
    
    % % % for k=1:n   
    % % % 
    % % %  %%%   x_vec = (k-1)*ones(1,N);   
    % % %     x_vec = repmat(k-1,1,N);   
    % % %     f_order(k,:) = binocdf(x_vec, n_vec, p_vec1) - binocdf(x_vec, n_vec, p_vec0);   
    % % % end
    
    
    %%% different coputation ....
    % k=1 : 
    x_vec = repmat(0,1,N);  
    f_order(1,:) = binocdf(x_vec, n_vec, p_vec1) - binocdf(x_vec, n_vec, p_vec0);  
    
    for k=2:n   
        x_vec = repmat(k-1,1,N);   
        f_order(k,:) = f_order(k-1,:) + binopdf(x_vec, n_vec, p_vec1) - binopdf(x_vec, n_vec, p_vec0);   
    end
    
    
    
    % Calculate the accumulative distributions of the order statistics
    % F is the accumulative distribution function of the order statistics.
    % F(k, x) is the probability : Pr(x_(k) <= x), where x_(k) is the k-th 
    % order statistic. (k=1 the minimum and k=n the maximum)
    F_order = cumsum(f_order, 2);
    
    
    timtim = cputime;
    
    
    % Prepare an inverse vector to increase performance ..
    inverse_vec = 1 ./ [N:-1:1];
    
    % % % pvals = 2*ones(n); 
    num_couples = binom(n-min_dist+1, 2);
    pvals_sorted = zeros( num_couples, 3);
    counter = 1;
    
    
    % Now do the main loop, go over all the couples
    for k=1:n-min_dist
        for m=k+min_dist:n
            y = vals(m) - vals(k);
            if(y < N*(m-k)/(n+1))  % If the 'window' is greater then expected there is nothing to check ...
                do_flag = 1; % for now we do ...
                for j=m+1:n
                    z =  vals(j) - vals(k);
                    if(y >= (m-k)*z/(j-k))
                        do_flag = 0;
                        break;
                    end
                end
                
                if(do_flag)   % Now check the other side for the same thing ..
                    for j=1:k-1
                        z =  vals(m) - vals(j);
                        if(y >= (m-k)*z/(m-j))
                            do_flag = 0;
                            break;
                        end
                    end 
                end
                
                if(do_flag)  % No point in checking 'windows' if they are containted in more significant windows ..
                    % Do the same in a vector form
                    %%%                x_vec = (m-k-1)*ones(1, N-y); n_vec = (n-k)*ones(1, N-y);
                    x_vec = repmat(m-k-1,1, N-y); n_vec = repmat(n-k,1, N-y);
                    %%%       p_vec = N+1-[1:N-y]; p_vec = (y+1) ./ p_vec;
                    p_vec = (y+1) .* inverse_vec(1:N-y); 
                    
                    F_m_given_k_vec = 1 - binocdf(x_vec, n_vec, p_vec);
                    %%%   pvals(k, m) = 1-F_order(k, N-y) + sum(f_order(k, 1:N-y) .* F_m_given_k_vec); 
                    
                    pvals_sorted(counter,1) = vals(k); pvals_sorted(counter,2) = vals(m); 
                    pvals_sorted(counter,3) = 1-F_order(k, N-y) + sum(f_order(k, 1:N-y) .* F_m_given_k_vec); counter=counter+1;
                    %%%pvals_sorted(counter,3) = pvals(k, m); 
                end
            end
        end
    end
    
    
    
    maintim = cputime - timtim;
    
    
    %%%pvals_sorted = pvals_sorted(1:counter-1,:);
    [sorted sort_perm] = sort(pvals_sorted(1:counter-1,3)); 
    
    pvals_sorted = pvals_sorted(sort_perm,:);
    
    
    
    % Now do the FDR procedure 
    if(FDR_thresh > -1)
        max_index = 0;
        for i=1:num_couples
            if(pvals_sorted(i,3) < (FDR_thresh*i/num_couples) )
                max_index = i;
            end
        end
        pvals_sorted = pvals_sorted(1:max_index,:);
    end    
    
    
    if(size(pvals_sorted, 1)==0 )                   
        pvals_sorted = zeros(1,3);
        pvals_sorted(1,1) = min(vals); pvals_sorted(1,2) = max(vals); pvals_sorted(1,3) = 0.999; 
    end
    
    
end 

function [m] = binom(n, k)
    
m = nchoosek(n, k);
    


