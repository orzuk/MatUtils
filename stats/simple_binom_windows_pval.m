function pvals_sorted = simple_binom_windows_pval(vals, N, min_dist, FDR_thresh)

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
    
    
    num_couples = binom(n-min_dist+1, 2);
    pvals_sorted = zeros( num_couples, 3);
    input_to_binom = zeros( num_couples, 2); % No need for n
    counter = 1;
    
    
    % Now do the main loop, go over all the couples
    for k=1:n-min_dist
        for m=k+min_dist:n
            y = vals(m) - vals(k);
            if(y*(n+1) < N*(m-k))  % If the 'window' is greater then expected there is nothing to check ...
                do_flag = 1; % for now we do ...                
                
                % an extra check from both sides. Is it usefull ?
                for j=k:-1:1
                    for jj=m:n
                        z = vals(jj) - vals(j) + epsilon;
                        if(y*(jj-j) >=  (m-k)*z)
                            do_flag = 0; 
                            break;
                        end
                    end
                    if(do_flag == 0)
                        break;
                    end
                end

            
                
                if(do_flag)  % No point in checking 'windows' if they are containted in more significant windows ..
                    
                    pvals_sorted(counter,1) = k; pvals_sorted(counter,2) = m; 
                    %%%input_to_binom(counter,1) = m-k+1;  
                    
                    input_to_binom(counter,2) = (y+1);
                    counter=counter+1;
% % %                      pvals_sorted(counter, 3) =  - binocdf(m-k+1, n, (y+1)/N) + binopdf(m-k+1, n, (y+1)/N); counter = counter+1;
                end
            end
        end
    end
            
    pvals_sorted = pvals_sorted(1:counter-1,:);
    input_to_binom = input_to_binom(1:counter-1,:);
    input_to_binom(:,1) = pvals_sorted(:,2) - pvals_sorted(:,1) + 1;  % m-k+1
    
    if(sum( size(pvals_sorted(:,1)) ==  size(vals(pvals_sorted(:,1))) ) == 2)  
        pvals_sorted(:,1) = (vals(pvals_sorted(:,1)));   pvals_sorted(:,2) = (vals(pvals_sorted(:,2)));  % Take the vals
    else
        pvals_sorted(:,1) = (vals(pvals_sorted(:,1)))';   pvals_sorted(:,2) = (vals(pvals_sorted(:,2)))';  % Take the vals
    end
    input_to_binom(:,2) = input_to_binom(:,2) ./ N;  % Normalize by N
    
 
    pvals_sorted(1:counter-1, 3)  = 1 - binocdf(  input_to_binom(:,1),n,input_to_binom(:,2) ) + ...
        binopdf(  input_to_binom(:,1),n,input_to_binom(:,2) );
 

    
    
    [sorted sort_perm] = sort(pvals_sorted(:,3)); 
    pvals_sorted = pvals_sorted(sort_perm,:);
    
    
 % Normalization !!!   
%%%    (length(find(big_simple_binom_pvals_table(:,10) <= myp))+0.5)/iters
    
    
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
    
    % correction to avoid empty output
    if(size(pvals_sorted, 1)==0 )                   
        pvals_sorted = zeros(1,3);
        pvals_sorted(1,1) = min(vals); pvals_sorted(1,2) = max(vals); pvals_sorted(1,3) = 0.999; 
    end
        
end 

function [m] = binom(n, k)
    
m = nchoosek(n, k);
    
