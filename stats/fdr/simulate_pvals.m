% Simulate a set of p-values
%
% Input:
% v - number of vectors of pvals
% m - number of pvals (hypothesis) in each vector
% m0 - number of null-pvals
% mu_gap - measure of difference between the null and non-null p-vals
% dependency_flag - what kind of dependency we want for the pvals (local/global)
% rho - measure of correlation between the null p-vals
% two_sided_flag - flag saying if the test is one-sided (default) or two sided 
% sort_flag - if this is on, we sort the p-values (default is false)
%
% Output:
% pvals_vec - a matrix of size v*m where each vector is a set of p-values
%
function pvals_vec = simulate_pvals(v, m, m0, mu_gap, dependency_flag, rho, two_sided_flag, sort_flag, varargin)

if(~exist('two_sided_flag', 'var'))
    two_sided_flag = 0;
end
if(isempty(two_sided_flag))
    two_sided_flag = 0;
end
if(~exist('sort_flag', 'var'))
    sort_flag = 0;
end

switch dependency_flag
    case 0  % global. Works only for positive rho for now
        x = simulate_correlated_gaussians(v,m,rho);
        x(:,m0+1:m) = x(:,m0+1:m) + mu_gap; % update the alternative hypothesis biased values
%        corr(x(:,1), x(:,2))
        if(two_sided_flag)
            pvals_vec=2*normcdf(-abs(x),0,1); % standard Z to p transofmation
        else
            pvals_vec=1-normcdf(x,0,1);   % one-sided transofmration keeping the sign
        end
            
    case 1 % here we've got one dependency for the null and one for the alternative, so rho must be of size 2!!!
        x = zeros(v,m);
        x(:,1:m0) = simulate_correlated_gaussians(v,m0,rho(1));
        x(:,m0+1:m) = simulate_correlated_gaussians(v,m-m0,rho(2)) + mu_gap;
        pvals_vec=1-normcdf(x,0,1); % one-sided transofmration keeping the sign (too lasy to list all configurations of rho's sign) 
    
    case 2 % Shlomo's idea for [0,1/2] [1/2 0] correlation. Work's currently only on two p-values
        x = rand(v,m); 

        if(two_sided_flag == 1)
            small_inds = find(x(:,1) < rho); big_inds = find(x(:,1) > 1-rho); med_inds = setdiff(1:v, union(small_inds,big_inds));
            x(big_inds,2) = 1 - x(big_inds,1);
            x(small_inds,2) = 1 - x(small_inds,1);
            %        x(med_inds,2) = (1-2.*rho) .* x(med_inds,2) + rho; % independent of others
            x(med_inds,2) = x(med_inds,1); % same as others
        else % one sided p-vals
            %            x(:,2) = 1-x(:,1);
            small_inds = find(x(:,1) < rho); big_inds = setdiff(1:v, small_inds);
            x(big_inds,2) = x(big_inds,1);
            x(small_inds,2) = rho .* x(small_inds,2);
            % % %             %        x(med_inds,2) = (1-2.*rho) .* x(med_inds,2) + rho; % independent of others
        end
        pvals_vec = x;
        
end

if(sort_flag) % sort pvalues
   pvals_vec = sort(pvals_vec,2); 
end
    
    

