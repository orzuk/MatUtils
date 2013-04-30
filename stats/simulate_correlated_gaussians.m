% Simulates a set of correlated gaussian random
% variables with zero mean, unit variance and correlation coefficient rho
%
% Input:
% n - number of vectors
% m - dimension
% rho - correlation coefficient
% mu - vector of means (optional. Default is zero)
% 
% Output
% x - random draws
%
function x = simulate_correlated_gaussians(n, m, rho, mu)

x=randn(n,m+1); % temporarily randomize uncorrelated i.i.d. gaussians
if(length(rho)==1)
    if(rho == 0)
        x = x(:,1:m);
    else
        if(rho > 0) % when rho>0 just take weighted average with a common gaussian
            x = sqrt(rho).*repmat(x(:,m+1),1,m)+sqrt(1-rho).*x(:,1:m);
        else % when rho<0, simulate the hard way - generate i.i.d. gaussians and multiply by a matrix
            A = eye(m) .* (1-rho) + repmat(rho,m); B = sqrtm(A); x = (B * x(:,1:m)')';
        end
    end
else % here rho is input covariance matrix
    B = sqrtm(rho); x = (B * x(:,1:m)')';
end

if(exist('mu', 'var')) % change mean
    if(length(mu) == 1)
        mu = repmat(mu, n, 1);
    end
    x = x + repmat(vec2row(mu), n, 1);
end

    
    

