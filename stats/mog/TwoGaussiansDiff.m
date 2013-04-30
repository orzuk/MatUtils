% Find the difference in the tail areas for two gaussians from a certain point x
%
% Input:
% x - threshold
% mu1 - first gauss. mean
% s1 - first gauss. sigma
% mu2 - second gauss. mean
% s2 - second gauss. sigma
% direction - should we take opposite directions (default, 0) right (1) or left (2)
%
% Output: 
% IntegralsDiff - the tail areas difference 
% 
%TH(i)=fzero(inline('sum(sum((M).*(exp(x).^(M)))./sum(exp(x).^(M)),2)-s','x','M','s'),...
% (score_vec(i)-clt_mean)/clt_var,optimset('Display','off'),M,score_vec(i));
function IntegralsDiff = TwoGaussiansDiff( x, mu1, s1, mu2, s2, direction, varargin)
if(~exist('direction', 'var'))
    direction = 0;
end

min_mu = min(mu1, mu2);
max_mu = max(mu1, mu2);
min_s = s1;
max_s = s2;
if(min_mu == mu2)
    min_s = s2;
    max_s = s1;
end
switch direction
    case 0 % opposite directions (default)
        IntegralsDiff = 1-normcdf(x,min_mu,min_s) - normcdf(x,max_mu,max_s);
    case 1 % take right tail
        IntegralsDiff = normcdf(x,min_mu,min_s) - normcdf(x,max_mu,max_s);
    case 2 % take left tail
        IntegralsDiff = -normcdf(x,min_mu,min_s) + normcdf(x,max_mu,max_s);
end
