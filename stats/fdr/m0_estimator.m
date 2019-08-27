% Estimate m0, the number of true null hypotheses, 
% from the vector of p-values.
%
% Input: 
% P - vector of p-values
% proc_str - string specifying type of estimattion procedure
%       Default procedure: 
%           'ibh_up' (Improved Benjamini Hochberg) 
%       Other Supported procedures: 
%           'ibh_down' (Improved Benjamini Hochberg step-down)
%           'bh95' (Standard Benhamini Hochberg)
%           'sth' (Storey's procedure with lambda=0.5)
%           'bky' (Benjamini-Krieger-Yekutieli)
%           'min_k' (Reject the lowest q fraction)
%
% (optional) lambda - additional procedure-specific parameters (e.g. lambda value of Storey's procedure)
% (optional) correction_factor - whether or not to correct the m0 estimator to satisfy certain bounds(e.g. m0 <= m)
%
% Output: 
% m0 - estimator of the number of true null hypothesis
%
%   Example: 
%   P = rand(1,1000); % generate a vector of p-values, of which a 100 are taken from the alternative hypothesis
%   P(1:100) = P(1:100) ./ 200; 
%   m0 = m0_estimator(P);  % estimate m0 (true m0 is 900) 
%
function m0 = m0_estimator(P, proc_str, lambda, correction_factor, varargin)

if(~exist('proc_str', 'var'))
    proc_str = 'ibh_up';
end  
if(isempty(proc_str))
    proc_str = 'ibh_up';
end
if(~exist('correction_factor', 'var'))
    correction_factor = 0;
end
switch lower(proc_str) % choose estimator
    case {'sum_pi', 'ibh_up', 'ibh_down', 'ibh'}
        m0 = m0_estimator_ibh(P, correction_factor);
    case {'bh', 'bh95'}
       m0 = size(P,2);
    case 'sth'
        m0 = m0_estimator_sth(P, lambda);
    case 'bky' % here there's no meaning to the estimator
        m0 = 9999;
    case {'sum_log_pi', 'ibh_log'}
        m0 = 2*correction_factor-sum(log(1-P));
    otherwise
        disp('Unknown procedure.')
end

function m0 = m0_estimator_ibh(P, correction_factor)
if(correction_factor)
    [n,m]=size(P);
    m_C_s=1e5 *[     0.001000000000000   0.000010306035633   0.000350000000000
        0.002000000000000   0.000010199148511   0.000550000000000
        0.003000000000000   0.000010156707816   0.000720000000000
        0.004000000000000   0.000010132670482   0.000860000000000
        0.005000000000000   0.000010117091093   0.000980000000000
        0.006000000000000   0.000010105538438   0.001090000000000
        0.007000000000000   0.000010096875092   0.001190000000000
        0.008000000000000   0.000010090004408   0.001290000000000
        0.009000000000000   0.000010084411050   0.001380000000000
        0.010000000000000   0.000010079681966   0.001470000000000
        0.020000000000000   0.000010054904854   0.002170000000000
        0.030000000000000   0.000010044257282   0.002720000000000
        0.040000000000000   0.000010038076350   0.003180000000000
        0.050000000000000   0.000010033856638   0.003590000000000
        0.060000000000000   0.000010030847440   0.003960000000000
        0.070000000000000   0.000010028436400   0.004300000000000
        0.080000000000000   0.000010026535868   0.004620000000000
        0.090000000000000   0.000010025022004   0.004910000000000
        0.100000000000000   0.000010023656804   0.005210000000000
        0.150000000000000   0.000010019218320   0.006450000000000
        0.200000000000000   0.000010016616056   0.007500000000000
        0.250000000000000   0.000010014818536   0.008430000000000
        0.300000000000000   0.000010013492425   0.009280000000000
        0.400000000000000   0.000010011679484   0.010770000000000
        0.500000000000000   0.000010010408058   0.012110000000000
        0.600000000000000   0.000010009486880   0.013320000000000
        0.700000000000000   0.000010008794243   0.014390000000000
        0.800000000000000   0.000010008212513   0.015430000000000
        0.900000000000000   0.000010007735493   0.016410000000000
        1.000000000000000   0.000010007339726   0.017310000000000];
    C=interp1(m_C_s(:,1),m_C_s(:,2),m,'spline');
    s=interp1(m_C_s(:,1),m_C_s(:,3),m,'spline');
    m0=C*min(m*ones(n,1),max([s*ones(n,1),2*sum(P,2)],[],2));    
else
    m0 = 2*sum(P,2);
end

 %Estimate m0 with Storey's approach
function m0 = m0_estimator_sth(P, lambda, varargin)

if(~exist('lambda', 'var')) % set default value lambda=0.5
    lambda = 0.5;
end
[n,m]=size(P);
m0 = zeros(n,1); 
for i=1:n
    m0(i)=(m+1-find(P(i,:)<(1-lambda),1,'last'))/lambda;
end

