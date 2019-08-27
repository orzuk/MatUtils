% General (vectorised) FDR procedure
%
% Input:
%  P - matrix of SORTED(!) p-values (each row is a vector)
%  q - scalar controlling FDR output
%  str - name of procedure (case insensetive)
%  u - number of auxillary p-values to generate (optional)
%  m0 - how many true null hypothesis we have (optional). If we know m0 the
%         function also outputs V - the number of false positives. Note: We
%         assume that thet first m0 pvals are the null
%  is_sorted - flag saying that the p-values are already sorted, to save  sorting time (default is 'false')
%
%  Output:
%  R - matrix containing indices of hypotheses that passsed selected procedure
%  V - the number of false positives for each test (only outupt if we know m0)
%
function [R, V] = FDR_mat_main(P, q, str, u, m0, is_sorted, correction_factor, varargin)

if(~exist('u', 'var')) % here add auxillary i.i.d. uniform p-values
    u = [];
end
if(~exist('correction_factor', 'var'))
    correction_factor = 0;
end

if(~isempty(u)) % here add auxillary i.i.d. uniform p-values
    P(:,end+1:end+u) = rand(size(P,1), u); % add auxillary p-values
    R = FDR_mat_main(P, q, str);
else
    if(~exist('is_sorted', 'var') || isempty(is_sorted))
        is_sorted = 0;
    end
    if(~is_sorted)
        [P, ranks]= sort(P,2);
    else
        ranks = repmat(1:size(P,2), size(P,1), 1);
    end
    switch lower(str) % choose which procedure to use
        case 'ibh_up'
            R = fdr_IBH_mat(P,q, 'up');
        case 'ibh_down'
            R = fdr_IBH_mat(P,q, 'down');
        case 'ibh_log_up'
            R = fdr_IBHlog_mat(P, q, 'up', correction_factor);
        case 'ibh_log_down'
            R = fdr_IBHlog_mat(P, q, 'down', correction_factor);            
        case 'bh95'
            R = fdr_BH_mat(P,q); % This one is always step-up
        case 'sth'
            R = fdr_STH_mat(P,q);
        case 'bky'
            R = fdr_BKY_mat(P,q);
        case 'min_k'
            R = fdr_MINK_mat(P,q);
        otherwise
            disp('Unknown procedure.')
    end
end
if(exist('m0', 'var')) % here calculate also V
    if(~isempty(m0))
        V = zeros(length(R),1);
        pos_inds = find(R > 0);
        for i=pos_inds'
            V(i) = sum(P(i,ranks(i,:) <= m0)  <= P(i,R(i)));
        end
        %    V = sum(P(ranks(:,1:m0) <= repmat(R,1,m0),2);
    end
end


%FDR procedure for the case where the p values are based on statistically
%independent tests. input: P is a matrix of sorted p-values (each row is a vector), q is a scalar which
%control the false discovery rate, lambda is the parameter of Storey's procedure.
% output: F vector contain the number of
%hypotheses that pass the procedure. Improve Benjamini Hochberg 1995
function F=fdr_STH_mat(P,q, lambda, varargin)

if(~exist('lambda', 'var')) % set default value for lambda parameters
    lambda = 0.5;
end
[n,m]=size(P);
m0_hat = zeros(n,1);
for i=1:n
    m0_hat(i)=(m+1-find(P(i,:)<(1-lambda),1,'last'))/lambda;
end
F = FDR_adaptive_mat(P, m0_hat, q, step_str); % % %step down from the smaller p-value and up


%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a matrix of sorted p-values (each row is a vector), q is a scalar which
%control the false discovery rate. output: F vector contain the number of rejected hypotheses
%according to Benjamini Krieger Yekutieli 2006
function F=fdr_BKY_mat(P,q)

[n,m]=size(P);
F=zeros(n,1);
fdr_line=[1:m].*q./(m+1-[1:m]*(1-q));
for i=1:n
    in_fdr=find((P(i,:)-fdr_line)>0,1,'first');
    if in_fdr>1
        F(i)=in_fdr-1;
    elseif isempty(in_fdr)==1
        F(i)=m;
    else
        F(i)=0;
    end
end



%FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a matrix of sortted p-values (each row is a vector), q is a scalar which
%control the false discovery rate. step_str says if we perform a step-down or step-up procedure
%output: F vector contain the number of hypotheses that pass the procedure. Improve Benjamini Hochberg 1995.
function F=fdr_IBH_mat(P,q, step_str)

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
m0_hat=C*min(m*ones(n,1),max([s*ones(n,1),2*sum(P,2)],[],2));

F=FDR_adaptive_mat(P, m0_hat, q, step_str);



%FDR procedure for the log-IBH m0 estimator
%input: P is a matrix of sortted p-values (each row is a vector), q is a scalar which
%control the false discovery rate. step_str says if we perform a step-down or step-up procedure
%output: F vector contain the number of hypotheses that pass the procedure. Improve Benjamini Hochberg 1995.
function F=fdr_IBHlog_mat(P, q, step_str, correction_factor)

if(~exist('correction_factor', 'var'))
    correction_factor = 0;
end

m0_hat = m0_estimator(P, 'ibh_log', [], correction_factor);
%m0_hat=C*min(m*ones(n,1),max([s*ones(n,1),2*sum(P,2)],[],2));

F=FDR_adaptive_mat(P, m0_hat, q, step_str);

%Classic BH FDR procedure for the case when the p values are based on statistically
%independent tests. input: P is a matrix of sorted p-values, each row is a vector, q is a scalar which
%control the false discovery rate. output: F is a vector of the nuber of
%rejected hypotheses in each row
function F=fdr_BH_mat(P,q)
[n,m]=size(P);
F=FDR_adaptive_mat(P, m+zeros(n,1), q, 'up');



% This is not really a 'meaningful' procedure, and just used for
% simulations and comparision. The procedure simply gives the smallest
% q/m p-values. So the 'meaning' of q here is different from usual, as
% it is the fraction of rejected hypothesis
function F=fdr_MINK_mat(P,q)
[n,m]=size(P);

F = repmat(floor(q*m), n, 1); % set how many to reject

% A generic adaptive procedure when we know the estiator m0_hat (this has
% (to be estimated from somewhere else outside the function )
function F=FDR_adaptive_mat(P, m0_hat, q, step_str)

[n,m]=size(P);
F=zeros(n,1);
fdr_line=[1:m]*q/m;
if(strmatch(lower(step_str), 'down')) % % %step down from the smaller p-value and up
    for i=1:n
        in_fdr=find((P(i,:)-fdr_line.*(m/m0_hat(i)))>0,1,'first');
        if in_fdr>1
            F(i)=in_fdr-1;
        elseif isempty(in_fdr)==1
            F(i)=m;
        else
            F(i)=0;
        end
    end
else % step-up
    for i=1:n
        in_fdr=find((P(i,:)-fdr_line.*(m/m0_hat(i)))<=0,1,'last');
        if in_fdr>0
            F(i)=in_fdr;
        elseif isempty(in_fdr)==1
            F(i)=0;
        else
            F(i)=m;
        end
    end
end


