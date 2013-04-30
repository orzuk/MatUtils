% Computes R (reject. #) for a given vector of p-values and FDR level q.
% For a given desired False-Discovery-Rate (FDR) level q, the number R of
% rejected hypothesis having such an FDR level is the output
% Input:
%   P - Matrix of p-values (each row is a vector)
%   q - Scalar controlling FDR output
%
%   (optional) proc_str - Name of procedure (case insensetive)
%       Default procedure: 
%           'ibh_up' (Improved Benjamini Hochberg) 
%       Other Supported procedures: 
%           'ibh_down' (Improved Benjamini Hochberg step-down)
%           'bh95' (Standard Benhamini Hochberg)
%           'sth' (Storey's procedure with lambda=0.5)
%           'bky' (Benjamini-Krieger-Yekutieli)
%           'min_k' (Reject the lowest q fraction)
%   (optional) u - number of auxillary uniform p-values to generate 
%   (optional) m0 - how many true null hypothesis we have (optional). If we know m0 the
%       function also outputs V - the number of false positives. 
%       Note: If m0 is supplied as input, we assume that thet first m0 pvals are the null
%   (optional) is_sorted - flag saying that the p-values are already sorted, to save sorting time - Default is 'false'
%
% Output:
%   R - Matrix containing indices of hypotheses that passsed selected procedure
%   (optional) V - The number of false positives for each test (computed only when m0 is given as output)
%
%   Example:
%   P = rand(1,1000); % generate a vector of p-values, of which a 100 are taken from the alternative hypothesis
%   P(1:100) = P(1:100) ./ 2;
%   R = FDR_q_to_R(P, 0.05)  % Compute the number of rejected hypothesis with FDR level 0.05
%
function [R V] = FDR_q_to_R(P, q, proc_str, u, m0, is_sorted, varargin)

if(~exist('proc_str', 'var'))
    proc_str = 'ibh_up';
end  
if(isempty(proc_str))
    proc_str = 'ibh_up';
end
if(~exist('u', 'var')) % here add auxillary i.i.d. uniform p-values
    u = [];
end
if(~isempty(u)) % here add auxillary i.i.d. uniform p-values
    P(:,end+1:end+u) = rand(size(P,1), u); % add auxillary p-values
end
if(~exist('is_sorted', 'var'))
    is_sorted = 0;
end
if(~is_sorted)
    [P ranks]= sort(P,2);
else
    ranks = repmat(1:size(P,2), size(P,1), 1);
end

switch lower(proc_str) % choose which procedure to use    
    case {'ibh_up', 'ibh'}
        R = fdr_IBH(P,q, 'up');
    case 'ibh_down'
        R = fdr_IBH(P,q, 'down');
    case {'bh95', 'bh', 'bh_up'}
        R = fdr_BH(P,q); % This one is always step-up
    case {'sth', 'sth_up'}
        R = fdr_STH(P,q);
    case 'bky'
        R = fdr_BKY(P,q);
    case 'min_k'
        R = fdr_MINK(P,q);
    otherwise
        disp('Unknown procedure.')
end
if(exist('m0', 'var')) % here calculate also V
    if(~isempty(m0))
        V = zeros(length(R),1);
        pos_inds = find(R > 0);
        for i=pos_inds'
            V(i) = sum(P(i,ranks(i,:) <= m0)  <= P(i,R(i)));
        end
    end
end


%Storey's FDR procedure. input: P is a matrix of sorted p-values (each row is a vector), 
% q is a scalar which controls the false discovery rate, lambda is the parameter of Storey's procedure.
% output: F vector contain the number of hypotheses that pass the procedure. 
function F=fdr_STH(P,q, lambda, varargin)

if(~exist('lambda', 'var')) % set default value for lambda parameters
    lambda = 0.5;
end
if(~exist('step_str', 'var')) % set default value for step 
    step_str = 'up';
end
m0_hat = m0_estimator(P, 'sth', lambda);

F = FDR_adaptive(P, m0_hat, q, step_str); % % %step down from the smaller p-value and up


%Benjamini-Krieger-Yekutieli FDR procedure. 
%input: P is a matrix of sorted p-values (each row is a vector), q is a scalar which
%control the false discovery rate. output: F vector contain the number of rejected hypotheses
%according to Benjamini Krieger Yekutieli 2006
function F=fdr_BKY(P,q)

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



%Our Imptoved Benjamini-Hochberg FDR procedure. 
%input: P is a matrix of sortted p-values (each row is a vector), q is a scalar which
%control the false discovery rate. step_str says if we perform a step-down or step-up procedure
%output: F vector contain the number of hypotheses that pass the procedure.
function F=fdr_IBH(P,q, step_str)

m0_hat = m0_estimator(P, 'ibh', [],1); % estimate m0 with the correction factor
F=FDR_adaptive(P, m0_hat, q, step_str);

%Classic BH FDR procedure. 
%input: P is a matrix of sorted p-values, each row is a vector, q is a scalar which
%control the false discovery rate. output: F is a vector of the nuber of
%rejected hypotheses in each row
function F=fdr_BH(P,q)
[n,m]=size(P);
F=FDR_adaptive(P, m+zeros(n,1), q, 'up');



% This is not really a 'meaningful' procedure, and just used for
% simulations and comparision. The procedure simply gives the smallest
% q*m p-values. So the 'meaning' of q here is different from usual, as
% it is the fraction of rejected hypothesis
function F=fdr_MINK(P,q)
[n,m]=size(P);

F = repmat(floor(q*m), n, 1); % set how many to reject

% A generic adaptive procedure when we know the estiator m0_hat (this has
% to be estimated from somewhere else outside the function )
function F=FDR_adaptive(P, m0_hat, q, step_str)

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


