% Clipp gaussian variables into discrete variables by finding the best
% threshold which seperates them.
% We assume mu1 < mu2 , so we take the right area of the 1st gaussian and
% the left area of the 2nd gaussian   

%TH(i)=fzero(inline('sum(sum((M).*(exp(x).^(M)))./sum(exp(x).^(M)),2)-s','x','M','s'),...
% (score_vec(i)-clt_mean)/clt_var,optimset('Display','off'),M,score_vec(i));
function Threshold = TwoGaussiansThreshold( mu1, s1, mu2, s2)

Threshold = fzero(@TwoGaussiansDiff, (mu1-mu2)/2, optimset('Display','off'), mu1, s1, mu2, s2);
