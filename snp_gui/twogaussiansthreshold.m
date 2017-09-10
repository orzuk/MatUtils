% Clipp gaussian variables into discrete variables 
% We assume mew1 < mew2 , so we take the right area of the 1st gaussian and
% the left area of the 2nd gaussian   

%TH(i)=fzero(inline('sum(sum((M).*(exp(x).^(M)))./sum(exp(x).^(M)),2)-s','x','M','s'),...
% (score_vec(i)-clt_mean)/clt_var,optimset('Display','off'),M,score_vec(i));
function Threshold = TwoGaussiansThreshold( mew1, s1, mew2, s2)

Threshold = fzero(@TwoGaussiansDiff, (mew1-mew2)/2, optimset('Display','off'), mew1, s1, mew2, s2);
