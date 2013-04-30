% Clipp gaussian variables into discrete variables 
% We assume mew1 < mew2 , so we take the right area of the 1st gaussian and
% the left area of the 2nd gaussian   

%TH(i)=fzero(inline('sum(sum((M).*(exp(x).^(M)))./sum(exp(x).^(M)),2)-s','x','M','s'),...
% (score_vec(i)-clt_mean)/clt_var,optimset('Display','off'),M,score_vec(i));
function IntegralsDiff = TwoGaussiansDiff( x, mew1, s1, mew2, s2)
minmew = min(mew1, mew2);
maxmew = max(mew1, mew2);
mins = s1; 
maxs = s2;
if(minmew == mew2)
    mins = s2; 
    maxs = s1; 
end
IntegralsDiff = 1-normcdf(x,minmew,mins) - normcdf(x,maxmew,maxs);
