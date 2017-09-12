% FDR procedure using Benjamini-Kreiger-Yekutieli method for matrices
% for the case when the p values are based on statistically
% independent tests. input: P is a matrix of sorted p-values (each row is a vector), q is a scalar which
% control the false discovery rate. output: F vector contain the number of rejected hypotheses
% according to Benjamini Krieger Yekutieli 2006
%
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

