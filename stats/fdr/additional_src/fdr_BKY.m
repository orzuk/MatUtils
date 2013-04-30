% FDR procedure using Benjamini-Kreiger-Yekutieli method 
% for the case when the p values are based on statistically
% independent tests. input: P is a vector of p-values, q is a scalar which
% control the false discovery rate. output: F the indices of the element
% that pass the procedure. Benjamini Krieger Yekutieli 2006
%
function F=fdr_BKY(P,q)

[m,n]=size(P);
if n>m; P=P'; end
[P_sor]=sort(P);
N=length(P);
in_fdr=find((P_sor-[1:N]'.*q./(N+1-[1:N]'*(1-q)))>0,1,'first');
if in_fdr>1
    p_d=P_sor(in_fdr-1);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[1:N];
else
    F=[];
end

