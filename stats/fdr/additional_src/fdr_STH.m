% FDR procedure with Storey's method
% for the case when the p values are based on statistically
% independent tests. input: P is a vector of p-values, q is a scalar which
% control the false discovery rate. output: F the indices of the element
% that pass the procedure. Storey, Taylor, Siegmund 2004
function F=fdr_STH(P,q)

[m,n]=size(P);
if n>m; P=P'; end
N=length(P);

m0=(N+1-length(find(P<0.5)))/0.5;

[P_sor]=sort(P);
% % %step down from the smaller p-value and up
in_fdr=find((P_sor-[1:N]'*q/m0)>0,1,'first');
if in_fdr>1
    p_d=P_sor(in_fdr-1);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[1:N];
else
    F=[];
end
