% FDR procedure with original BH procedure 
% for the case when the p values are based on statistically
% independent tests. input: P is a vector of p-values, q is a scalar which
% control the false discovery rate. output: F the indices of the element
% that pass the procedure
%
function F=fdr_proc(P,q)

[m,n]=size(P);
if n>m; P=P'; end
N=length(P);
[P_sor]=sort(P);
in_fdr=find((P_sor-[1:N]'*q/N)<=0,1,'last');
if in_fdr>0
    p_d=P_sor(in_fdr);
    F=find(P<=p_d);
elseif isempty(in_fdr)==1
    F=[];
else
    F=[1:N];
end

