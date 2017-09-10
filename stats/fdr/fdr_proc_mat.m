% FDR procedure with original BH procedure for matrics 
% for the case when the p values are based on statistically
% independent tests. input: P is a matrix of sorted p-values, each row is a vector, q is a scalar which
% control the false discovery rate. output: F is a vector of the nuber of
% rejected hypotheses in each row
%
function F=fdr_proc_mat(P,q)

[n,m]=size(P);
F=zeros(n,1);
fdr_line=[1:m]*q/m;
for i=1:n
    in_fdr=find((P(i,:)-fdr_line)<=0,1,'last');
    if in_fdr>0
        F(i)=in_fdr;
    elseif isempty(in_fdr)==1
        F(i)=0;
    else
        F(i)=m;
    end
end

