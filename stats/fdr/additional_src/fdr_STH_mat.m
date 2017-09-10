% FDR procedure with Storey's method for matrices
% for the case when the p values are based on statistically
% independent tests. input: P is a matrix of sortdet p-values (each row is a vector), q is a scalar which
% control the false discovery rate. output: F vector contain the number of
% hypotheses that pass the procedure. Improve Benjamini Hochberg 1995
%
function F=fdr_STH_mat(P,q)

[n,m]=size(P);



% % %step down from the smaller p-value and up
F=zeros(n,1);
fdr_line=[1:m]*q/m;
for i=1:n
    m0_hat=(m+1-find(P(i,:)<0.5,1,'last'))/0.5;
    in_fdr=find((P(i,:)-fdr_line*m/m0_hat)>0,1,'first');
    if in_fdr>1
        F(i)=in_fdr-1;
    elseif isempty(in_fdr)==1
        F(i)=m;
    else
        F(i)=0;
    end
end
