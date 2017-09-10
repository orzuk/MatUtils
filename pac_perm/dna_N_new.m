% Normalize genes in an array 
function [dat1n,mean_genes,norm_genes]=dna_N_new(dat1)

S=size(dat1,2);
mean_genes=mean(dat1,2);
dat1n=dat1-repmat(mean_genes,1,S);
clear dat1;
norm_genes=sqrt(sum(dat1n.^2,2));
norm_genes(find(norm_genes==0))=1;
dat1n=dat1n./repmat(norm_genes,1,S);

