% A temp script preparing gene-expression data for plotting
cd /a/fs-01/pdoc/mrosenz/Liat/GeneCorr/data
% cd /a/fs-01/pdoc/mrosenz/Liat
KNN = 10; % parameter for knn missing values algorithm
fraction=0.01;
saf=0.05;
% load WANG_DATA
% data=gene_expression_log2;
% samples=Prognosis;
% remove_absent_genes
load GISETTE_data_file
data=T;
clear T
samples=single(Labels==-1);
clear Labels
% clear all
% load Rosetta_data.mat
% clear samples
% samples=real(Labels<100);
% remove_absent_genes

