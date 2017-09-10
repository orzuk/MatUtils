% Written by Or Zuk 6/2007
%
% This function plots an image showing, for each of the 36
% different states, what does it represent, in terms of
% copy number and genotype for each chromosome, and for A and B
%
function DrawHMMStatesImage()


x_dim = 3; x_dim_eq = 2*x_dim-1;

% Tables for indexing the A and B copy 
A_copy_tab = [ 0, 0, 1, 0, 2, 0, 0, 0, 1, ...
               0, 2, 0, 1, 1, 2, 1, 3, 1, ...
               0, 0, 1, 0, 2, 0, 2, 2, 3, ...
               2, 4, 2, 0, 0, 1, 0, 2, 0]; 
B_copy_tab = [ 0, 0, 0, 1, 0, 2, 0, 0, 0, ...
               1, 0, 2, 0, 0, 0, 1, 0, 2, ...
               1, 1, 1, 2, 1, 3, 0, 0, 0, ...
               1, 0, 2, 2, 2, 2, 3, 2, 4];
AB_copy_tab = A_copy_tab + x_dim_eq*B_copy_tab; % here put A in the lsb (base 5) and B in the  msb (another 5)

alpha_copy_tab = repmat([0 0 1 1 2 2], 1, 6)';
beta_copy_tab = [0 0 0 0 0 0 0 0 0 0 0 0 ...
                 1 1 1 1 1 1 1 1 1 1 1 1 ...
                 2 2 2 2 2 2 2 2 2 2 2 2]';
alpha_geno_tab = repmat(0:1, 1,18)'; 
beta_geno_tab = repmat([0 0 0 0 0 0 1 1 1 1 1 1], 1, 3)';

M = zeros(36,7);

M(:,1) = A_copy_tab; M(:,2) = B_copy_tab; M(:,3) = A_copy_tab + B_copy_tab;
M(:,4) = alpha_copy_tab;   M(:,5) = beta_copy_tab;
M(:,6) = alpha_geno_tab;   M(:,7) = beta_geno_tab;


A_copy_tab2 = alpha_copy_tab.*(1-alpha_geno_tab) + beta_copy_tab.*(1-beta_geno_tab);
B_copy_tab2 = alpha_copy_tab.*alpha_geno_tab + beta_copy_tab.*beta_geno_tab;

A_copy_tab2'-A_copy_tab
B_copy_tab2'-B_copy_tab
alpha_copy_tab'  + beta_copy_tab' - A_copy_tab - B_copy_tab

figure; imagesc(M); colorbar; title('A copy   B copy   Tot. copy  \alpha copy.  \beta copy  \alpha geno.  \beta geno');

