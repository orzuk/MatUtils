A package for False-Discovery Rate estimation.

This package accompanies the paper: 
FDR Control with adaptive procedures and FDR monotonicity (2009), A. Zeisel, O. Zuk and E. Domany

List of files:


Installation: 

Simply download the file fdr.tgz, unzip to a directory of your choice and add to your Matlab path.

The file FDR_mat_main is used as the main FDR procedure, given a set of p-values

The function compute_IBH_constants is used to compute the C and S constants used in the IBH procedure.

Gene expression data p-values are available in gene_expression_pvals.mat.

A counter example showing that the FDR may not be monotonic is given in: FDR_monotonic_counterexample.m

The function simulate_pvals.m generates several p-values distributions (including dependencies, null and non-null p-values etc.) 

The figures produced in the paper were done using the scripts: 
Pr_VoverR_leq_B_maps.m
Pr_VoverR_leq_B_section.m 
