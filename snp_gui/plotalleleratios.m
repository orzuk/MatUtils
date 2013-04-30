% Plot two allele copy numbers
% We show the intensity levels of two alleles along with their jointly
% determined genotype calls
function PlotAlleleRatios(copy_num_vec, allele_ratio_vec, genotype_vec, title_str, new_figure_flag)

[copy_A, copy_B] = RatioToCopyMats(copy_num_vec, allele_ratio_vec); 

if(nargin == 5)
    PlotABCopyNums(copy_A, copy_B, genotype_vec, title_str, new_figure_flag);
else
    PlotABCopyNums(copy_A, copy_B, genotype_vec, title_str);
end    
    
