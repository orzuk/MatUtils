% Display a PCA of the samples.
% One can choose either to look at genotypes or on copy number
% 
% The inputs are : 
% ALL_MAT - data matrix (including copy # and genotypes for all chromosomes)
% data_flag - 0 is by copy number, 1 - is by genotypes
% chrom - show a PCA of a specific chromosome. -1 is to use all chromomosmes.
% Labels - The label for each sample. The coloring is done according to these labels
function DisplaySamplesPCA(ALL_MAT, data_flag, data_labels, chrom, NumRandSNPs)
num_snps = length(ALL_MAT.snp_ids);

%% NumRandSNPs = 1000; % Choose some random SNPs. Default is 1000.
%% RandSNPs = 1:NumRandSNPs; % Choose the first ..
if(chrom == -1) % Pick randomly from the whole genome
    RandSNPs = randperm(num_snps); RandSNPs = RandSNPs(1:NumRandSNPs);  % Choose random
else % Pick a specific chromosome
    [chr_starts chr_ends] = get_chrs_limits(1:24);
    chr_ends = cumsum(chr_ends);
    [ALL_MAT.chr_vec, ALL_MAT.locs, ALL_MAT.pres_ind] = get_chrom_locs(ALL_MAT.snp_ids, ALL_MAT.chip_type);
    chrom_inds = find(ALL_MAT.chr_vec == chrom);
    NumRandSNPs = min(NumRandSNPs, length(chrom_inds)); % Could be too few SNPs on the chrom. 
    RandSNPs = randperm(length(chrom_inds)); RandSNPs = chrom_inds(RandSNPs(1:NumRandSNPs));  % Choose random SNPs from this chromosome
    
    % Old way - not good and assumes that SNPs are already ordered
    %     start_snp = find(ALL_MAT.chr_vec == chrom, 1 );
    %     end_snp = find(ALL_MAT.chr_vec == chrom, 1, 'last' );
    %     NumRandSNPs = min(NumRandSNPs, end_snp-start_snp+1); % Could be too few SNPs on the chrom.
    %     RandSNPs = randperm(end_snp-start_snp+1); RandSNPs = RandSNPs(1:NumRandSNPs) + start_snp-1;  % Choose random
end

% Check which labels share the same colors
[legend_labels I Colors] = unique(data_labels);

if(data_flag == 0) % Show copy number
    [V Y] = PachterPolytopeAnalysis( (ALL_MAT.data_A(:,RandSNPs) + ALL_MAT.data_B(:,RandSNPs))' , data_flag); %%% SKIP this for now ....
    data_flag_str = 'Copy Number';
else  % show genotypes
    [V Y] = PachterPolytopeAnalysis(ALL_MAT.data_genotype_AB(:,RandSNPs)', data_flag); %%% SKIP this for now ....
    data_flag_str = 'Genotypes';
end

% Show the plot according to the labels ..
figure; hold on; 
% Prepare matrix
PlotX = zeros(length(legend_labels), length(data_labels));
PlotY = zeros(length(legend_labels), length(data_labels));
PlotX(:) = nan; PlotY(:) = nan;

for i=1:length(legend_labels)
    PlotX(i, Colors == i) = Y(1,Colors == i);
    PlotY(i, Colors == i) = Y(2,Colors == i);    
end

plot(PlotX', PlotY', '*', 'LineWidth',5);
title(['First Components Projection of ' data_flag_str]); 
legend(legend_labels); % Need to solve this legend problem ..

%% figure;  scatter(Y(1,:), Y(2,:), 5, Colors);  % We replaced scatter by a loop with plots.

% HAPMAP Sample Names
% SampleNames = {'NA06985_Hind_B5_3005533affybench', ...
% 'NA07029_Hind_A9_4000092affybench',	...
% 'NA07357_Hind_B10_3005533affybench', ...
% 'NA06991_Hind_B6_3005533affybench',	...
% 'NA07034_Hind_B1_4000092affybench',	...
% 'NA10830_Hind_E6_4000092affybench',	...
% 'NA06993_Hind_B4_4000092affybench',	...
% 'NA07048_Hind_B3_4000092affybench',	...
% 'NA10831_Hind_E9_4000092affybench',	...
% 'NA06994_Hind_A7_3005533affybench',	...
% 'NA07055_Hind_B2_4000092affybench',	...
% 'NA10835_Hind_E12_4000092affybench', ...
% 'NA07000_Hind_A8_3005533affybench',	...
% 'NA07056_Hind_A11_4000092affybench', ...	
% 'NA07019_Hind_A12_4000092affybench', ...	
% 'NA07345_Hind_B11_3005533affybench', ...
% 'NA07022_Hind_A10_4000092affybench', ...
% 'NA07348_Hind_B12_3005533affybench'};
% 

% Leukemia Sample Names
SampleNamesLeukemaia = {'TEL74_5_d', 'TEL74_5_n',	'TEL76_7_d', 'TEL76_7_n',	...
'HD78_9_d',	'HD78_9_n',	'TEL80_1_d', 'TEL80_1_n',	...
'HD82_3_d',	'HD82_3_n',	'TEL84_5_d', 'TEL84_5_n',	...
'HD86_7_d',	'HD86_7_n',	'TEL88_9_d', 'TEL88_9_n',	...
'HD90_1_d',	'HD90_1_n',	'HD92_3_d',	'HD92_3_n',	...
'TEL94_5_d', 'TEL94_5_n',	'TEL96_7_d', 'TEL96_7_n'};	
% 
