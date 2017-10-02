% Get metadata: file names, directories etc.
% Input:
% exome_data - string representing exome datasets
% Output:
% exome_struct - structure with meta-data for this dataset
%
function exome_struct = get_exome_data_info(exome_data)

exome_data_str = {'Nelson', 'Tennensen', 'TennensenSmall', 'ESP', 'ExAC'}; % name of different datasets
num_datasets = length(exome_data_str);
exome_ind = strfind_cell(exome_data_str, exome_data);

spectrum_data_files = {'/Nelson_Science_2012/Nelson_Science_data_table_S2.txt', ...
    '/Tennessen_Science_2012/all_chr_ESP6500.snps.vcf', ...
    '/Tennessen_Science_2012/all_chr_ESP6500.snps.small_gene_list.vcf', ...
    '/ESP/ESP6500SI*.snps_indels.vcf', ...
    '/ExAC/ExAC.chunks001.vcf', ...  %  '/ExAC/ExAC.Small.sites.vep.vcf', ...  %    '/ExAC/exome_aggregation.snps_indels.vcf', ... % NEW! exome aggregation data (much richer!!!)
    []}; % NEW: add gnomad dataset  % Data file with exome sequencing results %  '/Tennessen_Science_2012/ESP6500.chr19.snps.vcf', []}; % S3

spectrum_files_prefix = {'Nelson', 'all_chr_ESP6500', 'all_chr_ESP6500.small_gene_list', ...
    'ESP6500', 'ExAC.chunks'}; % 'ESP6500*.chr', 'ExAC'};


populations_vec = cell(num_datasets, 1); pop_str_vec = cell(num_datasets, 1);
populations_vec{1} = {'European', 'African'}; pop_str_vec{1} = {'EA_AC', 'AA_AC'};
populations_vec{2} = {'European', 'African'}; pop_str_vec{2} = {'EA_AC', 'AA_AC'};
populations_vec{3} = {'European', 'African'}; pop_str_vec{3} = {'EA_AC', 'AA_AC'};
populations_vec{4} = {'European', 'African'}; pop_str_vec{4} = {'EA_AC', 'AA_AC'}; % ESP
populations_vec{5} = {'European', 'African', 'Finish', 'SouthAsian', 'EastAsian', 'LatinoAmerican', 'Other'};
pop_str_vec{5} = {'AC_NFE', 'AC_AFR',  'AC_FIN', 'AC_SAS', 'AC_EAS', 'AC_AMR', 'AC_OTH'}; % ExAC


spectrum_labels = {'Nelson et al.', 'Tennessen et al.', 'Tennessen et al. few genes', 'Keinan et al.', 'Exome Aggregation'};
target_length_vec = [ 50*10^6/100,  50*10^6, 50*10^6/100, 50*10^6, 50*10^6]; % estimated # of nucleotides in target

exome_struct.data_str = exome_data; % copy name
exome_struct.spectrum_data_file = spectrum_data_files{exome_ind};
exome_struct.sub_dir_str = dir_from_file_name(exome_struct.spectrum_data_file);
exome_struct.spectrum_label = spectrum_labels{exome_ind};
exome_struct.target_length = target_length_vec(exome_ind);
exome_struct.populations = populations_vec{exome_ind};
exome_struct.pop_str = pop_str_vec{exome_ind};
exome_struct.prefix = spectrum_files_prefix{exome_ind};
exome_struct.triplet_mutations_file = 'scone_hg17_m3_64x64_rates.txt'; % 'triplets_file_human_chimp_baboon.mat'; % file with mutation rates from all 64 codons (estimated from human-chimp-baboon alignment)
exome_struct.mutation_rates_file = 'human_genes_mutation_rates_hg18.mat'; % output file with triplet mutation rates
exome_struct.exons_file = 'hg18_exons.mat'; % where to save exons

% New: set file with exome data (should fit ExAC dataset)
switch exome_data
    case 'ExAC'
        exome_struct.spectrum_data_files_str = fullfile(exome_data, 'AllPop', 'ExAC.chunks00*_AllPop.mat'); % ...
%            ['{''' fullfile(exome_data, 'AllPop', 'ExAC.chunks00*_AllPop.mat') '''}'];  % File with exome data (but only one chunk!!)
    otherwise
        exome_struct.spectrum_data_files_str = '{';  % File with exome data
end
