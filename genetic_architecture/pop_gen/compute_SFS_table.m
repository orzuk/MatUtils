%
% Compute different SFS statistics
% Input: 
% A - structure with SFS data
% MutationRateTable - mutation rate for each gene and each class
% MutationTypes - 
% output_file_name - 
% Output:
% R - cell array with alleles information 
% 
function R = compute_SFS_table(A, MutationRateTable, MutationTypes, exome_struct, output_file_name)

Assign24MammalsGlobalConstants;

TotalMutationTargetSize = sum(MutationRateTable);
num_classes = length(MutationTypes);
num_populations = 2; % tmp!! this shouldn't be hard-coded but given as input!!!
output_file_dir = dir_from_file_name(output_file_name);
R = cell(20,5); % cell type
R{1,1} = ['Summary' exome_struct.data_str ' Statistics for All Genes'];
R{2,1} = 'Num Genes:'; R{2,2} = num2str(A.num_genes);
ctr=3;

R{ctr,1} = ['Mutation Rates and Target Size for Different Classes (' exome_struct.data_str ' Data):']; ctr=ctr+1;
R{ctr,1} = '----------------------------------------------------------------'; ctr=ctr+1;
R{ctr,1} = 'Class:';
R{ctr+1,1} = 'Target Size (\mu):';
R{ctr+2,1} = 'Target Size (%):';
for i=1:num_classes
    R{ctr, 1+i*num_populations} = genome_types{MutationTypes(i)}; % string with mutation type
    R{ctr+1, 1+i*num_populations} = TotalMutationTargetSize(i); % total genomic target size
    R{ctr+2, 1+i*num_populations} = 100 * TotalMutationTargetSize(i) ./ ...
        sum(TotalMutationTargetSize); % percent genomic target size
end
ctr = ctr+5; % jump to next

R{ctr,1} = '#Alleles:'; % start filling
R{ctr+1,1} = 'Mean DAF'; % derived allele freq.
for j=1:length(MutationTypes) % loop on diffferent allele types
    for i=1:length(A.upper_freq_vec) % loop on thresholds
        cur_mutation_type_ind = A.good_allele_inds{4}(j); % index in list of alleles  (here take sub-classes)
        cur_mutation_type_rate_ind = ... % index in list of genomic mutation types
            find(A.allele_types_ind(cur_mutation_type_ind) == MutationTypes); % find index representing mutation
        R{ctr+5+j,1} = ['Cumulative freq. % (<' num2str(100*A.upper_freq_vec(j), 3) '%):'];
        R{ctr+5+j,num_populations*i+2} =  100*sum(A.total_freq_per_gene_mat{j}(cur_mutation_type_ind, :));
        R{ctr+5+length(A.upper_freq_vec)+j,1} = ['Heterozygisity %(<' num2str(100*A.upper_freq_vec(j), 3) '%):'];
        R{ctr+5+length(A.upper_freq_vec)+j,num_populations*i+2} =  100*sum(A.total_heterozygosity_per_gene_mat{j}(cur_mutation_type_ind, :));
    end
end

my_mkdir(output_file_dir); 
savecellfile(R, fullfile(output_file_dir, [exome_struct.data_str '_exome_statistics.txt']));
