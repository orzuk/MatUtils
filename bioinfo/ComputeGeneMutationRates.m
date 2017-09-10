% Calculate total mutation rates for different types of coding mutations for the human genome
% 
% Input: 
% TripletsMutationTable - a 64x64 table with mutation rates from each triplet to each triplet 
% GeneStruct - structure contatining all genes (names, positions, sequences etc.)
%
% Output: 
% MutationRateTable - a table with total mutation rates for each gene. Size: (#genes X #mutation-types)
% MutationTypes - indices with different types of mutations 
% MutationTypesStr - strings with different types of mutations 
% 
function [MutationRateTable MutationTypes MutationTypesStr] = ComputeGeneMutationRates(TripletsMutationTable, GeneStruct)



Assign24MammalsGlobalConstants;
[mutation_types_table codons] = get_mutation_types(); % get all possible mutation types 

num_mutation_types = 4; % synonymous, missesne, stop-gained, stop-lost

% MutationRateTable = zeros(num_genes, num_mutation_types); 

TripletsMutationTable = TripletsMutationTable - diag(diag(TripletsMutationTable)); % get rid of self mutations
CodonMutationRateTable = zeros(64, num_mutation_types); % table of codon-by substitution type
MutationTypes = [SYNONYMOUS MISSENSE STOP_GAINED STOP_LOST]; % get labels (mutation types) 
MutationTypesStr = genome_types(MutationTypes); % get labels (mutation types) 
for i=1:64 % loop on possible codons 
    ctr=1;
    for j = MutationTypes % loop on different muation types
        CodonMutationRateTable(i, ctr) = sum( TripletsMutationTable(i,:) .* (mutation_types_table(i,:) == j) ); 
        ctr=ctr+1; 
    end
end


if(ischar(GeneStruct)) % Get sequences and compute total mutation rates by class
    GeneStruct = load(GeneStruct, ...
        'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm');
end
num_genes = length(GeneStruct.gene_names); 
codon_count = zeros(num_genes, 64);
for i=1:num_genes % Scan all possible genes
    if(mod(i,50) == 0)
        run_gene = i
    end
    [~, tmp_codon_count] = codoncount(GeneStruct.seqs{i});
    codon_count(i,:) = reshape(tmp_codon_count, 1, 64); 
end
MutationRateTable = codon_count * CodonMutationRateTable; 




