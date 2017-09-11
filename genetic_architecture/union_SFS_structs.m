% Unite to allele freq. structures (from different chromosomes, genome chunks etc.)
function A = union_SFS_structs(A1, A2, unite_field_names)

if(~exist('unite_field_names', 'var') || isempty(unite_field_names))
    unite_field_names = {'XXX_VARIANT_COUNT_', 'XXX_REF_ALLELE_COUNT_', 'XXX_FEATURE_', 'GENE', 'XXX_CHROM', ...
        'POS',  'ALLELE_FREQ',   'GENE_INDS', 'unique_genes'};
end
unite_field_names = intersect(fieldnames(A2), unite_field_names);
A=A1;
for j=1:length(unite_field_names)  % concatenate all chromosomes to one file (is this for population file? or one file for all populations?)
    unite_str = ['A.' unite_field_names{j} ' = [A.' unite_field_names{j} ''' A2.' unite_field_names{j} ''']'';'];
    eval(unite_str)
end
[~, I_types, J_types] = intersect(A.allele_types, A2.allele_types); % Here unite cell-array. Problem! for different chromosomes might have different #allele_types and their encoding
for j=1:length(I_types)
    A.n_vec{I_types(j)} = [A.n_vec{I_types(j)}' A2.n_vec{J_types(j)}']';
    A.f_vec{I_types(j)} = [A.f_vec{I_types(j)}' A2.f_vec{J_types(j)}']';
    A.count_vec{I_types(j)} = [A.count_vec{I_types(j)}' A2.count_vec{J_types(j)}']';
end
[~, I_diff_types] = setdiff(A2.allele_types, A.allele_types); % Here get different #allele_types and their encoding
for j=1:length(I_diff_types)% New alleles
    A.n_vec{A.num_allele_types+j} = A2.n_vec{I_diff_types(j)};
    A.f_vec{A.num_allele_types+j} = A2.f_vec{I_diff_types(j)};
    A.count_vec{A.num_allele_types+j} = A2.count_vec{I_diff_types(j)};
    A.allele_types{A.num_allele_types+j} = A2.allele_types{I_diff_types(j)}; % add allele types
end
A.num_allele_types = length(A.allele_types); % update # of allele types
