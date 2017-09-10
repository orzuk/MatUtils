%function dist_cell = gene_dist_vec_into_dist_cell(dist_vec)
function dist_cell = gene_dist_vec_into_dist_cell(dist_vec)

dist_cell = cell(size(dist_vec));

dist_cell(find(dist_vec==-3)) = {'CDS'};
dist_cell(find(dist_vec==-2)) = {'intron'};
dist_cell(find(dist_vec==-1)) = {'UTR'};

non_gene_ind = find(dist_vec>0);
for i = 1:length(non_gene_ind)
    dist_cell{non_gene_ind(i)} = num2str(dist_vec(non_gene_ind(i)));
end