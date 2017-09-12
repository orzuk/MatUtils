%function snp_gene_descr = snp_gene_symbols_into_descr(snp_gene_symbols, gene_symbols, gene_descr)
function snp_gene_descr = snp_gene_symbols_into_descr(snp_gene_symbols, gene_symbols, gene_descr)

ind = find_snp_gene_symbols_ind(snp_gene_symbols, gene_symbols);
snp_gene_descr = cell(size(snp_gene_symbols)); snp_gene_descr(:) = {''};
snp_gene_descr(find(ind~=0)) = gene_descr(ind(find(ind~=0)));
