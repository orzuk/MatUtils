%function ind = find_snp_gene_symbols_ind(snp_gene_symbols, gene_symbols)
function ind = find_snp_gene_symbols_ind(snp_gene_symbols, gene_symbols)

[unique_snp_gene_symbols,I,J] = unique(snp_gene_symbols); 
%    that B = A(I) and A = B(J) (or B = A(I,:) and A = B(J,:)).
    
[c, ai, bi] = intersect(unique_snp_gene_symbols, gene_symbols);

ind = zeros(length(unique_snp_gene_symbols), 1);
ind(ai) = bi;
ind = ind(J);
