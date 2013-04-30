function combine_refgene_reflink(refgene_file,reflink_file,genome_build)

raw_DB=fullfile('..','raw_database');

refgene_mat=loadCellFile(fullfile(raw_DB,refgene_file)); %take only 11 columns
reflink_mat=loadCellFile(fullfile(raw_DB,reflink_file));

[dum idx]=ismember(refgene_mat(:,2),reflink_mat(:,3));

idx=idx(idx>0);

gene_table=cell(length(idx)+1,13);
gene_table(1,:)={'name' 'chrom' 'strand' 'txStart' 'txEnd' 'cdsStart' 'cdsEnd' 'exonCount' 'exonStarts' 'exonEnds' 'name' 'product' 'name'};
gene_table(2:end,1:10)=refgene_mat(dum,2:11);
gene_table(2:end,11:13)=reflink_mat(idx,1:3);

saveCellFile(gene_table,fullfile(raw_DB,['refgenes_' genome_build '.txt']));

save(fullfile(raw_DB,['refgenes_' genome_build '_cell.mat']),'gene_table');


