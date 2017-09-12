function make_segal_genesets_DB

load ORF_gene_symbol.mat

% genesets_mat=cell(0,6);
% files=dir;
% for i=1:length(files)
%     if length(files(i).name)>4 && strcmpi(files(i).name(end-3:end),'.gxa')
%         disp([num2str(i) ': ' files(i).name]);
%         genesets_mat=[genesets_mat; read_gxa_file(files(i).name)];
%     end   
% end
% 
% save tmp_segal_gene_sets.mat genesets_mat


union_gene_symbols={};
ORF_double=cell2mat(ORF_gene_symbol(:,1));
for i=1:size(genesets_mat,1)
    gene_symbols={};
    for j=1:length(genesets_mat{i,6})
        idx=find(ORF_double==genesets_mat{i,6}(j));
        for k=1:length(idx)
            if ischar(ORF_gene_symbol{idx(k),2})
                gene_symbols=[gene_symbols; ORF_gene_symbol{idx(k),2}];
            end
        end
    end
    
    gene_symbols=unique(gene_symbols);
    genesets_mat{i,5}=length(gene_symbols);
    genesets_mat{i,6}=gene_symbols;
    union_gene_symbols=union(union_gene_symbols,gene_symbols);
    if mod(i,100)==0
        disp(i);
    end
end


genesets_vs_symbols=sparse(size(genesets_mat,1),length(union_gene_symbols));
for i=1:size(genesets_vs_symbols,1)
    [dum idx]=ismember(genesets_mat{i,6},union_gene_symbols);
    genesets_vs_symbols(i,idx)=1;
    genesets_mat{i,6}=idx;
end

save segal_gene_sets.mat genesets_mat union_gene_symbols genesets_vs_symbols





