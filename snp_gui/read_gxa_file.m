function genesets_mat=read_gxa_file(gxa_file)

theStruct=parseXML(gxa_file);

if ~strcmp(theStruct.Children(2).Name,'GeneXPressAttributes') || ~strcmp(theStruct.Children(4).Name,'GeneXPressObjects')
    error(['Error in format in ' gxa_file]);
end

if ~strcmp(theStruct.Children(2).Children(2).Attributes(4).Name,'Name')
    error(['Error in format in ' gxa_file]);
end
geneset_type=theStruct.Children(2).Children(2).Attributes(4).Value;

if ~strcmp(theStruct.Children(2).Children(2).Attributes(5).Name,'NumAnnotations')
    error(['Error in format in ' gxa_file]);
end
num_of_genesets=str2num(theStruct.Children(2).Children(2).Attributes(5).Value);

if ~strcmp(theStruct.Children(2).Children(2).Attributes(6).Name,'Organism')
    error(['Error in format in ' gxa_file]);
end
geneset_organism=theStruct.Children(2).Children(2).Attributes(6).Value;

num_of_genes=(length(theStruct.Children(4).Children(2).Children)-1)/2;
ORF_genes=zeros(num_of_genes,1);
genes_exist_in_geneset=zeros(num_of_genes,num_of_genesets);

for i=1:num_of_genes
    if ~strcmp(theStruct.Children(4).Children(2).Children(i*2).Name,'Gene')
        error(['Error in format in ' gxa_file]);
    end
    
    if ~strcmp(theStruct.Children(4).Children(2).Children(i*2).Attributes(2).Name,'ORF')
        error(['Error in format in ' gxa_file]);
    end
    
    ORF_num=str2num(theStruct.Children(4).Children(2).Children(i*2).Attributes(2).Value);
    if isempty(ORF_num)
        ORF_num=0;
    end
    ORF_genes(i)=ORF_num(1);
    
    eval(['exist_in_geneset=[ ' theStruct.Children(4).Children(2).Children(i*2).Children(2).Attributes(3).Value '];']);
    if length(exist_in_geneset)~=num_of_genesets
        error(['Error in format in ' gxa_file]);
    end
    genes_exist_in_geneset(i,:)=exist_in_geneset;
end




genesets_mat=cell(num_of_genesets, 6);
for i=1:num_of_genesets
    if str2num(theStruct.Children(2).Children(2).Children(2*i).Attributes(2).Value)<i
        error(['Error in format in ' gxa_file]);
    end
    if ~strcmp(theStruct.Children(2).Children(2).Children(2*i).Attributes(3).Name,'Name')
        error(['Error in format in ' gxa_file]);
    end
    genesets_mat{i,1}=theStruct.Children(2).Children(2).Children(2*i).Attributes(3).Value;
    genesets_mat{i,2}=geneset_type;
    genesets_mat{i,3}=genesets_mat{i,1};
    genesets_mat{i,4}=geneset_organism;
    genesets_mat{i,6}=ORF_genes(find(genes_exist_in_geneset(:,i)));
    genesets_mat{i,5}=length(genesets_mat{i,6});
    
end


        


