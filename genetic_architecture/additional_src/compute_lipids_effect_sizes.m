% A script for computing Lipids effect sizes

% Parse file
L = ReadDataFile('../../common_disease_model/data/lipids/lipids_revised_with_maf.txt'); % read table from Eliana
%load('Y:\public_html\data\common_disease_broad_catalog\disease_data.mat'); % read current data
load('../../common_disease_model/data/gwas_mit_turk_catalog/GWAS_catalog_parameters_2011_03_01_processed.mat');
L_chr = ReadDataFile('../../common_disease_model/data/lipids/110302223217_CHR.txt',[],0);

L_chr_inds = strmatch('GRCh37', L_chr.version); % CRA_TCAGchr7v2, GRCh37, Celera HuRef
L_chr = struct_by_inds(L_chr, L_chr_inds);
L_chr.rs_ = num2str_cell(L_chr.rs_);
for i=1:length(L_chr.rs_)
    L_chr.rs_{i} = ['rs' L_chr.rs_{i}];
end

%manual_tab =

num_snps = length(L.MarkerName)
L.Gene = cell(num_snps,1); L.chr = zeros(num_snps,1); L.positions = zeros(num_snps,1); 
good_inds = strmatch('Lipid', data.Trait)
data.SNP_ID_ = strrep_cell(data.SNP_ID_, ' ', '');
[L.chr L.positions] = rs_ids_to_genomic_coordinates(L.MarkerName, 'hg18'); % add position information
for trait = {'LDL', 'HDL', 'TG', 'TC'} % loop on 3 lipids
    trait_inds = strmatch(trait, L.Trait)
    [~, I, J] = intersect(L.MarkerName(trait_inds), data.SNP_ID_) % (good_inds))
    L.Gene(trait_inds(I)) = data.Gene(J); % (good_inds(J));
    L.chr(trait_inds(I)) = data.chr(J);
    
    [~, I, J] = intersect(L.MarkerName(trait_inds), L_chr.rs_); % Complete positions
    L.positions(trait_inds(I)) = cell2mat(empty_cell_to_numeric_val(str2num_cell(L_chr.pos(J)),0));
    L.chr(trait_inds(I)) = cell2mat(empty_cell_to_numeric_val(str2num_cell(L_chr.chr(J)),0));
end


WriteDataFile(L, '../../common_disease_model/data/lipids/lipids_revised_with_maf_processed.txt');










