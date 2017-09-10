% Save data in PLINK format
% 
% Input: 
% snp_mat - matrix of SNPs values for each individual
% snp_ids - names of SNPs (rs ids) 
% people_ids - names of individuals
% allele_freq - risk-allele frequencies 
% PLINK_file - where to save output 
%
function save_genotypes_in_PLINK_format(snp_mat, snp_ids, people_ids, ...
    allele_freq, phenotype_vec, PLINK_file)

snp_type = 'binary'; 
[num_snps num_individuals] = size(snp_mat);

if(~exist('A_allele_vec', 'var') || isempty(A_allele_vec))
    A_allele_vec = repmat('A', num_snps, 1);
    B_allele_vec = repmat('B', num_snps, 1);
end

if(~exist('snp_ids', 'var') || isempty(snp_ids))
    snp_ids = mat2cell([repmat('rs', num_snps, 1) num2str((1:num_snps)')], ones(num_snps,1));
    snp_ids =  strrep_cell(snp_ids, ' ', '');     
end
if(~exist('chr_vec', 'var') || isempty(chr_vec))
    chr_vec = ones(num_snps, 1);
end
if(~exist('pos_vec', 'var') || isempty(pos_vec))
    pos_vec = 10 * (1:num_snps)';
end
if(~exist('morgans_vec', 'var') || isempty(morgans_vec))
    morgans_vec = repmat(0.01, num_snps, 1);
end


R = {}; % Save .MAP file (one line per each SNP)
% Format:  chromosome (1-22, X, Y or 0 if unplaced), rs# or snp identifier, Genetic distance (morgans), Base-pair position (bp units)

R = [num2cell(chr_vec) snp_ids num2cell(morgans_vec) ...
    num2cell(pos_vec) cellstr(A_allele_vec) cellstr(B_allele_vec)];
R = num2str_cell(R, 3);
switch snp_type
    case 'binary' % exclude half snps
        R = R(1:2:end,:);
end
savecellfile(R, [PLINK_file, '.bim']); % Save as .bim file 
savecellfile(R(:,1:4), [PLINK_file, '.map']); % .map contains only first 4 columns


R = {}; % Save .PED file (one line per each individual)
% Format:  Family-ID, Individual-ID, Paternal-ID, Maternal-ID, Sex (1=male; 2=female; other=unknown), Phenotype

if(~exist('people_ids', 'var') || isempty(people_ids))
    people_ids = (1:num_individuals)';
end

if(~exist('family_ids', 'var') || isempty(family_ids))
    family_ids = ones(num_individuals, 1);
end
if(~exist('sex', 'var') || isempty(sex))
    sex = ones(num_individuals, 1); % default: male 
end

if(~exist('paternal_ids', 'var') || isempty(paternal_ids))
    paternal_ids = zeros(num_individuals, 1);
    maternal_ids = zeros(num_individuals, 1);
end

R = [family_ids people_ids paternal_ids maternal_ids sex phenotype_vec]; % assume everything is 'double'
R = num2str_cell(num2cell(R), 3);
%R = {}; % Save .fam file (one line per each individual. These are the first 6 columns of .ped file)
savecellfile(R, [PLINK_file, '.fam']);
snp_mat_str = num2str(snp_mat');
remove_space_inds = 2+(0:(num_snps/2))*6;
keep_space_inds = setdiff(1:size(snp_mat_str,2), remove_space_inds);
snp_mat_str = mat2cellRows(snp_mat_str(:,keep_space_inds)); 

snp_mat_str = strrep_cell(snp_mat_str, '0', 'A');
snp_mat_str = strrep_cell(snp_mat_str, '1', 'B');
%snp_mat_str = strrep_cell(snp_mat_str, '  ', ' ')

R = [R snp_mat_str];
savecellfile(R, [PLINK_file, '.ped']);

savecellfile(R(:, [1 2 6]), [PLINK_file, '.phen']);% save only family-id, people-id, phenotype to .phen file




