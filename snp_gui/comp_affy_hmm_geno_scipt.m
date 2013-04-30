function comp_affy_hmm_geno_scipt()

affy_geno_path = 'E:\Libi\tools\SNP_tool\data\Leukemia_dchip_norm\';

hmm_geno_path = 'E:\Libi\tools\SNP_tool\data\Leukemia_dchip_norm\display\';
% Load Chip Annotations - We need the strand from here !!!
genome_assembly = get_genome_assembly();
load(['..\database\Hind_annot_data_' genome_assembly '.mat'], 'snp_ids', 'strand');

% Here we need to reverse the relevant SNPs according to strand
HindStrandSigns = zeros(1,length(strand));
HindStrandSigns(strmatch('-', strand)) = 1;


hind_info_file_table = loadCellFile([affy_geno_path 'Leukemia_samples_info_file_Hind.txt']);
first_line = hind_info_file_table(1, :);
samples_column = 1;

hind_samples = hind_info_file_table(2:end, samples_column);
num_samples = size(hind_samples, 1);
affy_geno_mat = zeros(0, num_samples);
hmm_geno_mat = zeros(0, num_samples);
for i = 1:num_samples
    sample = char(hind_samples{i});
    % load affy
    Affy = load([affy_geno_path sample '_hind.mat']);
    % load hmm
    HMM = load([hmm_geno_path sample '_hind_disp.mat']);
    t=1;
    affy_snp_ids = Affy.snp_ids;
    affy_geno = Affy.genotype_vec;

    ind = 1;
    for j = 1:length(HMM.DispStruct.Chrom)
        if(length(HMM.DispStruct.Chrom{j})>0)
            hmm_chr_snp_ids = HMM.DispStruct.Chrom{j}.SNPsIDs';
            hmm_chr_geno = HMM.DispStruct.Chrom{j}.Genotypes;
            [C, IA, IB] = intersect(hmm_chr_snp_ids, affy_snp_ids);
            [temp, IC, ID] = intersect_order_by_first_gr(C, snp_ids);
            if(length(temp) ~= length(C))
                b=1
            end
            HindStrandSigns_chr = HindStrandSigns(ID);
            num_chr_snps = length(IA);
            geno_vec = hmm_geno_into_old(hmm_chr_geno(IA));
            hmm_geno_mat(ind:ind+num_chr_snps-1, i) = hmm_geno_into_affy(geno_vec, HindStrandSigns_chr);
            affy_geno_mat(ind:ind+num_chr_snps-1, i) = affy_geno(IB);
            ind = ind+num_chr_snps;
        end
    end
    % plot figure
    num_diff_geno = length(find(hmm_geno_mat(:,i) ~= affy_geno_mat(:,i)));
    figure; plot(hmm_geno_mat(:,i), affy_geno_mat(:,i), '.');
    hold on; title([erase_backslash(sample) ':Num Diff.: ' num2str(num_diff_geno)]); ylabel('Affy genotype'); xlabel('HMM genotype');
end

