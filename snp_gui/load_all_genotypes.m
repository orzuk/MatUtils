%function [genotypes, snp_ids] = load_all_genotypes(samples, work_dir, chip)
function [genotypes, snp_ids] = load_all_genotypes(samples, work_dir, chip)

num_samples = length(samples);
genotypes = zeros(0,num_samples);
snp_ids = cell(0,1);
for i = 1:num_samples
    i
    sample = samples{i};
    load(fullfile(work_dir, 'display', [sample '_' lower(chip) '_disp.mat']));
    snp_ind=1;
    for j = 1:23
        num_snps_chr = length(DispStruct.Chrom{j}.Genotypes);
        genotypes(snp_ind:snp_ind+num_snps_chr-1,i) = DispStruct.Chrom{j}.Genotypes;
        if(i==1)
            snp_ids(snp_ind:snp_ind+num_snps_chr-1) = DispStruct.Chrom{j}.SNPsIDs';
        end
        snp_ind = snp_ind+num_snps_chr;
    end
end

