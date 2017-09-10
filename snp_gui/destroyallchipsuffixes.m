% A short utility that tries to eliminate all the '_hind' created by older
% version of the program WITHOUT running everything again
function SpecialInds = DestroyAllChipSuffixes(cur_dir)

%% cur_dir = 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\CEU';


FileNames = dir(fullfile(cur_dir, '*xba.mat')); % also hind needed
% FileNamesHind = dir(fullfile(cur_dir, '*hind.mat')); % also hind needed


SpecialInds =[];

for i=1:length(FileNames)
    clear snp_id_hind copy_num_vec_hind allele_ratio_vec_hind genotype_vec_hind snp_id_xba copy_num_vec_xba allele_ratio_vec_xba genotype_vec_xba;

    i
    load(fullfile(cur_dir,FileNames(i).name));

    if(exist('snp_id_hind'))
        snp_ids = snp_id_hind;
        if(exist('copy_num_vec_hind') & exist('allele_ratio_vec_hind'))

            copy_num_vec = copy_num_vec_hind;
            allele_ratio_vec = allele_ratio_vec_hind;
            chip = 'Hind';
            if(exist('genotype_vec_hind'))
                genotype_vec = genotype_vec_hind;
                save(fullfile(cur_dir,FileNames(i).name), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec', 'genotype_vec', 'chip');
            else
                save(fullfile(cur_dir,FileNames(i).name), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec', 'chip');
            end
        else % What the fuck do we do here ?
            SpecialInds = [SpecialInds i]
            save(fullfile(cur_dir,FileNames(i).name), 'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB', ...
                'SampleNames', 'chip', 'snp_ids');
        end


    end


    if(exist('snp_id_xba'))
        if(exist('copy_num_vec_xba') & exist('allele_ratio_vec_xba'))
            snp_ids = snp_id_xba;
            copy_num_vec = copy_num_vec_xba;
            allele_ratio_vec = allele_ratio_vec_xba;
            chip = 'Xba';
            if(exist('genotype_vec_hind'))
                genotype_vec = genotype_vec_xba;
                save(fullfile(cur_dir,FileNames(i).name), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec', 'genotype_vec', 'chip');
            else
                save(fullfile(cur_dir,FileNames(i).name), 'snp_ids', 'copy_num_vec', 'allele_ratio_vec', 'chip');
            end
        else
            SpecialInds = [SpecialInds i]
        end

    end


end