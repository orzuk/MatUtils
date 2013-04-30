% This function converts the output of dchip into a mat file for each sample in the
% hapmap hind samples. The files we produce will be used as input for the hmm program.
function HapMapDChipOutputTxtIntoMatFiles(probe_intens_path)

AssignAllGlobalConstants(); 

% Read the .txt file, extract all information and delete it's variable
AllFileStuff = loadcellfile(fullfile(probe_intens_path, 'hapmap_hind_raw_copy_num_dchip.txt'));
AllSnpsIDs = AllFileStuff(3:end,1);
AllRSIDs = AllFileStuff(3:end,2);
AllChrVec = char(AllFileStuff(3:end,3));
AllChrVec(find(AllChrVec(:,1) == 'X'),1) = '2'; AllChrVec(find(AllChrVec(:,1) == 'X'),2) = '3'; % chrom 'X' is 23
AllChrVec(find(AllChrVec(:,1) == 'Y'),1) = '2'; AllChrVec(find(AllChrVec(:,1) == 'Y'),2) = '4';
AllChrVec = str2num(AllChrVec);
AllCopyNum = str2num(char(AllFileStuff(3:end,7:end-1)));
AllCopyNum = reshape(AllCopyNum, length(AllCopyNum)/18, 18);
clear AllFileStuff; % Get rid - heavy memory ..


save 'E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\hapmap_hind_raw_copy_num_dchip.mat' ...
    'AllSnpsIDs' 'AllRSIDs' 'AllChrVec' 'AllCopyNum';

load('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\hapmap_hind_raw_copy_num_dchip.mat');
hind_txt_file_names = loadcellfile(fullfile(probe_intens_path, 'hind_txt_file_names.txt'));

for i=1:length(hind_txt_file_names)   % 2nd loop: Normalize and transfer to application input format
    doing_i = i

    % Load everything except the copy number
    eval(['load ' probe_intens_path hind_txt_file_names{i}(1:end-4) 'affybench_hind.mat;']); % load the allele_ratio_vec

    % Save back with the previous WRONG copy number - Do this only once !!!
     eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) ...
         'affybench_hind_bad_normalization2.mat snp_ids genotype_vec copy_num_vec allele_ratio_vec;']);

    % Calculate and update the copy number
    [snps_intersect I J] = intersect(snp_ids, AllSnpsIDs);
    snp_ids = snp_ids(I); % Now make everything smaller
    genotype_vec = genotype_vec(I);
    allele_ratio_vec = allele_ratio_vec(I);
    copy_num_vec = AllCopyNum(J,i);   % before it was  copy_num_vec(I) = AllCopyNum(J,i);

    
%%    [snps_minus I_minus] = setdiff(snp_ids, AllSnpsIDs); % This is just to see what snps are missing ...

    % Save back with the copy number
    eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) ...
        'affybench_hind.mat snp_ids genotype_vec copy_num_vec allele_ratio_vec;']);
end




