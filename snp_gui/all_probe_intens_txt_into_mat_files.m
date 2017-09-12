% This function reads several files of the hapmap hind samples and convert
% them into .mat files.
% The input files format is a txt file which is equivalent to the .cell files.
% The output file is in the format to be used by our HMM
function all_probe_intens_txt_into_mat_files(probe_intens_path, mat_f_name_path)

AssignAllGlobalConstants(); 

FIRST_TIME = 0;
if(FIRST_TIME) % Do very heavy stuff just once (reading .txt files)
    snp_ids = loadcellfile([ probe_intens_path 'Hind_SNP_IDs.txt']);
    eval(['save ' probe_intens_path 'Hind_SNP_IDs.mat snp_ids;']);

    genotypes_hind = loadcellfile([probe_intens_path 'Hind_Genotypes.txt']);
    eval(['save ' probe_intens_path 'Hind_Genotypes.mat genotypes_hind;']);

    %% genotypes_hind_only = genotypes_hind(2:end,[7:2:end-1]);
    %% genotypes_hind_num = zeros(size(genotypes_hind_only, 1), size(genotypes_hind_only, 2));
    %% genotypes_hind_num( strmatch('AA', genotypes_hind_only) ) = AA; % Set AA  (1)
    %% genotypes_hind_num( strmatch('AB', genotypes_hind_only) ) = AB; % Set AB  (2)
    %% genotypes_hind_num( strmatch('BB', genotypes_hind_only) ) = BB; % Set BB  (3)
    %% genotypes_hind_num( strmatch('NoCall', genotypes_hind_only) ) = NoCall; % Set NoCall (4)
    %% eval(['save ' probe_intens_path 'Hind_Genotypes_Only.mat genotypes_hind_only;']);
    %% eval(['save ' probe_intens_path 'Hind_Genotypes_Num.mat genotypes_hind_num;']);
else
    eval(['load(''' probe_intens_path 'Hind_SNP_IDs.mat'');']);
    RRR = 999;
end
eval(['load(''' probe_intens_path 'Hind_Genotypes_Num.mat'');']);

hind_txt_file_names = loadcellfile([probe_intens_path 'hind_txt_file_names.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert .txt to .mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CONVERT_TXT_TO_MAT_FILES = 0; % If 1 it takes a long time, since we convert all .txt files to .mat files!!!  (half a day !!!)
if(CONVERT_TXT_TO_MAT_FILES)
    median_vec = zeros(1,length(hind_txt_file_names));
    % Convert all text files to mat files - Very time consuming !!! (Can skip this)
    for i=1:length(hind_txt_file_names)
        txt_to_mat_i = i
        tmp_name = hind_txt_file_names{i}(1:end-4)
        median_vec(i) = probe_intens_txt_into_mat(probe_intens_path, hind_txt_file_names{i}, mat_f_name_path, ...
            [tmp_name '.mat'], tmp_name);
    end
    % Perform normalization:
    baseline_med = median(median_vec);
    [baseline_dummy baseline_ind] = min( (median_vec - baseline_med).^2 );
else
    baseline_ind = 15; % Set a fixed ind !
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start Normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[baseline_pm_column, pa_vec, ma_vec, pb_vec, mb_vec] = ...
    get_probe_pm_mm_columns(probe_intens_path, ...
    [hind_txt_file_names{baseline_ind}(1:end-4) '.mat'], hind_txt_file_names{baseline_ind}(1:end-4));

for i=1:length(hind_txt_file_names) % 1st loop we don't know yet normalization
    first_loop_i = i
    [cur_pm_column, probe_intens_mat, pa_vec, ma_vec, pb_vec, mb_vec] = ...
        get_probe_pm_mm_columns(probe_intens_path, ...
        [hind_txt_file_names{i}(1:end-4) '.mat'], hind_txt_file_names{i}(1:end-4));

    allele_ratio_vec  = calc_probe_allele_ratio(probe_intens_mat, 1, ...
        pa_vec, ma_vec, pb_vec, mb_vec);
    eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) ...
        'affybench_hind.mat allele_ratio_vec;']);
    % Call Matlab normalization (currently this normalization is problematic ..)
    [NormData MedStructure] = affyinvarsetnorm([baseline_pm_column' cur_pm_column'], 'baseline', 1);
    NormData = mean( (reshape(NormData(:,2), length(NormData) / length([pa_vec pb_vec]), length([pa_vec pb_vec])))' );
    eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) 'Norm.mat NormData;']);


    if(i == 1)
        MeanNormData = NormData; % Just copy to the correct name
    else
        MeanNormData = MeanNormData + NormData; % Just copy to the correct name
    end
end

MeanNormData = MeanNormData ./ length(hind_txt_file_names); % calc mean
eval(['save ' probe_intens_path 'MeanNormData.mat MeanNormData;']);

eval(['load ' probe_intens_path 'MeanNormData.mat;']); % just

for i=1:length(hind_txt_file_names)   % 2nd loop: Normalize and transfer to application input format
    i

    eval(['load(''' probe_intens_path hind_txt_file_names{i}(1:end-4) 'Norm.mat'');']);
    copy_num_vec = (2 .* NormData ./ MeanNormData)'; % Just copy to the correct name and transpose!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute A/B ratio
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% x = copy_num_vec'.*allele_ratio_vec ./ (1+allele_ratio_vec);
    %%    eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) 'Ratio.mat ratio_vec;']);

    genotype_vec = genotypes_hind_num(:,i);    % Get the genotype

    eval(['load ' probe_intens_path hind_txt_file_names{i}(1:end-4) 'affybench_hind.mat;']); % load the allele_ratio_vec
    if(size(copy_num_vec,1) == 1)
        copy_num_vec = copy_num_vec';
    end
    eval(['save ' probe_intens_path hind_txt_file_names{i}(1:end-4) ...
        'affybench_hind.mat snp_ids genotype_vec copy_num_vec allele_ratio_vec;']);
    
    PlotAlleleRatios(max(min(copy_num_vec,5),0), allele_ratio_vec, genotype_vec, 'allele probe intens Check');
end




