% A script for computing BMI effect sizes
% Parsing supp. tab. 8 from the paper [Speliotes, et al., Nature Genetics, 2010]

% Parse file
S = ReadDataFile('../../common_disease_model/data/BMI/Speliotes_Standard_Effect_Size_supp_info_tab8.txt');
M = ReadDataFile('../../common_disease_model/data/BMI/Speliotes_Standard_Effect_Size_main_tab1.txt');
T = ReadDataFile('../../common_disease_model/data/BMI/BMILociFromThreeTraitsFile.txt');
T.Effect_Size_Units = empty_cell_to_empty_str(T.Effect_Size_Units);
bmi_inds = strmatch('kg', T.Effect_Size_Units');
M.VarExplained = cell2mat(str2nums_cell(M.VarExplained))./100; % transfer var. explained to numeric and fraction (from %) 
M.ComputedBeta = sqrt(M.VarExplained ./ (2.*M.RAF .* (1-M.RAF)));
figure; hold on; plot(M.Beta(:,1), M.ComputedBeta, '.'); % title('Reported vs. standardized compute \beta'); 
xlim([0 0.4])
fit_params = polyfit(M.Beta(:,1), M.ComputedBeta, 1); 
plot(M.Beta(:,1), M.Beta(:,1).*fit_params(1) + fit_params(2), 'r'); 



[intersect_snps I J] = intersect(S.SNP, M.SNP)
figure; plot(M.Beta(J,1), S.MeanDifference(I), '.'); title('Reported effect sizes: main vs. supp. info.');
figure; hold on; plot(M.ComputedBeta(J,1), S.MeanDifference(I), '.'); title('Computed effect sizes: main vs. supp. info.');
plot(0:0.01:0.1, 0:0.01:0.1, 'r');


[intersect_snps I J] = intersect(S.SNP, T.SNP_ID_(bmi_inds))
replication_effect_sizes =  cell2mat(T.Effect_Size_Followup(bmi_inds(J)))

figure; plot(S.MeanDifference(I), replication_effect_sizes(:,1), '.'); title('Compare Effect Sizes');

study_MAF = cell2mat(T.Risk_Allele_Frequency(bmi_inds(J)))
figure; plot(min(S.MAF(I), 1-S.MAF(I)), min(study_MAF, 1-study_MAF), '.'); title('Compare RAF');


% Compute var. explained based on S (supp. table)
V = sum(2 .* S.MAF .* (1-S.MAF) .* S.MeanDifference.^2)


2 .* S.MAF .* (1-S.MAF) .* S.MeanDifference.^2 - S.StandardEffectSize
