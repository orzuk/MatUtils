% Compute Crohns variance explained from from Daly's paper: 

R = loadcellfile('../../common_disease_model/data/Crohn_Mark_Daly/TableAllLociFromMarkDaly.txt');

Y = loadcellfile('../../common_disease_model/data/Crohn_Mark_Daly/var_explained.txt');

for i=1:length(Y)
    YY{i} = str2word(']',Y{i}, 2);
end


X = str2nums_cell(R);

num_snps = length(X)-1;

for i=1:num_snps
    snp_ids{i} = str2word(' ', R{i+1}, 1);
    OR(i) = X{i+1}(2);
    CAF(i) = X{i+1}(3);
    set_ind(i) = X{i+1}(4);
end

prev = 0.002;
GRR = odds_ratio_to_genetic_relative_risk(OR, CAF, prev);
[lambda_s_vec ...
    lambda_s_add lambda_mz_add h_add V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab h_liab_vec] = ...
    genetic_relative_risk_to_heritability(mat2vec([CAF' CAF']'), mat2vec([OR' OR']'), prev); 
total_h = h_liab
total_h = h_liab*2

[lambda_s_vec ...
    lambda_s_add lambda_mz_add h_add V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab h_liab_vec] = ...
    genetic_relative_risk_to_heritability(mat2vec([CAF' CAF']'), mat2vec([GRR' GRR']'), prev); 
total_h_GRR = h_liab*2


RR = [snp_ids' num2cell([OR' CAF' set_ind']) YY'];
RR = [{'SNP-ID' 'OR' 'CAF' 'DISCOVERY', 'Var-Expl.'}' RR']'; 
savecellfile(RR, '../../common_disease_model/data/Crohn_Mark_Daly/TableAllLociFromMarkDaly_tab.txt'); 

savecellfile(RR, '../../common_disease_model/data/Crohn_Mark_Daly/CrohnsLoci_with_var_expl.txt'); 

    