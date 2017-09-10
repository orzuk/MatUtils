% Analyze the rare-alleles appearing in the paper on MTNR1B rare-variants on T2D

%function parse_T2D_MTNR1B()

Ji_BloodPressure_rare_file = '../../common_disease_model/data/rare/Ji_Lifton_BloodPressure/SuppInfoTable2.txt';

R = ReadDataFile(Ji_BloodPressure_rare_file, [], [], [], 9); 
R.AlleleFreq = strrep_cell(R.AlleleFreq, '>', '');
R.AlleleFreq = str2num_cell(R.AlleleFreq);
R.AlleleFreq = empty_cell_to_numeric_val(R.AlleleFreq, 0); 
R.AlleleFreq = my_cell2mat(R.AlleleFreq);

MAF = 0.1; % maximum allele frequency 
gene_names = {'NCCT', 'NKCC2', 'ROMK'};
%             SLC12A3 SLC12A1 KCNJ1
%R.Carriers = cell2mat(R.Carriers);

for i=1:length(gene_names)
    gene_inds = strfind_cell(R.Gene, gene_names{i});
    gene_inds = intersect(gene_inds, find(1-R.Synonymous)); 
    gene_inds = intersect(gene_inds, find(R.AlleleFreq <= MAF)); % restrict allele frequency 
    num_alleles(i) = length(gene_inds);
    num_carriers(i) = sum(R.Carriers(gene_inds));
    num_null_alleles(i) = sum(R.IsNull(gene_inds));
    num_null_carriers(i) = sum(R.Carriers(gene_inds) .* R.IsNull(gene_inds));    
end
    

%function [S R] = ReadDataFile(data_file, output_file, cell_to_mat, skip_lines, delimiter, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NEW: Get Table 3 from the Ji et al. paper: 
AssignGeneralConstants;
carriers_mat = [3 1 6 3 10 3 13 5 20 7; ...
     2 0 2 1 3 1 5 1 9 2; ...
     4 0 5 1 6 1 9 1 11 1];
noncarriers_vec = [304 312 613 621 920 934 1225 1245 1633 1443]; 
% quantiles_vec = [10 90 20 80 30 70 40 60 50 50] ./ 100; 
quantiles_vec = [10 10 20 20 30 30 40 40 50 50] ./ 100; 

carriers_mat = [carriers_mat' sum(carriers_mat)']';
numcases_vec = noncarriers_vec + carriers_mat(4,:);

num_samples = sum(numcases_vec(9:10))
total_carriers_per_gene = sum(carriers_mat(:,9:10),2)


MAF = 0.5 .* sum(0.5 .* carriers_mat(:, 9:10) ./ repmat(noncarriers_vec(9:10), 4, 1), 2); % for each gene


R_genes = {'n', 'SLC12A3', 'SLC12A1', 'KCNJ1', 'All 3 Genes'}'; 
R_header = {'Gene', 'Cases' , 'Unaffecteds', 'CAF',  '$1 + \lambda$', 'V'};
%R_mat = repmat(R_genes, 10, 1)



num_carriers_cases_vec = mat2vec([numcases_vec' carriers_mat']')
num_carriers_controls_vec = mat2vec([num_samples-numcases_vec' ...
    repmat(total_carriers_per_gene, 1, 10)' - carriers_mat']')  % solve this later 
numcontrols_vec = num_samples - numcases_vec;
MAF_vec = repmat([-1 MAF']', 10, 1);

lambda_mat = (0.5*carriers_mat ./ repmat(numcases_vec, 4, 1)) ./ repmat(MAF, 1, 10) - 1
lambda_inv_mat = (0.5 .* (repmat(total_carriers_per_gene, 1, 10) - carriers_mat) ./ ...
    repmat(numcontrols_vec, 4, 1) ) ./ repmat(MAF, 1, 10) - 1; 
%V_mat = lambda_mat.^2 .* repmat(2.*MAF, 1, 10) .* repmat(quantiles_vec ./ (1-quantiles_vec), 4, 1) 
V_mat = lambda_inv_mat.^2 .* repmat(2.*MAF, 1, 10) .* repmat((1-quantiles_vec) ./ quantiles_vec, 4, 1)



new_lambda_vec = (0.5.*num_carriers_controls_vec ./ mat2vec(repmat(num_samples-numcases_vec, 5, 1))) ./ MAF_vec
new_lambda_cases_vec = (0.5.*num_carriers_cases_vec ./ mat2vec(repmat(numcases_vec, 5, 1))) ./ MAF_vec
new_lambda_vec([1:5 11:15 21:25 31:35 41:45]) = new_lambda_cases_vec([1:5 11:15 21:25 31:35 41:45]); % copy cases 
new_lambda_vec=new_lambda_vec-1;

one_plus_lambda_vec = new_lambda_vec+1; 
%one_plus_lambda_vec(one_plus_lambda_vec < 1) = -1;
%one_plus_lambda_vec = mat2vec([repmat(-2, 10, 1)  lambda_inv_mat']')+1;
V_vec = mat2vec([repmat(-1, 10, 1)  V_mat']'); 

prev_vec = repmat(quantiles_vec, 5, 1); prev_vec(:,[2 4 6 8 10]) = 1-prev_vec(:,[2 4 6 8 10]);
prev_vec = mat2vec(prev_vec);
V_vec2 =  2.*MAF_vec .* new_lambda_vec.^2 .* prev_vec ./ (1-prev_vec); 
V_vec2(1:5:end) = -0.01; 
R_mat = [num2str_cell(num2cell([num_carriers_cases_vec num_carriers_controls_vec])) ...
    [num2cell(MAF_vec.*100) num2str_cell(num2cell(one_plus_lambda_vec), 3) num2cell(V_vec2.*100)]]; 


R_cell = [repmat(R_genes, 10, 1) R_mat]; 
R_cell = [R_header' R_cell']'; 

R_cell = num2str_cell(R_cell, 2); 
R_cell = strrep_Cell(R_cell, '-1', '');
R_cell = strrep_Cell(R_cell, '-0.5', '');
R_cell = strrep_Cell(R_cell, 'e+002', '');

for i=1:10 
    if(mod(i,2) == 1)
        R_cell{5*i-3, 2} = [R_cell{5*i-3, 2} ' (bottom $' num2str(quantiles_vec(i)*100) '\%$)'];
        R_cell{5*i-3, 3} = [R_cell{5*i-3, 3} ' (top $' num2str((1-quantiles_vec(i))*100) '\%$)'];
    else
        R_cell{5*i-3, 2} = [R_cell{5*i-3, 2} ' (top $' num2str(quantiles_vec(i)*100) '\%$)'];
        R_cell{5*i-3, 3} = [R_cell{5*i-3, 3} ' (bottom $' num2str((1-quantiles_vec(i))*100) '\%$)'];
    end
    
    for j=1:4
       R_cell{5*i+j-3, 4} = [ '$' R_cell{5*i+j-3, 4} '\%$' ]; 
       R_cell{5*i+j-3, 6} = [ '$' R_cell{5*i+j-3, 6} '\%$' ]; 
    end
end
savecellfile(R_cell, 'Ji_et_al_tab.txt'); 

%R_cell = R_cell([1 7:11 17:21 27:31 37:41 47:51],:); % Take only top half - why? 
R_cell = R_cell([1 7:11 17:21 27:31 37:41 47:51 ... 
    2:6 12:16 22:26 32:36 42:46],:); % Permute. First top, then bottom half.  


R_cell_latex = latex(R_cell, 2, precision);
R_cell_latex = mat2cell(R_cell_latex, ones(size(R_cell_latex,1),1));
R_cell_latex{1} = strrep(R_cell_latex{1},  '|c', '|r'); % align to right
R_cell_latex{1} = strrep(R_cell_latex{1},  '{|r', '{|l'); % align to left

for i=1:(floor(length(R_cell_latex)/5)-1) % Put lines between different percentiles 
    R_cell_latex{i*5+2} = [R_cell_latex{i*5+2} ' \hline'];
end


tab_header = {'\begin{table}[h!] % table with all precentiles', ...
    '\begin{center}', '{', ...
    ['\caption['']' ...
    '{}'], ...
    R_cell_latex{1}};
tab_footer = {'}',  ...
    '\label{tab:Ji_three_genes}', ...
    '\end{center}', '\end{table}', '', ''};
% S_latex = [tab_header  R_latex(2:end)' tab_footer]';
R_cell_latex = split_latex_table(R_cell_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 25); % Split to two pages


savecellfile(R_cell_latex, 'Ji_et_al_tab_latex.txt', [], 1);



     