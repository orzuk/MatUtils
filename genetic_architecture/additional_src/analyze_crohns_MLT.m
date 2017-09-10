% Analyze crohn's data as if it was from an MLT model
%
% Input:
% N - # of pathways in MLT model
% crohns_mu -  assumed prevalence
% crohns_loci_file - file with loci
% crohns_output_file - where to save output
%
function analyze_crohns_MLT(N, crohns_mu, crohns_loci_file, crohns_output_file)

AssignGeneralConstants;
[~, ~, html_outdir] = get_machine_type();

if(~exist('crohns_mu', 'var') || isempty(crohns_mu))
    crohns_mu = 0.002; % estimated prevalence of crohnss
end
if(~exist('crohns_loci_file', 'var') || isempty(crohns_loci_file))
    crohns_loci_file = fullfile(html_outdir, ...
        'common_disease_broad_catalog\Crohn_s_disease\Crohn_s_disease_table.txt');
    % % %     switch machine
    % % %         case PC
    % % %             crohns_loci_file = 'C:\public_html\data\common_disease_broad_catalog\Crohn_s_disease\Crohn_s_disease_table.txt';
    % % %         case UNIX
    % % %             crohns_loci_file = '~orzuk/public_html/data/common_disease_broad_catalog/Crohn_s_disease/Crohn_s_disease_table.txt';
    % % %     end
end
if(~exist('N', 'var') || isempty(N))
    N = 3; % number of pathways in the MLT model
end
if(~exist('crohns_output_file', 'var') || isempty(crohns_output_file))
    %    crohns_output_file = 'temp_crohns_file.txt';
    crohns_output_file = ...
        '../../common_disease_model/docs/pnas/genetic_interactions/tables/crohns_data_overestimation.txt';
end

crohns_h = 0.5; % estimated heritability of crohns (under the LT model)
crohns_lambda_s = 10; % estimated from sib study
%crohns_h = familial_risk_to_heritability(crohns_lambda_s, 'liability', crohns_mu, 0.5);
k = 1; % number of liabilities we need to pass
isoheritability_flag=1;


R = loadcellfile(crohns_loci_file,1, 9); % separate by tabs
R{1,1}='\#';
R{1,7} = 'GRR'; % change OR to GRR
R = strrep_cell(R, 'rs_', 'rs');
R = strrep_cell(R, ' ', '');
R = strrep_cell(R, '"', '');
for i=2:length(R)-1
    R{i,5} = str2word('(', str2word(',', R{i,5}, 1), 1);
end


crohns_loci_RAF = cell2mat(R(2:end-1,6));
crohns_loci_GRR = cell2mat(R(2:end-1,7));
crohns_var_explained = cell2mat(str2num_cell(strrep(R(2:end-1,15), '%', '')));

num_snps = length(crohns_loci_RAF);
[new_var_explained_total new_var_explained] = ...
    genetic_relative_risk_to_variance_explained( ...
    crohns_loci_RAF, crohns_loci_GRR, crohns_mu, 'diploid')
% [new_var_explained_logit_total new_var_logit_explained] = ...
%     genetic_relative_risk_to_variance_explained( ...
%     crohns_loci_RAF, crohns_loci_GRR, crohns_mu, 'diploid', 'logit')

options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance


num_epidemiological=4;
% Summary with many different estimates
W_header = {'$\mu (%)$', '$\lambda_{MZ}$', '$\lambda_s$', ... % '$h_{pop}^2  (%)$', ...
    'Model', '$h_{pathway}^2 (%)$', '$c_R (%)$', ...
    '$h_s^2 (%)$', '$h_{all}^2  (%)$', '$V_c (%)$', ...
    '$\pi_{phantom} (%)$', '$\pi_{explained} (%)$'}; % new! add shared environment!!! 
W = cell(num_epidemiological*2, length(W_header)); % each epi params have two models, A_delta and LP_delta

for crohn_epi_params = 1:num_epidemiological % try many different parameters
    switch crohn_epi_params % different epidemiological parameters for Crohn's disease, obtained from different sources
        case 1
            main_example_ind = crohn_epi_params;
            crohns_mu = 0.001; % estimated prevalence of crohns
            crohns_h = 0.5; % estimated heritability of crohns (under the LT model)
            %crohns_lambda_s = 10; % estimated from sib study
            crohns_lambda_s = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 0.5)
            crohns_lambda_mz = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 1)
            
        case 2
            crohns_mu = 0.001; % estimated prevalence of crohns
            crohns_lambda_s = 27; % estimated from sib study
            crohns_h = familial_risk_to_heritability(crohns_lambda_s, 'liability', crohns_mu, 0.5); % estimated heritability of crohns (under the LT model)
            crohns_lambda_s = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 0.5)
            crohns_lambda_mz = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 1)
            
        case 3
            crohns_mu = 0.002; % estimated prevalence of crohns
            crohns_h = 0.5; % estimated heritability of crohns (under the LT model)
            %crohns_lambda_s = 10; % estimated from sib study
            crohns_lambda_s = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 0.5)
            crohns_lambda_mz = heritability_to_familial_risk(crohns_h, 'liability', crohns_mu, 1)
            
        case 4
            crohns_mu = 0.002; % estimated prevalence of crohns
            crohns_lambda_s = 35; % estimated from sib study
            crohns_lambda_mz = 250; % corresponds for MZ risk of one half
            crohns_h = familial_risk_to_LP_parameters( ... % try model with shared environment
                [crohns_lambda_mz crohns_lambda_s], 1, [1 0.5], 'binary', crohns_mu); % estimated heritability of crohns (under the LT model)
    end
    [new_var_explained_total new_var_explained] = ...
        genetic_relative_risk_to_variance_explained( ...
        crohns_loci_RAF, crohns_loci_GRR, crohns_mu, 'diploid'); % compute variance explained
    h_s = new_var_explained_total;
    
    
    r_MZ = familial_risk_to_heritability(crohns_lambda_mz, 'liability', crohns_mu, 1);
    r_DZ = familial_risk_to_heritability(crohns_lambda_s, 'liability', crohns_mu, 0.5)/2;
    h_pathway_additive(N) = 2*(r_MZ-r_DZ);
    c_R_additive(N) = (2*r_DZ-r_MZ) / h_pathway_additive(N);
    
    crohns_h_shared_env = c_R_additive(N) * (1-h_pathway_additive(N)); 
    if(crohns_h_shared_env < 0.00001)
        crohns_h_shared_env = 0;
    end
    
    W{crohn_epi_params*2-1, 1} = num2str(crohns_mu*100,3); % prevalence
    W{crohn_epi_params*2-1, 2} = num2str(crohns_lambda_mz,3); % MZ risk
    W{crohn_epi_params*2-1, 3} = num2str(crohns_lambda_s,3); % DZ risk
    W{crohn_epi_params*2-1, 4} = '$A_{\Delta}$'; % model string
    W{crohn_epi_params*2, 4} = '$LP_{\Delta}(3)$'; % model string
    
    crohns_r_mz = familial_risk_to_heritability(crohns_lambda_mz, 'liability', crohns_mu, 1); % convert to liability scale 
    
    for ancor_lambda={'MZ', 'DZ'} % decide which lambda_s to use as 'anchor'
        for N=3:3 % 1:5 % loop on different LP models
            crohns_mu_l = fminbnd(@(x) abs(binocdf(k-1, N, x)-(1-crohns_mu)), 0, 1, options); % find mu_l that keeps the PREVALENCE
            
            switch ancor_lambda{1} % which lambda_s is 'anchor' for heritability calculations
                case 'MZ'
                    crohns_h_x_one_liab = fminbnd(@(x) abs(crohns_lambda_mz - ...
                        compute_k_of_N_liabilities_statistics( ...
                        N, k, crohns_mu_l, x, 0, 1)), 0, 1, options); % compute only for MZ to save time
                case 'DZ'
                    crohns_h_x_one_liab = fminbnd(@(x) abs(crohns_lambda_s - ...
                        compute_k_of_N_liabilities_statistics( ...
                        N, k, crohns_mu_l, x, 0, -2)), 0, 1, options); % compute only for MZ to save time
            end
            
            [lambda_R S_stats] = compute_k_of_N_liabilities_statistics( ...
                N, k, crohns_mu_l, crohns_h_x_one_liab, 0, []); % h_x_one_liab
            lambda_mz(N) = lambda_R(1); lambda_s(N) = lambda_R(2);
            crohns_r_mz_LP(N) = familial_risk_to_heritability(lambda_mz(N), 'liability', crohns_mu, 1); % convert to liability scale

            h_liab_loci(N) = S_stats.h_liab_loci;
            h_liab_twins(N) = S_stats.h_liab_twins;
            %            c_R_LP(N) = S_stats.h_shared_env / (1-S_stats.h_liab_twins);
            [h_pathway_fitted_LP(N) c_R_LP(N)] = ...
                familial_risk_to_LP_parameters([crohns_lambda_mz crohns_lambda_s], ...
                N, [1 0.5], 'binary', crohns_mu); % Fit LP model
            [lambda_R2 S_stats2] = compute_k_of_N_liabilities_statistics( ...
                N, k, crohns_mu_l, h_pathway_fitted_LP(N), 0, []); % h_x_one_liab again. Here there is no shared environment
            crohns_r_mz_no_shared_env_LP(N) = familial_risk_to_heritability(lambda_R2(1), 'liability', crohns_mu, 1); % convert to liability scale
            crohns_h_shared_env_LP(N) = crohns_r_mz_LP(N) - crohns_r_mz_no_shared_env_LP(N);
            
            
            h_all_LP(N) = S_stats2.h_liab_loci;
            
            if(abs(c_R_LP(N)) < 0.000001)
                c_R_LP(N)=0;
            end
            if(abs(c_R_additive(N)) < 0.000001)
                c_R_additive(N)=0;
            end
            
            [lambda_R_100 S_stats_100] = compute_k_of_N_liabilities_statistics(...
                N, k, crohns_mu_l, 1+0*crohns_h_x_one_liab, 0, []); % assume heritability of 100% % h_x_one_liab
            h_phantom_100 = S_stats_100.h_liab_unexplained_gap;
            h_phantom(N) = (crohns_h-h_all_LP(N))/crohns_h;
            
            S_stats2.h_liab_unexplained_gap;
            if(N==1)
                h_phantom(N)=0;
            end
        end
        if(crohn_epi_params == main_example_ind) % copy 'correct' values of var. explained !
            crohns_var_explained = new_var_explained;
            crohns_example_h_additive = crohns_h;
            crohns_example_h_LP = h_all_LP(N);
        end
        
        h_explained = sum(new_var_explained) ./ h_liab_loci;
        T = [(1:N)' repmat(100*crohns_mu,N,1) lambda_mz' lambda_s' 100*h_liab_twins' 100*h_liab_loci' ...
            100*h_phantom' 100*h_explained'];
        
        if(~exist('T_all', 'var'))
            T_all = [{'N', '\mu (%)', '\lambda_MZ', '\lambda_s', 'h_{pop} (%)', 'h_{all-loci} (%)', ...
                '\pi_{phantom} (%)', 'h_{explained} (%)'}' num2str_cell(num2cell(T'), 3)]'
        else
            T_all = [T_all' repmat({''}, 8, 1) num2str_cell(num2cell(T'), 3)]'
        end
        ctr=4;    % Additive model
        %        W{crohn_epi_params*2-1, ctr} = num2str(crohns_h *100,3); ctr=ctr+2; % h_pop^2
        ctr=ctr+1;
        W{crohn_epi_params*2-1, ctr} = num2str(h_pathway_additive(N)*100,3); ctr=ctr+1; % fitted h_pathway
        W{crohn_epi_params*2-1, ctr} = num2str(c_R_additive(N)*100,3); ctr=ctr+1; % fitted c_R
        W{crohn_epi_params*2-1, ctr} = num2str(h_s*100,3); ctr=ctr+1; % h_s^2 of the 74 loci
        W{crohn_epi_params*2-1, ctr} = num2str(crohns_h*100,3); ctr=ctr+1; % h_all^2
        W{crohn_epi_params*2-1, ctr} = num2str(crohns_h_shared_env*100,2); ctr=ctr+1; % NEW! h-shared environment
        W{crohn_epi_params*2-1, ctr} = '0'; ctr=ctr+1; % pi_phantom
        W{crohn_epi_params*2-1, ctr} = num2str((h_s/crohns_h)*100,3); % pi_explained
        
        ctr=4; % LP(3) model
        %        W{crohn_epi_params*2, ctr} = num2str(crohns_h*100,3);
        ctr=ctr+1; % 2; % h_pop^2
        W{crohn_epi_params*2, ctr} = num2str(h_pathway_fitted_LP(N)*100,3); ctr=ctr+1; % fitted h_pathway
        W{crohn_epi_params*2, ctr} = num2str(c_R_LP(N)*100,3); ctr=ctr+1; % fitted c_R
        W{crohn_epi_params*2, ctr} = num2str(h_s*100,3); ctr=ctr+1; % h_s^2 of the 74 loci
        W{crohn_epi_params*2, ctr} = num2str(h_all_LP(N)*100,3); ctr=ctr+1; % h_all^2
        W{crohn_epi_params*2, ctr} = num2str(crohns_h_shared_env_LP(N)*100,3); ctr=ctr+1; % NEW! h-shared environment
        W{crohn_epi_params*2, ctr} = num2str(h_phantom(N)*100,3); ctr=ctr+1; % pi_phantom
        W{crohn_epi_params*2, ctr} = num2str(h_s/h_all_LP(N)*100,3); ctr=ctr+1; % h explained
                
    end % loop on anchor lambda
    if(crohn_epi_params < num_epidemiological)  % add a dashed line
        W{crohn_epi_params*2, end} = [W{crohn_epi_params*2, end} '   \hdashline'];
    end
    savecellfile(T_all, ['tmp_crohns_lambda_analysis_epi_' ...
        num2str(crohn_epi_params) '.txt'],[],1);
end % loop on different epidemiological parameters




crohns_var_explained_liability = new_var_explained;
% for i=1:1 % num_snps
%     convert_i = i
%     crohns_var_explained_liability(i) = ...
%     heritability_scale_change_MLT(new_var_explained(i), 1, N, crohns_mu, 'MLT');
%     MLT_ratio = crohns_var_explained_liability(i) / new_var_explained(i);
% end
MLT_ratio = 1/(1-h_phantom_100);
crohns_var_explained_liability(2:end) = new_var_explained(2:end) * MLT_ratio;


temp_heritability_one_liab = familial_risk_to_heritability(lambda_R(1), ...
    'liability', crohns_mu, 1);
lambda_S_one_liab = heritability_to_familial_risk(temp_heritability_one_liab, ...
    'liability', crohns_mu, 0.5);


[~, sort_perm] = sort(crohns_var_explained, 'descend'); % sort loci by variance explained

for i=2:length(R)-1
    for j=[6 7]
        R{i,j} = sprintf('%0.2f', R{i,j});
    end
end

S = [R(1:end-1, [2 5]) num2str_cell(R(1:end-1,  [6 7]), 2) ...
    ['$V_i (%)$' num2cell(100*crohns_var_explained)']' ... % variance explained
    ['$% h_{all}^2$ expl. (\LT)' num2cell(100*crohns_var_explained./crohns_h)']' ...
    ['$% h_{all}^2$ expl. (\LPd)' num2cell(100*crohns_var_explained./h_all_LP(N))']' ];
% ['Var-expl. (LT)' num2str_cell(num2cell(100*new_var_explained), 2, [], '%')']' ...
% ['Var-expl. (LP)' num2str_cell(num2cell(100*crohns_var_explained_liability), 2, [], '%')']']; % take snp name, gene, allele freq, grr, h^2 explained, beta (?)
S = strrep_cell(S, '', '');

for i=2:length(S) % loop on all snps
    for j=[5 6 7] % var. explained
        S{i,j} = sprintf('%0.2f', S{i,j});
    end
end



pathway_vec = repmat(1:N, 1, ceil(num_snps/N)); pathway_vec = pathway_vec(1:num_snps);
for i=1:num_snps  % split SNP to pathway by their name
    run_snp = i
    first_chr = upper(R{i+1,5});
    first_chr = first_chr(1:min(1,length(first_chr)));
    switch first_chr
        case {'A','B','C','D','E','F','G','H','I'}
            pathway_vec(i)=1;
        case {'J','K','L','M','N','O','P','Q'}
            pathway_vec(i)=2;
        case  {'R','S','T','U','V','W','X','Y','Z'}
            pathway_vec(i)=3;
        otherwise
            pathway_vec(i) = ceil(rand(1)*3);
    end % switch
end
pathway_vec = ['Pathway' num2cell(pathway_vec)]';

S = [S pathway_vec]; % add pathway-assignment, GRR for each 'sub-disesase', heritability for each subdisease
S = S([1 1+sort_perm']',:); % sort according to var. explained
S = [R(1:end-1, 1) S]; % add index

% Add a line at the bottom summing up the contribution
S_bottom = {'Total (%):', '', '', '', '', ...
    [num2str(100*sum(crohns_var_explained), 3)] ...
    [num2str(100*sum(crohns_var_explained./crohns_example_h_additive), 3)] ...
    [num2str(100*sum(crohns_var_explained./crohns_example_h_LP), 3)], ''};

S = [S' S_bottom']'; % add total
S{end-1,end} = [num2str(S{end-1,end}) '   \hdashline']; % add dashed before total



%[num2cell(crohns_loci_RAF) num2cell(crohns_loci_GRR) num2]; % save loci

savecellfile(S, [remove_suffix_from_file_name(crohns_output_file) '_tab.txt']); % save results - this should be used in supp. info.

S_latex = latex(S, 2, precision);
S_latex(end-2,:) = strrep(S_latex(end-2,:), '\hdashline \\', '\\ \hdashline'); % add dashed before total

S_latex = mat2cell(S_latex, ones(size(S_latex,1),1));
S_latex = strrep_cell(S_latex, '%', '\%');
S_latex{1} = strrep(S_latex{1},  '|c', '|r'); % align to left
S_latex{1} = strrep(S_latex{1},  '{|r|r|r|', '{|l|l|l|'); % align first column to left

tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', '{\scriptsize', ...
    ['\caption[Phantom heritability for Crohn''s disease]{Assumed and actual variance explained by associated SNPs ' ...
    ' for Crohn''s disease. (Source: Franke et al. \cite{franke2010genome}). Epidemiological Parameters: ' ...
    '$\mu=' num2str(100*crohns_mu) '\%, h_{pop}^2=' num2str(crohns_h*100,2) '\%,' ...
    ' \lambda_{MZ}=' num2str(crohns_lambda_mz,3) ', \lambda_s=' num2str(crohns_lambda_s,3) '$.' ...
    ' \label{table:crohns_disease}}'] S_latex{1}};
tab_footer = {'}',  ...
    '\end{center}', '\end{table}', '', '\newpage', ''};



% tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
%     '\begin{center}', '{\scriptsize', ...
%     ['\caption[Model parameters for disease traits]{Model parameters for ' ...
%     ' the LP model $\archd(k,\hx,c,\mu)$ for disease traits. \label{table:fig2_data_disease_traits}}'], ...
%     S_latex{1}};


S_latex = split_latex_table(S_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 50, -1);
savecellfile(S_latex, crohns_output_file, [], 1); % save results - this should be used in supp. info.



% Summary with many different estimates
W = [W_header' W']'; % concatenae header

W_latex = latex(W, 2, precision); % New: add a summary table for different heritability estimates

W_latex = mat2cell(W_latex, ones(size(W_latex,1),1));
W_latex = strrep_cell(W_latex, '%', '\%');
W_latex{1} = strrep(W_latex{1},  '|c', '|r'); % align to left
W_latex{1} = strrep(W_latex{1},  '{|r|r|r|r|', '{|r|r|r|c|'); % align fourth (fifth) column to center

tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', '{\scriptsize', ...
    ['\caption[Phantom heritability for Crohn''s disease for different epidemiological parameters]{Assumed and actual variance explained by all 74 associated SNPs ' ...
    ' for Crohn''s disease. (Source: Franke et al. \cite{franke2010genome})' ...
    ' \label{table:crohns_disease_heritability_different_epidemiological_parameters}}']};
tab_footer = {'}',  ...
    '\end{center}', '\end{table}', '', '\newpage', ''};

for i=1:length(W_latex)
    W_latex{i} = strrep(W_latex{i}, '\hdashline \\', '\\ \hdashline'); % add dashed before total
end
W_latex = [tab_header W_latex' tab_footer]';

savecellfile(W_latex, ...
    [remove_suffix_from_file_name(crohns_output_file) '_summary.txt'], [], 1); % save results - this should be used in supp. info.


%


