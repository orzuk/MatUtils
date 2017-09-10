% Take sib-risk values from Wray's paper
function [disease_vec disease_short_vec ...
    lambda_s_vec lambda_mz_vec prevalence_vec ...
    liab_assum_lambda_s_vec liab_assum_lambda_mz_vec ...
    h_from_sibs_vec h_from_mz_vec lambda_s_overestimation_vec] = wray_lambda_table(figs_flag)

AssignGeneralConstants;

output_dir = '../../common_disease_model/docs/pnas/genetic_interactions/tables/';

if(exist('lambda_s_vec', 'var') || exist('wray_disease_lambda_list.mat', 'file'))
    compute_params_flag = 0;
    load('wray_disease_lambda_list.mat');
else
    compute_params_flag = 1;
end

%compute_params_flag = 1;
    figs_dir = '../../common_disease_model/figs/';
if(compute_params_flag)
    
    compute_params_flag = 0; % make sure next time
    lambda_s_vec = [1.3 2.1 3.2 2.2 3.5 3.4 3.5 7 8.6 14 20 64 82 29 65 ...
        5.7786 1.5102 2.1085]; % last 3 values are for MLT
    lambda_mz_vec = [2 4.7 4.6 4.1 10.4 6.6 12.2 60 52.1 79 190 600 630 29^2 774 ...
        32.0225 2.4140 4.9239];  % last 3 values are for MLT
    prevalence_vec = [0.24 0.12 0.056 0.036 0.028 0.19 0.01 0.01 ...
        0.0085 0.005 0.001 0.001 0.001 0.001 0.0003 ...
        0.01 0.1153 0.1];  % last 3 values are for MLT
    
    disease_vec = {'Depression', ...
        'AMD',  ...
        'Myocardial infarction', ...
        'Breast cancer', ...
        'Type 2 diabetes',  ...
        'Asthma',  ...
        'Rheumatoid Arthritis', ...
        'Bipolar Disorder',  ...
        'Schizophrenia',  ...
        'Type 1 diabetes',  ...
        'Multiple sclerosis',  ...
        'Crohn�s',  ...
        'Ankylosis spondylitis', ...
        'Lupus^{(1)}', ...
        'Lupus^{(2)}', ...
        '{\bf MLT1}', '{\bf MLT2}', '{\bf MLT3}'}
    
    disease_short_vec = {'Dep.', ...
        'AMD',  ...
        'MI', ...
        'BC', ...
        'T2D',  ...
        'Ast.',  ...
        'RA', ...
        'BD',  ...
        'Sch.',  ...
        'T1D',  ...
        'MS',  ...
        'Cro.',  ...
        'AS', ...
        'Lup.^{(1)}', ...
        'Lup.^{(2)}', ...
        'MLT1', 'MLT2', 'MLT3'}
    %    '{\bf MLT1}', '{\bf MLT2}', '{\bf MLT3}'}
    
    
    
    N = length(lambda_s_vec)
    h_from_mz = familial_risk_to_heritability(50, ...
        'liability', 0.99, 1)
    lambda_s_from_mz = heritability_to_familial_risk(h_from_mz, ...
        'liability', 0.99, 0.5)
    
    h_from_sibs_vec = zeros(1,N); h_from_mz_vec = zeros(1,N);
    liab_assum_lambda_mz_vec = zeros(1,N); liab_assum_lambda_s_vec = zeros(1,N);
    for i=1:N
        h_from_sibs_vec(i) = familial_risk_to_heritability(lambda_s_vec(i), ...
            'liability', prevalence_vec(i), 0.5)
        h_from_mz_vec(i) = familial_risk_to_heritability(lambda_mz_vec(i), ...
            'liability', prevalence_vec(i), 1)
        liab_assum_lambda_mz_vec(i) = heritability_to_familial_risk(h_from_sibs_vec(i), ...
            'liability', prevalence_vec(i), 1)
        liab_assum_lambda_s_vec(i) = heritability_to_familial_risk(h_from_mz_vec(i), ...
            'liability', prevalence_vec(i), 0.5)
        % % % %     k=1;
        % % % %     for try_mu = 0.001:0.01:(1/ lambda_mz_vec(i))
        % % % %      h_sib(k) = familial_risk_to_heritability(lambda_s_vec(i), ...
        % % % %         'liability', try_mu, 0.5);
        % % % %      h_mz(k) = familial_risk_to_heritability(lambda_mz_vec(i), ...
        % % % %         'liability', try_mu, 1)
        % % % %      k=k+1
        % % % %     end
        % % % %     figure; plot(h_sib, h_mz, '.');
        % % % %     hold on; plot(0:0.01:1, 0:0.01:1, 'r');
        % % % %     xlabel('h sib'); ylabel('h mz');
    end
    h_from_ACE_vec = 2*h_from_mz_vec - h_from_sibs_vec; % compute heritability based on ACE model 
    
    
    inferred_parameter = zeros(8, N);
    for i=1:8 % infer a parameter from two others
        for j=1:N
            switch i % which variable are we trying to infer
                case 1 % infer lambda_s from lambda_mz, prevalence
                    inferred_parameter(i,j) = liab_assum_lambda_s_vec(j);
                case 2
                    inferred_parameter(i,j) = ...
                        100 * (liab_assum_lambda_s_vec(j) - lambda_s_vec(j)) / lambda_s_vec(j);
                case 3 % infer lambda_mz from lambda_s, prevalence
                    inferred_parameter(i,j) = liab_assum_lambda_mz_vec(j);
                case 4
                    inferred_parameter(i,j) = ...
                        100 * (liab_assum_lambda_mz_vec(j) - lambda_mz_vec(j)) / lambda_mz_vec(j);
                case 5 % infer prevalence from lambda_s, lambda_mz
                    inferred_parameter(i,j) = 0;
                case 6 % heritability inferred from sib risk  (assume no shared environment)
                    inferred_parameter(i,j) = 100*h_from_sibs_vec(j);
                case 7 % heritability inferred from mz risk  (assume no shared environment)
                    inferred_parameter(i,j) = 100*h_from_mz_vec(j);
                case 8 % heritability inferred from mz and sib risk  using the ACE model
                    inferred_parameter(i,j) = 100*min(1,h_from_ACE_vec(j));
            end % switch j
        end
    end
    lambda_s_ratio = 100*(liab_assum_lambda_s_vec ./ lambda_s_vec - 1); % ratio in percentage
    lambda_s_overestimation_vec = (lambda_s_vec - liab_assum_lambda_s_vec) ./ liab_assum_lambda_s_vec
    
    save('wray_disease_lambda_list.mat', ...
        'disease_vec',  'disease_short_vec', ...
        'lambda_s_vec', 'lambda_mz_vec', 'prevalence_vec', ...
        'liab_assum_lambda_s_vec', 'liab_assum_lambda_mz_vec', ...
        'lambda_s_overestimation_vec', 'inferred_parameter', ...
        'h_from_sibs_vec', 'h_from_mz_vec'); %  'lambda_s_overestimation_vec');
else
    load('wray_disease_lambda_list.mat');     
end % compute parameters
R = [disease_vec' num2cell([lambda_mz_vec' lambda_s_vec' (prevalence_vec.*100)' ...
    inferred_parameter' ])]; % lambda_s_ratio'
for i=1:size(R,1)
    for j=2:size(R,2)
        R{i,j} = sprintf('%.1f', R{i,j}); 
    end
end
%R = num2str_cell(R, 3);
R = [{'Disease', 'lambda_mz', 'lambda_s', 'prevalence (%)', ...
    'lambda_s (exp.)', 'deviation (%)', 'lambda_mz (exp.)', 'deviation (%)', 'prevalence (liab.) ',  ...
    'h^2 (% liab., lambda_s)', 'h^2 (% liab., lambda_mz)', 'h_{pop}^2 (%)' }' R']'; % , 'lambda_s deviation'
savecellfile(R, fullfile(output_dir, 'wray_disease_lambda_list_tab.txt')); % save in tab-delimited .txt file

R = strrep_cell(R, 'prevalence', '$\mu$');
R = strrep_cell(R, '%', '\%');
R = strrep_cell(R, 'h^2', '$h^2$');
R = strrep_cell(R, 'h_{pop}^2', '$h_{pop}^2$');
R = strrep_cell(R, '^{(1)}', '$^{(1)}$');
R = strrep_cell(R, '^{(2)}', '$^{(2)}$');
R = strrep_cell(R, 'lambda_s', '$\lambda_s$');
R = strrep_cell(R, 'lambda_mz', '$\lambda_{MZ}$');
R_latex = latex(R(1:end-3,[1:6 end]), 2, precision); % remove artificial diseases. Take only relevant fields
R_latex = mat2cell(R_latex, ones(size(R_latex,1),1));
R_latex{1} = strrep(R_latex{1},  '|c', '|r'); % align to right  
R_latex{1} = strrep(R_latex{1},  '{|r|', '{|l|'); % align first column to left  

R_latex = [{'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', '{\scriptsize', ...
    '\caption[Epidemiological parameters for disease traits]{Epidemiological parameters for ' ...
    ' a list of diseases. (Source: Wray et al. \cite{Wray:2010ys}) \label{table:wray_diseases}}'} R_latex' ... % reported in [Wray et al.]. 
    {'}',  ...
    '\end{center}', '\end{table}', '', '\newpage', ''}]'; 
% tab_header = {'\begin{table}[h!] % put tables summarizing all architectures', ...
%     '\begin{center}', '{\scriptsize', 'GGGG'};
% tab_footer = {'\end{tabular}', '}', '\caption{Architectures Summary Statistics}', ...
%     '\label{tab:additive_vs_interaction_effects}', ...
%     '\end{center}', '\end{table}', '', '\newpage', ''};
savecellfile(R_latex, fullfile(output_dir, 'wray_disease_lambda_list.txt'), [], 1);  % save also in latex format



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(figs_flag)
    bad_diseases = {'{\bf MLT2}', 'Depression'};
    [~, I, J] = intersect(disease_vec, bad_diseases);
    I = setdiff(1:length(disease_vec), I); % keep these inds
    lambda_s_vec = lambda_s_vec(I); lambda_mz_vec = lambda_mz_vec(I);
    liab_assum_lambda_s_vec = liab_assum_lambda_s_vec(I); disease_vec = disease_vec(I);
    
    full_figure(0); % don't use hold on.
    plot_wray_lambda_table(lambda_s_vec, lambda_mz_vec, liab_assum_lambda_s_vec, disease_vec, 100);
    my_saveas(gcf, fullfile(figs_dir, 'lambda_s_vs_lambda_mz_various_diseases_log_scale'), format_fig_vec);
    % Add a blow-up panel for low values of lambda
    % hp = uipanel('Title','Main Panel','FontSize',12,...
    %     'BackgroundColor','white',...
    %     'Position',[.6 .12 .3 .4]);
    % figure(hp); % set current figure to panel
    
    full_figure(0);
    plot_wray_lambda_table(lambda_s_vec, lambda_mz_vec, liab_assum_lambda_s_vec, disease_vec, 5); % plot inside panel
    %xlim([0 25]); ylim([0 5]);
    
    my_saveas(gcf, fullfile(figs_dir, 'lambda_s_vs_lambda_mz_various_diseases_log_scale_zoom_in'), format_fig_vec);
    
    mean(lambda_s_overestimation_vec)
    median(lambda_s_overestimation_vec)
    std(lambda_s_overestimation_vec)    
end
