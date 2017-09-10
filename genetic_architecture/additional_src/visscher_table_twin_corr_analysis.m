% Read and analyse data from r_MZ and r_DZ table from Visscher
% function visscher_table_twin_corr_analysis()
AssignGeneralConstants;
plot_figs = 0;
visscher_table_file = '../../common_disease_model/data/epidemiology/visscher_QTL_twin_correlation_statistics.txt';
output_dir = '../../common_disease_model/docs/pnas/genetic_interactions/tables/';

R = ReadDataFile(visscher_table_file,[],[],1);

R = struct_by_inds(R, 2:87); % take 86 traits
R.all_r_DZ = cell2mat(R.all_r_DZ);
R.all_r_MZ = cell2mat(R.all_r_MZ);
R.rMZ_2rDZ = cell2mat(R.rMZ_2rDZ);
R.V_r_DZ = cell2mat(R.V_r_DZ);
R.V_r_MZ = cell2mat(R.V_r_MZ);
R.V_rMZ_2rDZ = cell2mat(R.V_rMZ_2rDZ);
R.SE = cell2mat(R.SE);
R.Phenotype = strrep_cell(R.Phenotype, '"', '');


S = [R.Phenotype  num2cell(100*[R.all_r_MZ R.all_r_DZ R.rMZ_2rDZ])]; % Save data in latex
[sorted_phenotypes sort_perm] = sort(R.Phenotype); 
S = S(sort_perm,:); % sort alphabetically 
S = [{'Trait', '$r_{MZ} (\%)$', '$r_{DZ}  (\%)$', '$r_{MZ}-2r_{DZ}  (\%)$'}' S']';
savecellfile(S, fullfile(output_dir, 'hill_quantitative_traits_list_tab.txt'), [], 1);

%S = strrep_cell(S, 'Serum', 'Serum'); 
S = strrep_cell(S, 'Alcohol consumption (combined model)', 'Alcohol consumption ');
S_latex = latex(S, 2, precision);
S_latex = mat2cell(S_latex, ones(size(S_latex,1),1));
S_latex{1} = strrep(S_latex{1}, '|c|c|c|c|', '|p{5.0cm}|r|r|r|'); 

%S = [S(1:44,:) S([1 45:end],:)]; % make two columns 
% S_latex{1} = strrep(S_latex{1}, '|c|c|c|c|c|c|c|c|', '|p{3.6cm}|r|r|r|p{3.6cm}|r|r|r|'); 

S_latex = strrep_cell(S_latex, 'Body Mass Index', 'BMI');

tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', ...
    ['\caption[Epidemiological parameters for quantitative traits]{Epidemiological parameters for ' ...
    ' a list of quantitative traits. (Source: Hill et al. \cite{hill2008data})}'], ...  %  reported in [Hill et al.].   
    '{\scriptsize', S_latex{1}};
tab_footer = {'}',  ...
    '\label{table:hill_quantitative_traits}', ...
    '\end{center}', '\end{table}', '', '\newpage', ''};
%S_latex = [tab_header  S_latex(2:end)' tab_footer]'; 
S_latex = split_latex_table(S_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 45, -1);

savecellfile(S_latex, fullfile(output_dir, 'hill_quantitative_traits_list.txt'), [], 1);  % save also in latex format



if(plot_figs)
    show_traits = {'Body Mass Index', 'Height',  'Height (clinically measured)', ...
        'Lipid',  'Total cholesterol', 'ApoE',   'Birth weight', ...
        'Age at menarche',   'Body Mass Index, age 20-29', ... % 'Age at menopause',
        'Body Mass Index, age 30-39',  'Systolic blood pressure', 'Diastolic blood pressure'};
    [~, show_inds] = intersect(R.Phenotype, show_traits);
    
    full_figure; plot(R.all_r_DZ, R.all_r_MZ./2, '.', 'linewidth', 2);
    plot(R.all_r_DZ(show_inds), R.all_r_MZ(show_inds)./2, '*k', 'linewidth', 2);
    plot(0:0.01:0.8, 0:0.01:0.8, 'r', 'linewidth', 2); xlabel('r_{DZ}'); ylabel('r_{MZ}/2');
    title('Comparison of Monozgous and Dyzygous twin correlations for 86 QTLs');
    % for i=show_inds
    %    text(R.all_r_DZ(i)+0.01, R.all_r_MZ(i)/2,  R.Phenotype{i}, 'fontweight', 'bold');
    % end
    %title('Deviation from additivity : r_{MZ} - 2r_{DZ} mean \pm st.d.');
    title('Deviation from additivity : r_{MZ}/2 vs. r_{DZ}');
    
    my_saveas(gcf, '../../common_disease_model/figs/r_s_vs_r_mz_various_traits_main', format_fig_vec);
    ylabel('r_{MZ} - 2r_{DZ}');
    
    [s sort_perm] = sort(R.rMZ_2rDZ); % Plot differences with error bars
    full_figure; errorbar(s, R.SE(sort_perm), 'linewidth', 2); plot(s, 'r', 'linewidth', 2);
    hold on; plot(1:86, zeros(1,86), 'k--');
    ylabel('r_{MZ} - 2 r_{DZ}');
    my_saveas(gcf, '../../common_disease_model/figs/r_mz_minus_2_r_dz_various_traits', format_fig_vec);
    
    
    
    for i=1:86
        [s_mu(i) s_sigma(i)] = rationormstat(R.all_r_DZ(i), sqrt(R.V_r_DZ(i)), ...
            R.all_r_MZ(i), sqrt(R.V_r_MZ(i)));
    end
    [s sort_perm] = sort(2.*R.all_r_DZ ./ R.all_r_MZ); % Plot ratio's with error bars
    full_figure; errorbar(1-s, s_sigma, 'linewidth', 2); plot(1-s, 'r', 'linewidth', 2);
    hold on; plot(1:86, zeros(1,86), 'k--');
    ylabel('1 - 2 r_{DZ}/r_{MZ}');
    title('Relative deviation in r_{DZ} from additivity');
    my_saveas(gcf, '../../common_disease_model/figs/r_mz_vs_2_r_dz_ratio_various_traits', format_fig_vec);
    
    figure; hist(100*(s-1),20); title(['Relative deviation in r_{DZ} from additivity. \delta = (2r_{DZ}-r_{MZ})/r_{MZ} \mu=' ...
        num2str(mean(s-1)*100,3) '%, median=' num2str(median(s-1)*100,3) '%, \sigma=' num2str(std(s-1)*100,3) '%']);
    xlabel('(2 r_{DZ} - r_{MZ})/r_{MZ} (%)'); ylabel('# Traits');
    my_saveas(gcf, '../../common_disease_model/figs/r_mz_vs_2_r_dz_ratio_various_traits_distribution', ...
        format_fig_vec);
end


