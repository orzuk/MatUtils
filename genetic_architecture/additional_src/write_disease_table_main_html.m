% Internal function: Save one page for all diseases
%
% Input:
% data - structure containing SNP-specific date
% disease_data - additional disease-specific data
% data_params - more data parameters
% disease_names - names of traits
% num_diseases - number of traits
% catalog_html_file - where to save catalog
% root_disease_dir - where to save figures
% good_familial_inds_vec - which traits have epidemiological parameters
% R_start - template html structure start
% R_end - template html structure end
%
% Output:
% S_main_table - table with trait-specific parameters (html format)
% S_main_table_cell - table with trait-specific parameters (cell format)
% disease_table_file - file where ..
% h_explained_total - heritability explained for each snp (?)
% lambda_s_explained_total -
% num_shadow_snps_vec - number of unidentified SNPs for each SNP found
% power_correct - vector with power corrections
% disease_power_inds - indices of correct power adjustment
%
function [S_main_table S_main_table_cell disease_table_file ...
    h_explained_total lambda_s_explained_total num_shadow_snps_vec power_correct ...
    disease_power_inds] = ...
    write_disease_table_main_html(data, disease_data, data_params, ...
    disease_names, num_diseases, catalog_html_file, ...
    root_disease_dir, good_familial_inds_vec, R_start, R_end, plot_main_figures_flag) % save the main website

AssignGeneralConstants;
machine = get_machine_type();
AssignGeneticArchitectureConstants();

correction_mode = 'exact'; % NEW TEMP
switch machine
    case UNIX
        bmp_fig_format = {'jpg'};
    case PC
        bmp_fig_format = {'bmp', 'jpg'};
end
if(~iscell(correction_type))
    correction_type = {correction_type};
end
num_corrections = length(correction_type); % enable more than one type of corrections
for i=1:num_corrections
    correction_str{i} = ['h<sup>2</sup> explained (pow. adj. ' correction_type{i} ')'];
end
S = cell(num_diseases+2, 16+num_corrections); % prepare table
S(1,:) = [{'#', 'Trait', 'Type', 'Prevalence', 'Heritability', ... % 1-5
    'Sib-Risk(&lambda;<sub>s</sub>)', 'Study', ... % 6-10 % removed 6, 7, 9
    '# Cases', '# Controls', 'Effective-Sample-Size', '# Loci', 'Var-Explained', ...
    'h<sup>2</sup>-explained', '&lambda;<sub>s</sub> explained', '# shadow-loci'} ...
    correction_str 'fit'];
% '# shadow-loci', 'h<sup>2</sup> explained (pow. adj.)', 'fit'}; % 11-15

% New: add different corrections


%    'Heritability-Scale', 'Heritability-Ref.', 'Sibling-Relative-Risk', 'Sibling-Relative-Risk-Ref.', 'Study', ... % 6-10

[~, ~, html_outdir] = get_machine_type();
num_cutoffs = length(data_params.alpha_vec);
heritability_total_bar = ones(num_diseases,4);  % make a large bar plot h^2
lambda_s_total_bar = ones(num_diseases,4); % make a large bar plot for lambda_s
all_trait_first_inds = zeros(num_diseases,1);
power_correct = cell(num_diseases, num_corrections);
num_shadow_snps_vec = cell(num_diseases, num_corrections);
disease_power_inds = zeros(num_diseases, 1); disease_table_file = cell(num_diseases,1);
h_explained_total = zeros(num_diseases, 1); lambda_s_explained_total = zeros(num_diseases, 1);
for i=1:num_diseases % loop on diseases
    run_i = i
    run_disease = disease_names{i}
    trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait
    all_trait_first_inds(i) = trait_inds(1);
    disease_data.num_loci(i) = length(trait_inds);
    disease_table_file{i} = fullfile(disease_names{i}, [disease_names{i} '_table.html']);
    h_explained_total(i) = sum(data_params.snp_h_liab(trait_inds)); % compute variance explained statistics
    lambda_s_explained_total(i) = prod(data_params.snp_lambda_s(trait_inds));
    cur_alpha_used = data.disease_alpha_vec(trait_inds(1));
    [~, cur_ind_used] = min(abs(data_params.alpha_vec - cur_alpha_used)); % choose with what cutoff to work for this disease
    disease_power_inds(i) = cur_ind_used;
    if(length(cur_ind_used) == 1) % didn't specify index for each correction
        cur_ind_used = [cur_ind_used ones(1,num_corrections-1)];
    end
    for j=1:num_corrections % loop on different corrections
        switch correction_type{j}
            case 'empirical'
                power_correct{i,j} = data_params.Power_empirical(trait_inds,:);
            case 'theoretical'
                power_correct{i,j} = data_params.Power(trait_inds,:);
            case 'pop-gen' % what do we put here ???
                if( (length(data_params.VarExplainedCorrected) >= i) && ~isempty(data_params.VarExplainedCorrected{i})) % check we actually computed sophisticated correction
                    power_correct{i,j} = data_params.VarExplained{i} ./ ...
                        data_params.VarExplainedCorrected{i}(:,correction_inds(j));% determine power as ratio of vars explained
                else % we didn't compute the special correction here - just set to one
                    power_correct{i,j} = ones(length(trait_inds), 1);
                end
        end
        switch correction_mode
            case 'floor' % use the same based on #shadow loci
                num_shadow_snps_vec{i,j} = floor(1 ./ max(MIN_POWER, power_correct{i,j}))-1; % count how many shadow SNPs are there for each SNP
            case {'round'} % use the same based on #shadow loci
                num_shadow_snps_vec{i,j} = round(1 ./ max(MIN_POWER, power_correct{i,j}))-1; % count how many shadow SNPs are there for each SNP
            case 'exact'
                num_shadow_snps_vec{i,j} = (1 ./ max(MIN_POWER, power_correct{i,j}))-1; % count how many shadow SNPs are there for each SNP
        end
    end
    data_params.h_liab_power_adjusted(:,i) = ... % h_explained_power_adjusted = ...
        100 * sum(repmat(data_params.snp_h_liab(trait_inds), 1, num_cutoffs) ./ ...
        max(MIN_POWER, power_correct{i})); % divide by power. Take minimum: 0.05
    
    h_explained_power_adjusted_vec = cell(num_corrections,1); 
    for j=1:num_corrections
        h_explained_power_adjusted_vec{j} = ...
            repmat(data_params.snp_h_liab(trait_inds), 1, size(num_shadow_snps_vec{i,j},2)) .* ...
            (num_shadow_snps_vec{i,j}+1);  %   max(MIN_POWER, data.Power(trait_inds,:));
    end
    
    %    switch disease_data.trait_type{i} % show sibling relative risk for binary traits
    %        case 'Binary'
    heritability_total_bar(i,2) = str2double(disease_data.h_familial{i})/100; % get the heritability estimate
    heritability_total_bar(i,3) = (data_params.h_liab_power_adjusted(i)/100); % get the inflated heritability estimate
    heritability_total_bar(i,4) =  h_explained_total(i); % data_params.h_liab(i); % why are the two different????
    %             familial_bar = [str2double(disease_data.h_familial{i})/100 1]; % str2double(disease_data.lambda_s_familial{i})]
    %             adjusted_bar = [(data_params.h_liab_power_adjusted(i)/100), ...
    %                 str2double(data_params.lambda_s_power_adjusted_percent{i})/100];
    %             explained_bar = [data_params.h_liab(i), str2double(lambda_s_explained_total_percent)/100]; %  data_params.lambda_s(i)];
    %             x_ticks = {'h^2', '\lambda_s'};
    %         case 'QTL'
    %             total_bar = 1;
    %             familial_bar = str2double(disease_data.h_familial{i})/100;
    %             adjusted_bar = (data_params.h_liab_power_adjusted(i)/100);
    %             explained_bar = data_params.h_liab(i);
    %             x_ticks = {'h^2'};
    %     end
    
    
    
    S{i+1,1} = num2str(i);
    S{i+1,2} = html_write_link(disease_names{i}(1:min(30,length(disease_names{i}))), ...
        disease_table_file{i}); % just get diseaes name (cut long names)
    S{i+1,3} = disease_data.trait_type{i};
    S{i+1,4} = disease_data.Prevalence{i};
    S{i+1,5} = html_write_link([disease_data.h_familial{i} '%'], ... %     [disease_data.h_familial{i} '%'];
        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.h_ref{i}]);
    %    S{i+1,6} = disease_data.h_scale{i};
    %    S{i+1,7} = html_write_link(disease_data.h_ref{i}, ...
    %        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.h_ref{i}]);
    ctr=6;
    S{i+1,ctr} = html_write_link(disease_data.lambda_s_familial{i}, ... % disease_data.lambda_s_familial{i};
        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.lambda_s_ref{i}]); ctr=ctr+1; % lambda_s with link to study
    %    S{i+1,9} = html_write_link(disease_data.lambda_s_ref{i}, ...
    %        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.lambda_s_ref{i}]);
    S{i+1,ctr} = html_write_link(disease_data.PUBMEDID{i}, ...
        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.PUBMEDID{i}]); ctr=ctr+1;
    S{i+1,ctr} = disease_data.discovery_num_cases(i); ctr=ctr+1; % number of cases
    S{i+1,ctr} = disease_data.discovery_num_controls{i}; ctr=ctr+1; % number of controls (for diseases)
    S{i+1,ctr} = num2str(round(data_params.effective_sample_size(trait_inds(1)))); ctr=ctr+1; % get proxy for power
    S{i+1,ctr} = disease_data.num_loci(i); ctr=ctr+1;  % number of loci found
    S{i+1,ctr} = [num2str(h_explained_total(i)*100, 4) '%']; ctr=ctr+1; % total variance explained
    S{i+1,ctr} = [num2str((h_explained_total(i) / ...
        str2double(disease_data.h_familial{i}))*100*100, 4) '%']; ctr=ctr+1; % total h^2 explained
    S{i+1,ctr} = [num2str(lambda_s_explained_total(i), 4) ]; ctr=ctr+1; % lambda_s explained
    S{i+1,ctr} = [num2str(sum(num_shadow_snps_vec{i}(:,cur_ind_used(1)))) ]; ctr=ctr+1; % # shadow loci
    
    for j=1:num_corrections % copy var. explained for different types of corrections. This DOESN'T MATCH!!!!!! 
        S{i+1,ctr+j-1} = [num2str((sum(h_explained_power_adjusted_vec{j}(:,cur_ind_used(j))) / ...
            str2double(disease_data.h_familial{i})) *100*100, 4) '%']; % h^2 explained power adjusted
    end
    ctr=ctr+num_corrections;
    %    fit_string{i} = '$alpha*n<sup>$beta</sup$';
    S{i+1,ctr} = 'fit_string'; % fit_string{i}; % display fit
    
    %    S{i+1,12} = [' <br>']; % Move to next line
end % loop on all traits
S{end,1} = '<b>Total:</b>';    % Create a sum table
for k = (14-3):(14-3) % For now only num loci ...
    S{end,k} = ['<b>' num2str(sum(cell2mat(S(2:end-1,k)))) '</b>'];
end


S = num2str_cell(empty_cell_to_empty_str(S));
S_tab = R_start;
ctr = length(S_tab)+1;
% ctr = ctr+1;
S_main_table_cell = S;
S_main_table = vec2row(html_write_table(S, 1)); % This is the main table (to be used later)
S_tab = [S_tab' S_main_table]'; % add table

S_tab{end+1,1} = ['Explanation on all entries in the table, including units, how measured etc. is ' ...
    'avaliable in the ' html_write_link('README', 'README.txt') ' file <br>'];
S_tab{end+1,1} = ['Some summary plots statistics for gwas data as well as theoretical curves are ' ...
    html_write_link('here', 'common_diseases_summary_figs.html') ' <br>']; % New: add figures below the table
S_tab{end+1,1} = ['A concatenated list of ALL loci is ' ...
    html_write_link('here (QTL)', 'common_diseases_all_loci_concatenated_QTL.html') ...
    ' and ' html_write_link('here (binary)', 'common_diseases_all_loci_concatenated_binary.html') ' <br>']; % New: add figures below the table


ctr = size(S_tab,1); % add stuff at the end to file
for i=1:length(R_end)
    S_tab{ctr+i,1} = R_end{i};
end
% S = [S' R_end']';
savecellfile(S_tab, catalog_html_file); % save everything to main home page


S_img = R_start; % Create webpage for figures
explanation_ind = strfind_cell(S_img, 'List of parameters');
S_img{explanation_ind} = [' <LI> List of theoretical curves for different genetic epidemiological parameters, ' ...
    ' and summary plots for data from the common disease catalog </LI>'];
strrep_cell(S_img, 'Parameters', 'Figures');
S_img{end+1,1} = html_write_img('Theoretical relation between h<sup>2</sup> and &lambda;<sub>s</sub> (under liability-threshold model)', ...
    'figs/h_params_relations.jpg');
S_img{end+1,1} = html_write_img('Theoretical relation between prevalence (&mu) and &lambda;<sub>s</sub> (under liability-threshold model)', ...
    'figs/h_params_relations_prevalence_vs_lambda_s.jpg');
S_img{end+1,1} = html_write_img('Theoretical relation between  &lambda;<sub>s</sub> and &lambda;<sub>MZ</sub>(under liability-threshold model)', ...
    'figs/h_params_relations_lambda_s_vs_lambda_mz.jpg'); % Theoretical plots
S_img{end+1,1} = html_write_img('Effect sizes (GRR or &beta) for all SNPs in all traits ', ...
    'figs/odds_ratio.jpg'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img('Risk Allele Frequency for all SNPs in all traits ', ...
    'figs/risk_allele_freq.jpg'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img('Effect sizes (GRR or &beta) vs. Risk Allele Frequency for all SNPs in all traits ', ...
    'figs/odds_ratio_vs_risk_allele_freq.jpg'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img(['Heritability explained for all traits. <br>' ...
    'Trait sorted by fraction of total variance already explained by known loci (shown in green)'], ...
    'figs/genetic_effect_explained_summary.bmp'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img(['Heritability explained vs. sample size. As expected, studies with larger <br> ' ...
    'sample size achieve a higher fraction of heritability explained, indicating that increasing sample size will <br> ' ...
    'lead to a larger fraction explained for other traits with currently under-powered studies.'], ...
    'figs/var_explained_vs_sample_size.bmp'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img('Variance explained vs. rank  for each SNP (all traits, semilog scale) ', ...
    'figs/variance_explained_vs_rank_summary.bmp'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img('Variance explained vs. rank  for each SNP (all traits, loglog scale) ', ...
    'figs/variance_explained_vs_rank_summary_loglog.bmp'); % Empirical plots for catalog data
S_img{end+1,1} = html_write_img('Power of found loci for chosen traits ', ...
    'figs/power_of_found_loci_vs_rank_summary.bmp'); % Empirical plots for catalog data

ctr = size(S_img,1); % add stuff at the end to file
for i=1:length(R_end)
    S_img{ctr+i,1} = R_end{i};
end
savecellfile(S_img, strrep(catalog_html_file, 'table', 'summary_figs')); % save everything to figures global page



special_traits_flag = 1;
if(special_traits_flag) % make some plots for a few special traits
    [~, special_inds] = intersect(disease_data.Trait, special_traits);
    special_inds = intersect(special_inds, find(disease_data.trait_type_num==0)); % take only binary traits for special plot
    special_traits_again = disease_data.Trait(special_inds);
    if(~isempty(special_inds))
        plot_heritability_parameters(cell2mat(disease_data.Prevalence(special_inds)), ...
            cell2mat(str2num_cell(disease_data.lambda_s_familial(special_inds))), ...
            cell2mat(str2num_cell(disease_data.h_familial(special_inds))) ./ 100, ...
            special_traits_again, ...
            fullfile(html_outdir, 'common_disease_broad_catalog/figs/h_params_relations')); % plot a curve matching heritability and lambda_s for binary traits
    end
    % Plot a bar graph showing how much heritability is explained for each disease
end

heritability_total_bar = min(heritability_total_bar, 1.1); % 1.1 indicates we explain too much
[~, sort_perm] = sort(heritability_total_bar(:,4), 'descend'); % sort by h explained
heritability_total_bar = heritability_total_bar(sort_perm, :); % sort by perm
if(plot_main_figures_flag)
    full_figure;
    bar(heritability_total_bar(:,1), 'k');
    bar(heritability_total_bar(:,2), 'r');
    bar(heritability_total_bar(:,3), 'b');
    bar(heritability_total_bar(:,4), 'g');
    %if(good_familial_inds_vec(i))
    legend_vec = {'Environmental', 'Genetic-unexplained', 'Explained (power adjusted)', 'Explained (known loci)'};
    %else
    %    legend_vec = {'Environmental+Genetic-unexplained', 'Explained (power adjusted)', 'Explained (known loci)'};
    %end
    %set(gca, 'XTick', [1 2]); % set labels
    %set(gca, 'XTickLabel', x_ticks); % set labels
    ylabel('h^2'); % Heritability
    legend(legend_vec);
    for i=1:num_diseases
        text(i, 0, str2title(disease_names{sort_perm(i)}), 'rotation', 90);
    end
    % switch disease_data.trait_type{i} % show sibling relative risk for binary traits
    %     case 'Binary'
    %         lambda_s_ticks = [0 1]; % log(str2double(disease_data.lambda_s_familial{i}))];
    %         %                lambda_s_ticks = [1:(lambda_s_ticks-1)/10: lambda_s_ticks;
    %         add_right_yticks(lambda_s_ticks, '\lambda_s'); % right-side scale is for sibling relative risk
    %         ylabels = get(gca, 'YTickLabel');
    %         tmp_lambda_familial = str2double(disease_data.lambda_s_familial{i});
    %         if(isempty(tmp_lambda_familial))
    %             tmp_lambda_familial = data_params.lambda_s_power_adjusted(i) % temp wrong! (only if no sibling relative risk available!)
    %         end
    %         ylabels = num2str(tmp_lambda_familial.^str2double(ylabels), 3);
    %         set(gca, 'YTickLabel', ylabels);
    %         %                my_loglog('y'); % set as exponential scale
    % end
    my_saveas(gcf, fullfile(root_disease_dir, 'figs', 'genetic_effect_explained_summary'), bmp_fig_format); % summary of heritability explained for all diseases
    
    full_figure; plot(data_params.effective_sample_size(all_trait_first_inds(sort_perm)), 100*heritability_total_bar(:,4), '.');
    title('Fraction of Variance Explained vs. Study Power');
    xlabel('Effective Sample Size'); ylabel('Frac. Var. Explained (%)');
    my_saveas(gcf, fullfile(root_disease_dir, 'figs', 'var_explained_vs_sample_size'), bmp_fig_format); % summary of heritability explained for all diseases
end

