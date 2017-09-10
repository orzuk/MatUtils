% Master script for parsing and computing statistics for SNPs from GWAS

% % % 1. The recurrence formula seems to parse incorrectly
% % % 
% % % 2. Unless I'm missing something, this seems straightforward. After $k$ tosses, 
% % % the probability that at least one coin didn't get 'heads', 
% % % is $g_k(n) = 1 - (1-\frac{1}{2^k})^n$. Therefore, $P_k(n) \equiv Prob(T(n) = k) = 
% % % g_{k-1}(n) - g_{k}(n) = (1-\frac{1}{2^{k}})^n - (1-\frac{1}{2^{k-1}})^n$. 
% % % From this, the expectation is $f(n) = 
% % % \sum_{k=1}^{\infty} k [(1-\frac{1}{2^{k}})^n - (1-\frac{1}{2^{k-1}})^n] = 
% % % -\sum_{k=1}^{\infty} (1-\frac{1}{2^{k}})^n$.
% % % 
% % % N = 50000; 
% % % 
% % % f = sum((1:N).* ((1-1./2.^(1:N)).^N - (1-1./2.^(0:(N-1))).^N ))
% % % f2 = 1+sum(1-((1-1./2.^(1:N)).^N))



% Read gwas database from http://genome.gov/gwastudies/
% function read_gwas_database(gwas_database_name)
% Note: for some studies we need to filter fields like odds-ratio

AssignGeneralConstants;
AssignGeneticArchitectureConstants;
[machine machine_delim html_outdir] = get_machine_type();
format_fig_vec = {'fig', 'epsc', 'jpg'};
controls_flag = 'population'; % assume general population ( no correction for computed odds-ratios)
disease_html_root_dir = fullfile(html_outdir, 'common_disease_broad_catalog');

run_database_computation = 1; % determine if to run power computations
run_from_filter = 1; % run filtering and compute additional parameters for data
output_figs_dir = fullfile(strdiff(html_outdir(1:end-1), 'data'), 'genetic_architecture/figures');

% disease_data_from_literature_file = '../../common_disease_model/data/disease_curated_literature_statistics.tsv';
disease_data_from_literature_file = '../../common_disease_model/data/common_disease_epidemiology_statistics.tsv';
% gwas_database_name  = '../../common_disease_model/data/GWASCatalog112608.txt'
% gwas_database_name  = '../../common_disease_model/data/gwas_database_catalog/gwascatalog_2010_04_05.txt'
% gwas_database_name  = '../../common_disease_model/data/gwas_database_catalog/gwascatalog_2010_07_21.txt'

%gwas_NHGRI_database_name  = '../../common_disease_model/data/gwas_database_catalog/gwascatalog_2011_02_21.txt'
gwas_NHGRI_database_name  = '../../common_disease_model/data/gwas_database_catalog/gwascatalog_2010_08_19.txt'





switch catalog_type
    case 'NHGRI'
        gwas_database_dir = '../../common_disease_model/data/gwas_database_catalog/';
        gwas_database_name  = gwas_NHGRI_database_name; % '../../common_disease_model/data/gwas_database_catalog/gwascatalog_2010_08_19.txt'
    case 'MIT-TURK'
        gwas_database_dir = '../../common_disease_model/data/gwas_mit_turk_catalog/';
        %        gwas_database_name  = '../../common_disease_model/data/gwas_mit_turk_catalog/GWAS_catalog_parameters_2010_08_26.txt';
        gwas_database_name  = '../../common_disease_model/data/gwas_mit_turk_catalog/GWAS_catalog_parameters_2011_08_02.txt';
end
unite_traits_flag = 1; % unite different studies of the same disease

if(run_database_computation) % here run all the stuff
    if(~exist(file_name_to_mat(gwas_database_name), 'file'))
        data = ReadDataFile(gwas_database_name, file_name_to_mat(gwas_database_name), -1); % don't convert cell to mat
    else
        data = load(file_name_to_mat(gwas_database_name));
    end
    
    
    switch catalog_type
        case 'NHGRI'
            unite_traits_flag = 1; % unite different studies of the same disease
            
            num_snps = length(data.SNPs)
            
            data.Trait = data.Disease_Trait;
            data.P_value = data.p_Value;
            data.Gene = data.Reported_Gene_s_;
            if(iscell(data.Initial_Sample_Size)) %            data.Initial_Sample_Size% Correct sample size
                data.discovery_num_cases = zeros(num_snps,1);
                data.discovery_num_controls = zeros(num_snps,1);
                data.replication_num_cases = zeros(num_snps,1);
                data.replication_num_controls = zeros(num_snps,1);
                
                for i=1:num_snps
                    tmp_sample_size = str2nums(data.Initial_Sample_Size{i});
                    if(~isempty(tmp_sample_size))
                        data.discovery_num_cases(i) = tmp_sample_size(1);
                    else
                        data.discovery_num_cases(i) = -1;
                    end
                    if(length(tmp_sample_size) > 1)
                        data.discovery_num_controls(i) = tmp_sample_size(2);
                    else
                        data.discovery_num_controls(i) = -1;
                    end
                end
            end
            
        case 'MIT-TURK'         % Figure out how much was already filled by the kids
            unite_traits_flag = 1; % only one study for each disease
            
            data.SNPs = data.SNP_ID_;
            num_snps = length(data.SNPs)
            data.trait_type = data.Type__binary_QTL_;
            data.trait_type_num = zeros(num_snps,1);
            data.trait_type_num(strmatch('QTL', empty_cell_to_empty_str(data.trait_type))) = 1; % get quantitative traits
            cur_trait_name = ''; cur_trait_ind = 0; % Fill trait's name and other variables
            for i=1:num_snps
                if( (isempty(data.Trait{i})) || ...
                        ( (cur_trait_ind > 0) && strcmp(data.Trait{i}, data.Trait{cur_trait_ind}) ) )% compare to previous trait
                    data.Trait{i} = cur_trait_name;
                    data.Effect_Size_Units{i} = data.Effect_Size_Units{cur_trait_ind};
                    data.trait_type{i} = data.trait_type{cur_trait_ind};
                    data.trait_type_num(i) = data.trait_type_num(cur_trait_ind);
                    data.Prevalence{i} = data.Prevalence{cur_trait_ind};
                    data.Reference__lambda_s{i} = data.Reference__lambda_s{cur_trait_ind};
                    data.Discovery_Sample_Size__Cases_(i) = data.Discovery_Sample_Size__Cases_(cur_trait_ind);
                    data.Discovery_Sample_Size__Controls_(i) = data.Discovery_Sample_Size__Controls_(cur_trait_ind);
                    data.Replication_Sample_Size__Cases_(i) = data.Replication_Sample_Size__Cases_(cur_trait_ind);
                    data.Replication_Sample_Size__Controls_(i) = data.Replication_Sample_Size__Controls_(cur_trait_ind);
                    
                    data.Person_annotating{i} = data.Person_annotating{cur_trait_ind};
                else
                    cur_trait_name = data.Trait{i};
                    cur_trait_ind = i; % set index
                end
            end
            data.Trait = strrep_cell(data.Trait, '"', '');
            
            
            %        data.Region = empty_cell_to_empty_str(cell(num_snps,1));
            data.discovery_num_cases = zeros(num_snps,1);
            data.discovery_num_controls = zeros(num_snps,1);
            data.replication_num_cases = zeros(num_snps,1);
            data.replication_num_controls = zeros(num_snps,1);
            data.PUBMEDID = data.Reference__Study_PubMed_ID;
            data.h_ref = data.Reference__Heritability;
            data.h_scale = data.Heritability_Scale__broad_narrow_liability_QTL_;
            data.lambda_s_ref = data.Reference__lambda_s;
            data = get_additional_NHGRI_parameters(data, ...
                file_name_to_mat(gwas_database_name), ...
                file_name_to_mat(gwas_NHGRI_database_name)); % intersect with NHGRI catalog to get misssing data
            for i=1:num_snps
                tmp_sample_size = str2nums(data.Discovery_Sample_Size__Cases_{i});
                if(~isempty(tmp_sample_size))
                    data.discovery_num_cases(i) = sum(tmp_sample_size);
                else
                    data.discovery_num_cases(i) = -1;
                end
                tmp_sample_size = str2nums(data.Discovery_Sample_Size__Controls_{i});
                if(~isempty(tmp_sample_size))
                    data.discovery_num_controls(i) = sum(tmp_sample_size);
                else
                    data.discovery_num_controls(i) = -1;
                end
                tmp_sample_size = str2nums(data.Replication_Sample_Size__Cases_{i});
                if(~isempty(tmp_sample_size))
                    data.replication_num_cases(i) = sum(tmp_sample_size);
                else
                    data.replication_num_cases(i) = -1;
                end
                tmp_sample_size = str2nums(data.Replication_Sample_Size__Controls_{i});
                if(~isempty(tmp_sample_size))
                    data.replication_num_controls(i) = sum(tmp_sample_size);
                else
                    data.replication_num_controls(i) = -1;
                end
            end
            
            
            % TEMP WRONG HACK FOR T2D:
            % ---------------------------
            %         xx = strmatch('Type 2 diabetes', empty_cell_to_empty_str(data.Trait), 'exact')
            %         for i=27:length(xx)
            %             data.Trait{xx(i)} = ['XXX_Type 2 diabetes_X'];
            %         end
            
            trait_inds = 1;
            for i=2:length(data.Trait)
                if(~strcmp(data.Trait{i}, data.Trait{i-1}))
                    if(~isempty(data.Trait{i}))
                        trait_inds = [trait_inds i];
                    end
                end
            end
            studies_done_inds = trait_inds(find(~isempty_cell(data.SNP_ID_(trait_inds))));
            num_studies_done = length(studies_done_inds); % sum(~isempty_cell(data.SNP_ID_(trait_inds)))
            special_comments_inds = find(~isempty_cell(data.Comments));
            num_additional_comments_done = length(setdiff(special_comments_inds, studies_done_inds));
            
            sibling_relative_risk_inds = find(~isempty_cell(data.Sibling_Relative_Risk__lambda_s_(trait_inds)));
            heritability_inds = find(~isempty_cell(data.Heritability(trait_inds)));
            prevalence_inds = find(~isempty_cell(data.Prevalence(trait_inds)));
            num_epidemiology_done = length(union(prevalence_inds, union(sibling_relative_risk_inds, heritability_inds)));
            num_loci_done = sum((~isempty_cell(my_unique(data.SNP_ID_))));
            
            if(~isfield(data, 'Risk_Allele_Frequency'))
                data.Risk_Allele_Frequency = data.Risk_Allele_Freq_;
            end
            data.OR_or_beta_discovery = data.Discovery_Effect_Size;
            data.OR_or_beta_replication = data.Replication_Effect_Size;
            data.OR_or_beta_combined = data.Combined_Effect_Size;
    end % switch catalog type
    
    save([gwas_database_name(1:end-4) '_processed.mat'], 'data');
    
end % run database computations

load([gwas_database_name(1:end-4) '_processed.mat'], 'data');

if(run_from_filter) % next stage
    output_gwas_database_name = [remove_suffix_from_file_name(gwas_database_name) '_complete.txt'];
    [data filtered_inds] = filter_gwas_data(data, output_gwas_database_name, catalog_type); num_filtered = length(filtered_inds);
    
    filter_type = 'keep_numeric';
    switch filter_type
        case 'keep_numeric'
            
        case 'only_binary'
            data = struct_by_inds(data, find(1-data.trait_type_num)); % take only quantitative traits - why?
    end
    num_filtered =length(data.Trait);
    
    % Fomrat: disease name, prevalence, lambda_s, heritability ???
    
    switch catalog_type % add epidemiological paraemters from other files
        case 'NHGRI'
            data = add_disease_familial_data(data, disease_data_from_literature_file); % add missing epidemiological variables
        case 'MIT-TURK' % shouldn't be done. All data should be in the google docs !!!
            base_prevalence = 0.01; base_lambda_s_familial = 2; base_h_familial = 1;
            data.Prevalence = str2nums_cell(data.Prevalence, 1, 1); % convert percentiles to nums
            data.base_prevalence_inds = isempty_cell(data.Prevalence);
            for i=vec2row(find(data.base_prevalence_inds))
                data.Prevalence{i} = base_prevalence;
            end
            data.Prevalence = cell2mat(data.Prevalence);
            
            for i=1:length(data.Prevalence)
                switch data.trait_type{i}
                    case {'binary', 'Binary'}
                        data.OR(i,:) = min(data.OR(i,:), 1 ./ data.Prevalence(i)); % avoid too high values !!! (error, QTL!!!)
                end
            end
            
            data.lambda_s_familial = str2nums_cell(data.Sibling_Relative_Risk__lambda_s_, 1, 1);
            data.base_lambda_s_inds = isempty_cell(data.lambda_s_familial); % find which have a sibling relative risk
            for i=1:length(data.lambda_s_familial)
                if(isempty(data.lambda_s_familial{i}))
                    data.lambda_s_familial{i} = base_lambda_s_familial;
                end
            end
            data.lambda_s_familial = cell2mat(data.lambda_s_familial);
            
            data.h_familial = str2nums_cell(data.Heritability, 1, 1);
            for i=1:length(data.h_familial)
                if(isempty(data.h_familial{i}))
                    data.h_familial{i} = base_h_familial;
                end
            end
            data.h_familial = cell2mat(data.h_familial);
    end
    switch controls_flag
        case 'non_cases' % here we need to correct odds ratio to genetic relative risk, and also control to population allele frequency
            [data.GRR data.RAF] = odds_ratio_to_genetic_relative_risk(data.OR, data.RAF, data.Prevalence);
        case 'population'
            data.GRR = data.OR; % data.MAF = data.RAF;
    end
    
    
    
    
    if(unite_traits_flag)
        data = unite_gwas_data_by_traits(data); field = 'Trait';
    else
        field = 'PUBMEDID';
    end
%     [data_params.lambda_s data_params.lambda_mz data_params.lambda_mz_add data_params.h_add data_params.V_add ...
%         data_params.h data_params.V data_params.unique_inds data_params.num_loci ...
%         data_params.max_penetrance data_params.h_liab data_params.Prevalence ...
%         data_params.binary_inds data_params.trait_name ...
%         data_params.Power data_params.Power_empirical data_params.non_centrality ...
%         data_params.std data_params.OR_min_confidence data_params.OR_max_confidence ... % New: add confidence parameters
%         data_params.snp_lambda_s data_params.snp_h_liab data_params.snp_h_liab_min data_params.snp_h_liab_max ...
%         data_params.effective_sample_size data_params.alpha_vec ...
%         data_params.beta_discovery_grid, data_params.beta_discovery_inv_grid] = ...
        
    data_params = compute_gwas_parameters(data, field);  % Compute explained heritability for each disease
    data_params.num_snps = length(data.SNPs);
    data_params.quant_inds = setdiff(1:data_params.num_snps, data_params.binary_inds);
    %    [data.MAF(data_params.binary_inds), data.OR_Directon(data_params.binary_inds)] = ... % problem: this doesn't flip the OR!!!
    %        flip_allele(data.MAF(data_params.binary_inds), data.OR(data_params.binary_inds), 'Binary', 1); % Set MAF always less than 0.5
    %    [data.MAF(data_params.quant_inds), ~] = ...
    %        flip_allele(data.MAF(data_params.quant_inds), data.OR(data_params.quant_inds), 'QTL', 1); % Set MAF always less than 0.5
    
    plot_data = 0;
    if(plot_data)
        plot_gwas_data_stats(data, data_params, output_figs_dir); % generate several plots for data
    end
    
    % Here save to latex
    save_latex = 0;
    if(save_latex)
        labels = {'disease', 'prev. (\%)', 'num loci', '$h^2 (\%)$', '${H^2}_{mult} (\%)$', ...
            '$H^2$ (fam.) (\%)', 'h-explained (\%)', '$\lambda_s$', ...
            '$\lambda_s$ (fam.)', '$\lambda_s$ explained (\%)'}; %   'max penet.'}; % don't take study
        h_familial = data_params.h; h_familial(:) = -1; % we don't know the h familial (need to fill)
        lambda_s_familial = h_familial;
        h_explained = max(0, max(data_params.h_add, data_params.h) ./ data.h_familial(data_params.unique_inds));
        lambda_s_explained = max(0, (data_params.lambda_s - 1) ./ (data.lambda_s_familial(data_params.unique_inds) - 1));
        R = [data.Trait(data_params.unique_inds) ...
            num2cell([data.Prevalence(data_params.unique_inds)*100 data_params.num_loci data_params.h_add*100 ...
            data_params.h*100 data.h_familial(data_params.unique_inds)*100 ... %  h_familial ...
            h_explained.*100 data_params.lambda_s], precision) num2str_cell(data.lambda_s_familial(data_params.unique_inds), precision, '-') ...
            num2cell([lambda_s_explained.*100])]; %  max_penetrance'])];
        R = [labels' sortrows(R)']';
        gwas_stats_table_file = ...
            '../../common_disease_model/data/gwas_database_catalog/gwas_stats_table.txt';
        for i=2:size(R,1)
            for j=2:size(R,2)
                R{i,j} = ['{\bf ' num2str(R{i,j}, precision) '}'];
            end
        end
        savecellfile(R, gwas_stats_table_file, [], 1);
        
        
        [intersect_vec I J] = intersect(lower(R(:,1)), lower(disease_data_from_literature_vec.Disease));
        for i=2:size(R,1)
            R{i,1} =  ['{\bf ' R{i,1} '}'];
        end
        
        R_small = R([1 I'],:);
        latex_R = latex(R_small, 2, precision);
        latex_R = mat2cell(latex_R, ones(size(latex_R,1),1));
        savecellfile(latex_R, file_name_to_new_suffix(gwas_stats_table_file, 'tex'), [], 1); % save also in latex
    end
    
    % save_expression_mat_as_text_generic(labels(1:4), R, [], gwas_stats_table_file, 1);
    
    
    % Temp for debugging:
    % f = 0.7; GRR = 4; mu = 0.05;
    % [lambda_s_vec lambda_s lambda_mz h_add V_add h] = ...
    %     genetic_relative_risk_to_heritability(f, GRR, mu)
    % beta = mu / (1+f*(GRR-1));
    % alpha = beta * (GRR-1);
    % hhh = alpha^2 * f * (1-f) / ( (beta+alpha*f) * (1-beta-alpha*f) )
    
    disease_html_template_file = fullfile(disease_html_root_dir, 'common_diseases_template.html');
    
    save(fullfile(disease_html_root_dir, 'disease_data.mat'), ...
        'data', 'data_params', 'disease_html_root_dir', 'disease_html_template_file');
    
end % if run database computation

load(fullfile(disease_html_root_dir, 'disease_data.mat'));
disease_html_root_dir = fullfile(html_outdir, 'common_disease_broad_catalog');
disease_html_template_file = fullfile(disease_html_root_dir, 'common_diseases_template.html');
write_disease_tables_html(data, data_params, disease_html_root_dir, disease_html_template_file, 0); % save everything to html


% test_maximize_allele_freq_likelihood
% % Tmp: determine effect of HLA locus in T1D
% GRR = 8.03;
% [lambda_s_vec ...
%     lambda_s_add lambda_mz_add h_add V_add ...
%     lambda_s_mult lambda_mz_mult] = ...
%     genetic_relative_risk_to_heritability([0.1 0.1]', [GRR GRR]', 0.004)
%


