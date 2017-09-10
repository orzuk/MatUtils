% Make a central web-page for all diseases and a disease-specific page for
% each one
%
% Input:
% data - a structure with stuff on each diseasese
% data_params - an additional structure with data on each diseasese
% root_disease_dir - where to save the html files
% html_template_file - template to start html files from
% disease_specific_plots - whether to create plots for each disease
%
function write_disease_tables_html(data, data_params, ...
    root_disease_dir, html_template_file, disease_specific_plots)

plot_main_figures_flag = 1; % Temp: don't plot ..
AssignGeneralConstants;
machine = get_machine_type();
switch machine
    case UNIX
        bmp_fig_format = {'jpg'};
    case PC
        bmp_fig_format = {'bmp', 'jpg'};
end
AssignGeneticArchitectureConstants;
if(~iscell(correction_type))
    correction_type = {correction_type};
end
num_corrections = length(correction_type); % enable more than one type of corrections
correction_str = cell(1, num_corrections);
for i=1:num_corrections
    correction_str{i} = ['h<sup>2</sup> explained (pow. adj. ' correction_type{i} ')'];
end


if(~exist('disease_specific_plots', 'var') || isempty(disease_specific_plots))
    disease_specific_plots = 1; % 1; % make also plots (otherwise just link to them)
end
curators_list = {'Diana', 'Diana J Zhu'; ...
    'Carlos', 'Carlos A Castellanos'; ...
    'Joanne', 'Joanne Huang'; ...
    {'Divya', 'divya', 'divya '}, 'Divya Chhabra'; ...
    'Eliana', 'Eliana Hechter'; ...
    'Or', 'Or Zuk'}; % Give full name (and email?) of curator
% 	Diana J Zhu <dianazhu@MIT.EDU>, Carlos A Castellanos <carlito@MIT.EDU>, Joanne Huang <joanne89@mit.edu>, Divya Chhabra <divyac@mit.edu>

for i=1:size(curators_list,1)
    data.Person_annotating = strrep_cell(data.Person_annotating, curators_list{i,1}, ...
        curators_list{i,2});
end

num_cutoffs = length(data_params.alpha_vec); % number of different power cutoffs
num_snps_total = length(data.SNPs);
num_diseases = length(data_params.trait_name);
disease_data = struct_by_inds(data, data_params.unique_inds);

R = loadcellfile(html_template_file); % load template file
start_ind = strfind_cell(lower(R), '<tr>');
R_start = R(1:start_ind(1)); R_end = R(start_ind(end):end); % Split to too
title_line_inds = strfind_cell(R, 'Template');

% Prepare the master file
%matlab_master_html_libs_file = fullfile(matlab_libs_dir, 'matlab_utils.html');

disease_names = strrep_cell(data_params.trait_name, {' ', '/', '\', ''''}, '_');
my_mkdir(root_disease_dir);

if(~isfield(disease_data, 'h_scale'))% Fill empty fields
    disease_data.h_scale = cell(num_diseases,1);
end
if(~isfield(disease_data, 'h_ref'))% Fill empty fields
    disease_data.h_ref = cell(num_diseases,1);
end
if(~isfield(disease_data, 'lambda_s_ref'))% Fill empty fields
    disease_data.lambda_s_ref = cell(num_diseases,1);
end
if(~isfield(disease_data, 'num_loci'))% Fill empty fields
    disease_data.num_loci = zeros(num_diseases,1);
end
if(~isfield(data_params, 'h_liab_power_adjusted'))
    data_params.h_liab_power_adjusted = repmat(vec2row(data_params.h_liab), num_cutoffs, 1);
end
if(~isfield(data_params, 'lambda_s_power_adjusted'))
    data_params.lambda_s_power_adjusted = repmat(vec2row(data_params.lambda_s), num_cutoffs, 1);
end


%disease_data.trait_type  = strrep_cell(disease_data.trait_type, 'Quantitative', 'QTL'); % Correct some names for easier display
disease_data.trait_type  = strrep_cell(disease_data.trait_type, 'QTL', 'Quantitative'); % Correct some names for easier display

quant_inds = find(disease_data.trait_type_num == 1);
binary_inds = find(disease_data.trait_type_num == 0);
disease_data.Prevalence = num2cell(disease_data.Prevalence);
disease_data.lambda_s_familial = num2str_cell(num2cell(disease_data.lambda_s_familial));
disease_data.discovery_num_controls = num2cell(disease_data.discovery_num_controls);
disease_data.h_familial = num2str_cell(num2cell(100*disease_data.h_familial)); % plot in percentage

for i= vec2row(quant_inds)
    disease_data.lambda_s_familial{i} = '-';
    disease_data.lambda_s_ref{i} = '-';
    disease_data.Prevalence{i} = '-';
    disease_data.discovery_num_controls{i} = '-';
end
for i=vec2row(find(data.base_prevalence_inds(data_params.unique_inds))) % don't show default arbitrarily set prevalence
    disease_data.Prevalence{i} = '-';
end
for i=vec2row(find(data.base_lambda_s_inds(data_params.unique_inds))) % don't show default arbitrarily set lambda_s
    disease_data.lambda_s_familial{i} = '-';
end

disease_data.h_familial = strrep_cell(disease_data.h_familial, {'-100', '-1'}, '-');
disease_data.lambda_s_familial = strrep_cell(disease_data.lambda_s_familial, '-1', '-');
data.disease_alpha_vec = repmat(5*10^(-8), num_snps_total, 1); % for now assume all have the same power.
data.disease_replication_alpha_vec = repmat(0.01, num_snps_total, 1); % for now assume all have the same power.
data.disease_combined_alpha_vec = repmat(5*10^(-8), num_snps_total, 1); % for now assume all have the same power.

for i=1:num_diseases % loop on diseases to set power configuration
    trait_inds = strmatch(data_params.trait_name{i}, data.Trait); % get all indices of current trait
    switch data_params.trait_name{i}
        case 'Breast Cancer'
            data.disease_alpha_vec(trait_inds) =  5*10^(-7); % special threshold for breast cancer
        case 'Crohns'
            %            data.disease_alpha_vec(trait_inds) =  5*10^(-6); % This 'kills'
            data.disease_combined_alpha_vec(trait_inds) =  5*10^(-8);
        case 'Height'
            
        case 'Lipids'
            
            
            
    end
end

good_familial_inds_vec = ~( isempty_cell(str2nums_cell(disease_data.h_familial)) & ...
    isempty_cell(str2nums_cell(disease_data.lambda_s_familial)) );

% Set names to appear in big table
catalog_html_file = fullfile(root_disease_dir, 'common_diseases_table.html');
[S_main_table S_main_table_cell disease_table_file ...
    h_explained_total lambda_s_explained_total ...
    num_shadow_snps_vec power_correct] = ... %  disease_power_inds2] = ...
    write_disease_table_main_html(data, disease_data, data_params, ...
    disease_names, num_diseases, catalog_html_file, root_disease_dir, ...
    good_familial_inds_vec, R_start, R_end, plot_main_figures_flag); % save the main website

if(~isfield(data, 'Power'))% Fill empty fields
    data.Power = cell(num_snps_total,num_cutoffs);
end
if(~isfield(data, 'lambda_s_snp'))% Fill empty fields
    data.lambda_s_snp = cell(num_snps_total,1);
end
if(~isfield(data, 'h_snp'))% Fill empty fields
    data.h_snp = cell(num_snps_total,1);
end
% for i=1:num_snps_total % copy power from literature. Exclude this for now!!!!
%     if(~isempty(data.Power{i}))
%         i_is = i
%         data_params.Power(i) = str2nums(data.Power{i}, 1, 1); % get first number. Convert percentage to nums
%     end
% end
data.Power = data_params.Power;
data.chr_str = chr_num2str(data.chr, 'hg18');
% Special: breast cancer

all_inds_vec = []; fit_string = {};
lambda_s_explained_total_percent = cell(num_diseases,1);
all_var_explained_vec = cell(1,num_diseases); % use a row vector
S_concatenated_QTL = []; S_concatenated_binary = [];
S_concatenated_QTL_tab = []; S_concatenated_binary_tab = [];
min_QTL_ind = find(disease_data.trait_type_num == 1, 1);
min_binary_ind = find(disease_data.trait_type_num == 0, 1);

for i=1:num_diseases % loop on diseases
    save_html_disease_i_is = i
    if(i==44) % Height
       TTTTTT = 132412341234 
    end
    my_mkdir(fullfile(root_disease_dir, disease_names{i})); % make disease-specific directory
    S_tab = strrep_cell(R_start, 'Common Diseases', disease_names{i});
    S_tab = strrep_cell(S_tab, 'common diseases', disease_names{i});
    S_tab = strrep_cell(S_tab, '<script src="sorttable.js">', '<script src="../sorttable.js">'); % link to script in root dir
    trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait
    num_snps = length(trait_inds);
    
    S = cell(num_snps+2,18+num_corrections);     % Set names to appear in disease-specific table
    switch disease_data.trait_type{i}
        case {'Binary', 'Disease'}
            effect_size_str = 'OR';
        case {'QTL', 'Quantitative'}
            effect_size_str = '&beta;'; % Need html to display beta
    end
    
    S(1,:) = [{'#', 'SNP-ID', 'chr', 'pos', 'Gene', ... % 1-5
        'RAF', effect_size_str, 'p-val', 'Power', 'Power (emp.)',  ... % 6-10
        'Non-Cent.-Par.','&lambda;<sub>s,i</sub>', 'Var-expl. (&lambda;<sub>s</sub>)', 'Var-expl. (&lambda;<sub>s</sub>, pow. adjust)', 'Var-expl.', ... % 11-15
        'h<sup>2</sup> expl.', 'Var-expl. (pow.-adjust.)'}, ...
        correction_str  '# shadow-loci'];
    %        'h<sup>2</sup> (pow.-adjust.)', '# shadow-loci'}; % 16-19
    %    S{1,10} = '<br>';
    
    cur_alpha_used = data.disease_alpha_vec(trait_inds(1));
    [~, cur_ind_used] = min(abs(data_params.alpha_vec - cur_alpha_used)); % choose with what cutoff to work for this disease
    disease_power_inds(i) = cur_ind_used;
    
    lambda_s_var_explained_vec = zeros(num_snps,1);
    
    correction_mode = 'floor'; correction_type = 'theoretical'; % 'empirical';
    %     switch correction_type
    %         case 'empirical'
    %             power_correct = data_params.Power_empirical(trait_inds,:);
    %         case 'theoretical'
    %             power_correct = data.Power(trait_inds,:);
    %     end
    %     switch correction_mode
    %         case 'floor' % use the same based on #shadow loci
    %             num_shadow_snps_vec = floor(1 ./ max(MIN_POWER, power_correct))-1; % count how many shadow SNPs are there for each SNP
    %         case {'round'} % use the same based on #shadow loci
    %             num_shadow_snps_vec = round(1 ./ max(MIN_POWER, power_correct))-1; % count how many shadow SNPs are there for each SNP
    %         case 'exact'
    %             num_shadow_snps_vec = (1 ./ max(MIN_POWER, power_correct))-1; % count how many shadow SNPs are there for each SNP
    %     end
    h_explained_power_adjusted_vec = ...
        repmat(data_params.snp_h_liab(trait_inds), 1, num_cutoffs) .* ...
        (num_shadow_snps_vec{i}+1);  %   max(MIN_POWER, data.Power(trait_inds,:));
    h_explained_power_adjusted_liberal_vec = repmat(data_params.snp_h_liab(trait_inds), 1, num_cutoffs) ./ ...
        max(MIN_POWER_LIBERAL, power_correct{i,1}); % a more aggresive power correction (make big corrections for loci with low power)
    num_shadow_snps_exact_vec = (1 ./ max(MIN_POWER, power_correct{i,1}))-1; % count how many shadow SNPs are there for each SNP
    
    for j=1:num_snps % loop on snps for a given disease
        S{j+1,1} = num2str(j);
        S{j+1,2} = data.SNPs{trait_inds(j)}; % SNP ID
        S{j+1,3} = data.chr_str{trait_inds(j)}; %   chr_num2str(data.chr(trait_inds(j)), 'hg18');
        S{j+1,4} = num2str(data.pos(trait_inds(j))); % genomic position
        S{j+1,5} = data.Gene{trait_inds(j)}; % closest reported gene
        S{j+1,6} = data.RAF(trait_inds(j)); % Risk Allele Frequency
        switch disease_data.trait_type{i} % Effect size. Take REPLICATION
            case {'Binary', 'Disease'}
                S{j+1,7} = data.OR(trait_inds(j),2);
            case {'QTL', 'Quantitative'}
                S{j+1,7} = data.Beta(trait_inds(j),2);
        end
        S{j+1,8} = data.P_value{trait_inds(j)}; % p-value in initial study
        S{j+1,9} = num2str(data.Power(trait_inds(j),cur_ind_used), 4); % analytical power
        S{j+1,10} = num2str(data_params.Power_empirical(trait_inds(j),cur_ind_used), 4); % empirical power
        S{j+1,11} = num2str(data_params.non_centrality(trait_inds(j)), 4); % rounded NCP
        
        S{j+1,12} = num2str(data_params.snp_lambda_s(trait_inds(j))); % sibling relative risk attributed to each snp
        cur_lambda_s_familial = str2double(disease_data.lambda_s_familial{i});
        if(~isempty(cur_lambda_s_familial)) % var. explained on sib-risk scale attributed to each snp
            S{j+1,13} = num2str(100*log(data_params.snp_lambda_s(trait_inds(j))) / log(cur_lambda_s_familial), 3);
            lambda_s_var_explained_vec(j) = 100*log(data_params.snp_lambda_s(trait_inds(j))) / log(cur_lambda_s_familial);
        else
            S{j+1,13} = '0';
            lambda_s_var_explained_vec(j) = 0;
        end
        S{j+1,15} = num2str(100*data_params.snp_h_liab(trait_inds(j)), 3); % variance explained for each SNP       S{j+1,10} = ' <br>'; % Move to next line
        S{j+1,16} = num2str(100*data_params.snp_h_liab(trait_inds(j))/ ...
            (str2double(disease_data.h_familial{i})/100), 3); % heritability explained for each SNP       S{j+1,10} = ' <br>'; % Move to next line
        
        for jj=1:num_corrections % add values for different corrections
            if( (length(data_params.VarExplainedCorrected) >= i) && ~isempty(data_params.VarExplainedCorrected{i})) % check we actually computed sophisticated correction
                S{j+1,17+jj} = num2str( ...
                    100*data_params.VarExplainedCorrected{i}(j,correction_inds(jj)) / ... %   h_explained_power_adjusted_vec(j,cur_ind_used) /...
                    (str2double(disease_data.h_familial{i})/100) , 3); % variance explained for each SNP (power adjusted)
            else % just copy OLD correction again and again
                S{j+1,17+jj} = num2str( ...
                    100*h_explained_power_adjusted_vec(j,cur_ind_used) /...
                    (str2double(disease_data.h_familial{i})/100) , 3); % variance explained for each SNP (power adjusted)
            end
        end
        switch correction_mode % here do all power corrections
            case 'exact'
                S{j+1,14} = num2str(lambda_s_var_explained_vec(j) / ...
                    max(MIN_POWER, power_correct{i,1}(j)), 3); % lambda_s power adjusted
                S{j+1,17} = num2str(100*h_explained_power_adjusted_vec(j,cur_ind_used), 3); % variance explained for each SNP (power adjusted)
                %                S{j+1,18} = num2str(100*h_explained_power_adjusted_vec(j,cur_ind_used) / ...
                %                    (str2num(disease_data.h_familial{i})/100) , 3); % variance explained for each SNP (power adjusted)
            case {'floor', 'round'}
                S{j+1,14} = num2str(lambda_s_var_explained_vec(j) * ...
                    (num_shadow_snps_vec{i}(j,cur_ind_used)+1), 3); % lambda_s power adjusted
                S{j+1,17} = num2str(100*h_explained_power_adjusted_vec(j,cur_ind_used), 3); % variance explained for each SNP (power adjusted)
                %                S{j+1,18} = num2str(100*h_explained_power_adjusted_vec(j,cur_ind_used) /...
                %                    (str2num(disease_data.h_familial{i})/100) , 3); % variance explained for each SNP (power adjusted)
                S{j+1,18+num_corrections} = num2str(num_shadow_snps_vec{i}(j,cur_ind_used)); % print # shadow loci
        end
    end
    h_explained_power_adjusted_different_cutoff_sensitivity_vec = 100*sum(h_explained_power_adjusted_vec);
    %    h_explained_total(i) = sum(data_params.snp_h_liab(trait_inds));
    % % %     if(disease_specific_plots) % Plot figures for each disease
    % % %         % Make some plots to debug power calculations
    % % %         figure; plot(data_params.non_centrality(trait_inds), data.Power(trait_inds,:), '.');
    % % %         xlabel('non-centrality'); ylabel('power');
    % % %         figure; plot(data_params.snp_h_liab(trait_inds), data_params.non_centrality(trait_inds), '.');
    % % %         xlabel('Var-explained'); ylabel('non-centrality');
    % % %         figure; plot(data_params.snp_h_liab(trait_inds), data.Power(trait_inds,:), '.');
    % % %         xlabel('Var-explained'); ylabel('Power');
    % % %         close all;
    % % %     end
    S{end,1} = '<b>Total:</b>';    % Create a sum table
    for k = 13:(18+num_corrections) % 11:11 % For now only variance explained ...
        %        h_explained_total = sum(cell2mat(S(2:end-1,k))) / 100;
        switch k
            case 13 % total sibling relative risk explained
                if(~isempty(cur_lambda_s_familial))
                    S{end,k} = ['<b>' num2str(sum(log(data_params.snp_lambda_s(trait_inds))) / ...
                        log(cur_lambda_s_familial)*100, 4) '% </b>'];
                else
                    S{end,k} = '<b> 0% </b>';
                end
            case 14 % total sibling relative risk explained (power adjusted)
                if(~isempty(cur_lambda_s_familial))
                    S{end,k} = ['<b>' num2str(sum(lambda_s_var_explained_vec .* ...
                        (num_shadow_snps_vec{i}(:,cur_ind_used)+1)), 4) '% </b>'];
                    %%                    ./   max(MIN_POWER, data.Power(trait_inds,cur_ind_used))), 4) '% </b>'];
                else
                    S{end,k} = '<b> 0% </b>';
                end
            case 15 % total variance explained
                S{end,k} = ['<b>' num2str(h_explained_total(i)*100, 4) '% </b>'];
                data_params.h_liab(i) = h_explained_total(i); % sum(cell2mat(S(2:end-1,k))) / 100; % re-calculate heritability
            case 16 % total h^2 explained
                S{end,k} = ['<b>' num2str(h_explained_total(i)*100 / (str2double(disease_data.h_familial{i})/100), 4) '% </b>'];
            case 17 % total variance explained (power adjusted)
                S{end,k} = ['<b>' num2str(100*sum(h_explained_power_adjusted_vec(:,cur_ind_used)), 4) '% </b>'];
            case 18+num_corrections % num shadow loci
                S{end,k} = ['<b>'  num2str(sum(num_shadow_snps_vec{i}(:,cur_ind_used))) '</b>'];
            otherwise % total h^2 explained  (power adjusted)
                if( (length(data_params.VarExplainedCorrected) >= i) && ~isempty(data_params.VarExplainedCorrected{i})) % check we actually computed sophisticated correction
                    S{end,k} = ['<b>' num2str(100*sum(data_params.VarExplainedCorrected{i}(:,correction_inds(k-17))) / ... %
                        (str2double(disease_data.h_familial{i})/100), 4) '% </b>'];
                else
                    S{end,k} = ['<b>' num2str(100*sum(h_explained_power_adjusted_vec(:,cur_ind_used)) / ...
                        (str2double(disease_data.h_familial{i})/100), 4) '% </b>'];
                end
        end
    end
    for k = 12:12 % 10:10 % For now only sibling relative risk ... take the product
        lambda_s_explained_total(i) = prod(cell2mat(str2num_cell(S(2:end-1,k))));
        S{end,k} = ['<b>' num2str(lambda_s_explained_total(i)) '</b>'];
        tmp_num = str2double(disease_data.lambda_s_familial{i});
        if(~isempty(tmp_num))
            lambda_s_explained_total_percent{i} = num2str(100*log(lambda_s_explained_total(i)) / ...
                log(tmp_num));
        else
            lambda_s_explained_total_percent{i} = '-';
        end
    end
    for j=1:num_snps % add '%' to all snps
        for k = 13:(17+num_corrections)
            S{j+1,k} = [num2str(S{j+1,k}) '%'];
        end
    end
    
    %    h_explained_power_adjusted = h_explained_total; % New divide by power
    
    
    switch correction_mode
        case 'exact'
            data_params.h_liab_power_adjusted(:,i) = ... % h_explained_power_adjusted = ...
                100 * sum(repmat(data_params.snp_h_liab(trait_inds), 1, num_cutoffs) ./ ...
                max(MIN_POWER, power_correct{i,1})); % divide by power. Take minimum: 0.05
            data_params.lambda_s_power_adjusted(:,i) = ... %    lambda_s_explained_power_adjusted =
                prod(repmat(data_params.snp_lambda_s(trait_inds), 1, num_cutoffs) .^ ...
                (1 ./ max(MIN_POWER, power_correct{i,1}))); % divide by power. Take minimum: 0.05
        case {'round', 'floor'} % use the same based on #shadow loci
            data_params.h_liab_power_adjusted(:,i) = ...
                100 * sum(repmat(data_params.snp_h_liab(trait_inds), 1, num_cutoffs) .* ...
                (num_shadow_snps_vec{i}+1)); % divide by power. Take minimum: 0.05
            data_params.lambda_s_power_adjusted(:,i) = ... %    lambda_s_explained_power_adjusted =
                prod(repmat(data_params.snp_lambda_s(trait_inds), 1, num_cutoffs) .^ ...
                (num_shadow_snps_vec{i}+1)); % divide by power. Take minimum: 0.05s            
            %        case 'floor'            
    end
    
    tmp_num = str2double(disease_data.lambda_s_familial{i});
    if(~isempty(tmp_num))
        data_params.lambda_s_power_adjusted_percent{i} = num2str(100*log(data_params.lambda_s_power_adjusted(i)) / ...
            log(tmp_num));
    else
        data_params.lambda_s_power_adjusted_percent{i} = '-';
    end
    switch disease_data.trait_type{i}
        case {'Binary', 'Disease'}
        case {'QTL', 'Quantitative'} % remove 10, 12,13,14
            S = S(:,[1:9 11 15:end]); % don't need lambda_s_i and empirical power (doesn't exist)
    end
    
    S_tab{end+1} = ['This is a <i>' disease_data.trait_type{i} '</i> trait. ' ... % Add trait's type
        'Associations reported in the following publication:<br>']; % Add some data
    S_tab{end+1} = ['<b> ' disease_data.Study{i} '</b>, <i>' disease_data.First_Author{i} ' et al. , ' ...
        disease_data.Journal{i} ', ' str2word('/', disease_data.Date{i}, 3) ...
        ' </i>' html_write_link(['pubmed: ' disease_data.PUBMEDID{i}], ...
        ['http://www.ncbi.nlm.nih.gov/pubmed/' disease_data.PUBMEDID{i}]) '<br> <br>'];
    
    S_tab{end+1} = 'Summary Trait''s Parameters: <br>'; % This should appear Before table !
    switch disease_data.trait_type{i} % remove irrelevant fields according to trait's type
        case {'Binary', 'Disease'}
            retain_inds = 1:size(S_main_table_cell,2);
        case {'QTL', 'Quantitative'}
            remove_fields = {'Prevalence',  '# Controls', 'Sib-Risk(&lambda;<sub>s</sub>)', ...
                'Effective-Sample-Size', '&lambda;<sub>s</sub> explained'};
            [~, ~, remove_inds] = intersect(remove_fields, S_main_table_cell(1,:));
            retain_inds = setdiff(1:size(S_main_table_cell,2), remove_inds);
    end
    cur_S_main_table = vec2row(html_write_table(S_main_table_cell(:,retain_inds), 1)); % This is the main table (to be used later)
    
    S_tab = [S_tab' cur_S_main_table([1 2 i+2 end])]'; % add line with main statistics for the specific trait
    S_tab{end+1} = ['<br>SNP-Specific Parameters: (' ...
        html_write_link('download', remove_dir_from_file_name(file_name_to_txt(disease_table_file{i}))) ...
        ' as tab-delimited .txt.) <br>']; % This should appear After table !
    
    %     S_tab{end+1} = ['Total genetic effect and fraction explained: <br>']; % Add two figures showing contribution of genetics
    %     S_tab{end+1} = '<a href="genetic_effect_explained_summary.jpg" target="_blank">'; % add link to image
    %     S_tab{end+1} = ['<img src="genetic_effect_explained_summary.jpg" ' ...
    %         'alt="Estimated and explained genetic effect size for &lambda;<sub>s</sub> ' ...
    %         'and h<sup>2</sup>"width="500" height="400"></img> </a><br>']; % Add figure showing effect size explained
    %     S_tab{end+1} = ['Cumulative effect size of all discovered loci: <br>'];
    %     S_tab{end+1} = '<a href="all_loci_cumulative_effect_size.jpg" target="_blank">'; % add link to image
    %     S_tab{end+1} = ['<img src="all_loci_cumulative_effect_size.jpg" ' ...
    %         'alt="Cumulative effect size"width="500" height="400"></img> </a><br>']; % Add figure showing effect size explained
    
    %     S_tab{end+1} = ['Effect size (var. explained) of all discovered loci: <br>'];
    %     S_tab{end+1} = '<a href="all_loci_variance_explained.jpg" target="_blank">'; % add link to image
    %     S_tab{end+1} = ['<img src="all_loci_variance_explained.jpg" ' ...
    %         'alt="Effect size (var. explained)"width="500" height="400"></img> </a><br>']; % Add figure showing effect size explained
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    if(disease_specific_plots) % Plot figures for each disease
    [var_explained_vec var_explained_extrapolated_vec ...
        num_loci_needed num_snps_inflated fit_string{i} ...
        true_beta_mu{i} true_beta_std{i}] = ...
        write_disease_specific_plots(disease_specific_plots, ...
        disease_data, data, data_params, i, ...
        lambda_s_explained_total_percent{i}, good_familial_inds_vec, ...
        num_shadow_snps_vec(i,:), h_explained_power_adjusted_vec, lambda_s_var_explained_vec, ...
        root_disease_dir, disease_names, power_correct(i,:), cur_ind_used);
    close all;
    S = [S ['Gaussian &mu;' num2cell(true_beta_mu{i})' '0']' ...
        ['Gaussian &sigma;' num2cell(true_beta_std{i})' '0']']; % add mean and st.d. at the end
    savecellfile(strrep_cell(strrep_cell(S, '<b>', ''), '</b>', ''), ... % save disease-specific data as .txt (tab-delimited, can be viewed in excel)
        fullfile(root_disease_dir, file_name_to_txt(disease_table_file{i})));
    
    cur_trait_vec = [];
    switch lower(disease_data.trait_type{i})
        case {'qtl', 'quantitative'}
            if(i == min_QTL_ind) % first Quantitative 
                S_concatenated_QTL_tab = S_tab(1:end-17);
                S_concatenated_QTL_tab{end-1} = '<H2> All QTL SNPs Parameters </H2>';
                S_concatenated_QTL_tab{end+1} = ['<br>Parameters for all QTL SNPs: (' ...
                    html_write_link('download', ...
                    remove_dir_from_file_name(file_name_to_txt(strrep(catalog_html_file, ...
                    'table', 'all_loci_concatenated_QTL')))) ...
                    ' as tab-delimited .txt.) <br>']; % This should appear After table !
                cur_trait_vec = {'Trait'};
            end
            cur_trait_vec = [cur_trait_vec repmat(disease_data.Trait(i), 1, num_snps)]';
            S_concatenated_QTL = [S_concatenated_QTL' ...
                [cur_trait_vec S(min(i-min_QTL_ind+1,2):end-1,:)]']'; % concatenate all traits (excel)
            
            %             S_concatenated_QTL_tab = [S_concatenated_QTL_tab' ...
            %                 vec2row(html_write_table(S(min(i-min_QTL_ind+1,2):end-1,:), 1))]'; % concatenate all traits (html)
        case {'binary', 'disease'}
            if(i == min_binary_ind)
                S_concatenated_binary_tab = S_tab(1:end-17);
                S_concatenated_binary_tab{end-1} = '<H2> All binary SNPs Parameters </H2>';
                S_concatenated_binary_tab{end+1} = ['<br>Parameters for all binary SNPs: (' ...
                    html_write_link('download', ...
                    remove_dir_from_file_name(file_name_to_txt(strrep(catalog_html_file, ...
                    'table', 'all_loci_concatenated_binary')))) ...
                    ' as tab-delimited .txt.) <br>']; % This should appear After table !
                cur_trait_vec = {'Trait'};
            end
            cur_trait_vec = [cur_trait_vec repmat({disease_data.Trait{i}}, 1, num_snps)]';
            S_concatenated_binary = [S_concatenated_binary' ...
                [cur_trait_vec S(min(i-min_binary_ind+1,2):end-1,:)]']'; % concatenate all traits (excel)
            %             S_concatenated_binary_tab = [S_concatenated_binary_tab' ...
            %                 vec2row(html_write_table(S(min(i-min_binary_ind+1,2):end-1,:), 1))]'; % concatenate all traits (html)
    end
    S = num2str_cell(empty_cell_to_empty_str(S));
    S_tab = [S_tab' vec2row(html_write_table(S, 1))]'; % add table
    S_tab{end+1} = '<br>'; % This should appear AFTER table !
    S_tab{end+1} = ['Total variance explained by all found loci: ' num2str(data_params.h_liab(i)*100) '% <br>'];
    S_tab{end+1} = ['Total variance explained adjusting for power: ' num2str(data_params.h_liab_power_adjusted(i)) '% <br>'];
    S_tab{end+1} = ['Estimated variance explained by genetics (heritability): ' disease_data.h_familial{i} '% <br>'];
    switch disease_data.trait_type{i} % show sibling relative risk for binary traits
        case 'Binary'
            S_tab{end+1} = ['Total &lambda;<sub>s</sub> explained by all found loci: ' ...
                num2str(lambda_s_explained_total(i)) ' (' ...
                lambda_s_explained_total_percent{i} '% of &lambda;<sub>s</sub>) <br>'];
            S_tab{end+1} = ['Total &lambda;<sub>s</sub> explained adjusting for power: ' ...
                num2str(data_params.lambda_s_power_adjusted(i)) ' (' ...
                data_params.lambda_s_power_adjusted_percent{i} '% of &lambda;<sub>s</sub>) <br>'];
            S_tab{end+1} = ['Estimated &lambda;<sub>s</sub> from sibbling studies: ' disease_data.lambda_s_familial{i} ' <br><br>'];
    end
    
    
    %    end % if plot
    all_var_explained_vec{i} = var_explained_vec{1};
    all_inds_vec = [all_inds_vec 1:length(var_explained_vec{1})]; % get how many snps are in each disease
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tab = strrep_cell(S_tab, 'fit_string', fit_string{i}); % Get the true string. Then insert disease-specific plots
    S_tab{end+1} = html_write_img('Effect Sizes and MAF of Associated Loci:', 'effect_size_vs_MAF.jpg');
    S_tab{end+1} = html_write_img('Variance explained for different effect sizes (density):', 'var_explained_by_effect_size_smoothed_density.jpg');
    S_tab{end+1} = html_write_img('Variance explained for different effect sizes (cumulative):', 'var_explained_by_effect_size_smoothed_cumulative.jpg');
    
    S_tab{end+1} = html_write_img('Total genetic effect and fraction explained:', 'genetic_effect_explained_summary.jpg');
    S_tab{end+1} = html_write_img('Cumulative effect size of discovered loci:', 'all_loci_cumulative_effect_size.jpg');
    S_tab{end+1} = html_write_img('Effect size (var. explained) of discovered loci:', 'all_loci_variance_explained.jpg');
    if(~isempty(var_explained_extrapolated_vec))
        S_tab{end+1} = 'Cumulative h<sup>2</sup> Explained:';
        n_loci_example_vec = [10 30 50 100 200 300 500 1000 2000 3000 5000];
        S_power_tab = [{'Num. loci', 'h<sup>2</sup> explained'}' ...
            num2str_cell(num2cell([n_loci_example_vec' var_explained_extrapolated_vec(n_loci_example_vec)']), 4)']';
        for j=2:length(S_power_tab)
            S_power_tab{j,2} = [S_power_tab{j,2} '%'];
        end
        S_tab = [S_tab' vec2row(html_write_table(S_power_tab, 4))]';   % add table
    end
    S_tab{end+1} = html_write_img(['Cumulative Var. Explained: ' num2str(num_snps) ' loci. ' ...
        num2str(num_snps_inflated(1)) ' loci (power-corrected). Need N=' ...
        num2str(round(num_loci_needed(1))) ' loci to reach 50% h<sup>2</sup>:'], ...
        'all_loci_cumulative_variance_explained.bmp'); % use just first power correction for display
    switch disease_data.trait_type{i}
        case 'Binary'
            S_tab{end+1} = html_write_img('Var. Explained on two scales', ...
                'var_explained_two_scales.jpg');
    end
    S_tab{end+1} = html_write_img('Sensitiviry of Var. Explained to p-value cutoff for power calculations:', ...
        'sensitivity_variance_explained_pval_cutoff.jpg');
    
    S_tab{end+1} = html_write_img('Var. Explained (log. scale):', 'all_loci_variance_explained_log_scale.bmp');
    S_tab{end+1} = html_write_img('Var. Explained (log. scale):', 'all_loci_variance_explained_log_scale_loglog.bmp');
    %  S_tab{end+1} = html_write_img(['Cumulative Var. Explained (extrapolated). Need N=' ...
    %      num2str(round(num_loci_needed(1))) ' loci to reach 50% h<sup>2</sup>:'], ...
    %      'all_loci_cumulative_variance_explained_extrapolation.jpg'); % report #loci needed based on power-law
    S_tab{end+1} = html_write_img('Cumulative Var. Explained vs. Sample Size', 'var_explained_vs_sample_size.bmp');
    
    
    
    switch disease_data.trait_type{i}
        case 'Binary'
            S_tab{end+1} = ['Effect size (GRR) of all discovered loci: <br>'];
            S_tab{end+1} = '<a href="all_loci_genetic_relative_risk.jpg" target="_blank">'; % add link to image
            S_tab{end+1} = ['<img src="all_loci_genetic_relative_risk.jpg" ' ...
                'alt="Effect size (GRR)"width="500" height="400"></img> </a><br>']; % Add figure showing effect size explained
        case {'QTL', 'Quantitative'}
            S_tab{end+1} = ['Effect size (\beta) of all discovered loci: <br>'];
            S_tab{end+1} = '<a href="all_loci_beta_effect_size.jpg" target="_blank">'; % add link to image
            S_tab{end+1} = ['<img src="all_loci_beta_effect_size.jpg" ' ...
                'alt="Effect size (\beta)"width="500" height="400"></img> </a><br>']; % Add figure showing effect size explained
    end
    
    %    switch disease_data.trait_type{i} % plot power computed via two methods
    %        case 'Binary'
    S_tab{end+1} = html_write_img('Power Empirical vs. Theoretical', 'power_empirical_vs_theoretical.jpg');
    %    end
    
    S_tab{end+1} = ['<i>Detailed Description: </i><b>' disease_data.Detailed_Description{i} '</b><br>'];    % Insert more comments:
    S_tab{end+1} = ['<i>Comments: </i><b>' disease_data.Comments{i} '</b><br>'];    % Insert more comments:
    S_tab{end+1} = ['<i>Curator: </i><b>' data.Person_annotating{i} '</b><br>'];
    
    
    ctr = size(S_tab,1); % add stuff at the end to file
    for j=1:length(R_end)
        S_tab{ctr+j,1} = R_end{j};
    end
    S_tab{end+1} = html_write_link('Back to catalog', ...
        fullfile('../', remove_dir_from_file_name(catalog_html_file)));
    %    dir_is = fullfile(root_disease_dir, disease_names{i})
    
    savecellfile(S_tab, fullfile(root_disease_dir, disease_names{i}, ...
        [disease_names{i} '_table.html']), [], 1); % save disease-specific table
    
end % loop on diseases




S = loadcellfile(catalog_html_file); % Patch: add some fields which we've just computed
first_ind = strfind_cell(S, disease_names{1}); first_ind = first_ind(1);
for i=1:num_diseases
    S{first_ind+i-1} = strrep(S{first_ind+i-1}, 'fit_string', fit_string{i});
end
savecellfile(S, catalog_html_file); % save again main file

% all_inds_vec = length_cell(all_var_explained_vec); % Get how many snps are in each locus
all_var_explained_vec = cell2vec(all_var_explained_vec); % Unite all
good_inds = find((all_var_explained_vec < 1) & (all_var_explained_vec > 0.00000001));

if(plot_main_figures_flag)
    figure; plot(all_inds_vec(good_inds), all_var_explained_vec(good_inds), '.');  % Plot all SNPs effect size on a linear scale
    figure; plot(log(all_inds_vec(good_inds)), log(all_var_explained_vec(good_inds)), '.');  % Plot all SNPs effect size on a log scale
    title('All SNPs effect size vs. rank');
    xlabel('Rank (in a given disease) (log)'); ylabel('Var Explained (log)');
    my_saveas(gcf, fullfile(root_disease_dir, 'figs', 'variance_explained_vs_rank_summary_loglog'), bmp_fig_format); % summary of heritability explained for all diseases
    figure; plot(all_inds_vec(good_inds), log(all_var_explained_vec(good_inds)), '.');  % Plot all SNPs effect size on a semi-log scale
    title('All SNPs effect size vs. rank');
    xlabel('Rank (in a given disease)'); ylabel('Var Explained (log)');
    my_saveas(gcf, fullfile(root_disease_dir, 'figs', 'variance_explained_vs_rank_summary'), bmp_fig_format); % summary of heritability explained for all diseases
    full_figure; % Plot power of all main diseases
    special_traits = strrep_cell(special_traits, {' ', '/', '\', ''''}, '_'); ctr=1;
    plot_inds = zeros(length(special_traits), 1);
    for i=1:min(6,length(special_traits))
        special_ind = strmatch(special_traits{i}, disease_names,  'exact')
        if(~isempty(special_ind))
            plot((1:length(power_correct{special_ind,1}(:,cur_ind_used))) ./ ...
                length(power_correct{special_ind,1}(:,cur_ind_used)), ...
                sort(power_correct{special_ind,1}(:,cur_ind_used)), color_vec(ctr));
            ctr=ctr+1, plot_inds(i) = 1
        end
    end
    title('Power of found loci for different traits');
    xlabel('Rank'); ylabel('Power');
    legend(str2title(special_traits(find(plot_inds))), 4);
    my_saveas(gcf, fullfile(root_disease_dir, 'figs', 'power_of_found_loci_vs_rank_summary'), bmp_fig_format); % summary of heritability explained for all diseases
end


[~, ~, special_inds] = intersect(special_traits, disease_names);
Special_S = strrep_cell(strrep_cell(S_main_table_cell, '<b>', ''), '</b>', '');
remove_fields = {'#', 'Study', 'Sib-Risk(&lambda;<sub>s</sub>)', ...
    '# Cases', '# Controls', ...
    'Effective-Sample-Size', '&lambda;<sub>s</sub> explained', 'fit'};
[~, keep_inds] = setdiff(Special_S(1,:), remove_fields); keep_inds = sort(keep_inds)
Special_S = (Special_S([1 special_inds+1],keep_inds)); % Save detailed summary only for special traits
for i=1:size(Special_S,1)
    for j=1:size(Special_S,2)
        Special_S{i,j} = str2word('">', Special_S{i,j}, 'end');
    end
end
Special_S = strrep_cell(Special_S, '</a>', '');
savecellfile(Special_S,  ... % save as .txt (tab-delimited, can be viewed in excel)
    [remove_suffix_from_file_name(catalog_html_file) '_special.txt'], [], 1);



S_concatenated_QTL = num2str_cell(empty_cell_to_empty_str(S_concatenated_QTL));
S_concatenated_QTL_tab = [S_concatenated_QTL_tab' ...
    vec2row(html_write_table(S_concatenated_QTL, 1))]'; % concatenate all traits (html)
S_concatenated_QTL_tab{end+1} = '<br>'; % This should appear AFTER table !
S_concatenated_binary = num2str_cell(empty_cell_to_empty_str(S_concatenated_binary));
S_concatenated_binary_tab = [S_concatenated_binary_tab' ...
    vec2row(html_write_table(S_concatenated_binary, 1))]'; % concatenate all traits (html)
S_concatenated_binary_tab{end+1} = '<br>'; % This should appear AFTER table !


ctr = size(S_concatenated_QTL_tab,1); % add stuff at the end to file
for j=1:length(R_end)
    S_concatenated_QTL_tab{ctr+j,1} = R_end{j};
end
S_concatenated_QTL_tab{end+1} = html_write_link('Back to catalog', ...
    fullfile(remove_dir_from_file_name(catalog_html_file)));
ctr = size(S_concatenated_binary_tab,1); % add stuff at the end to file
for j=1:length(R_end)
    S_concatenated_binary_tab{ctr+j,1} = R_end{j};
end
S_concatenated_binary_tab{end+1} = html_write_link('Back to catalog', ...
    fullfile(remove_dir_from_file_name(catalog_html_file)));




savecellfile(S_concatenated_binary_tab, strrep(catalog_html_file, ...
    'table', 'all_loci_concatenated_binary'), [], 1); % Save a table of ALL loci concatenated
savecellfile(S_concatenated_QTL_tab, strrep(catalog_html_file, ...
    'table', 'all_loci_concatenated_QTL'), [], 1); % Save a table of ALL loci concatenated
savecellfile(strrep_cell(strrep_cell(S_concatenated_binary, '<b>', ''), '</b>', ''), ... % save disease-specific data as .txt (tab-delimited, can be viewed in excel)
    file_name_to_txt(strrep(catalog_html_file, ...
    'table', 'all_loci_concatenated_binary')), [], 1);
savecellfile(strrep_cell(strrep_cell(S_concatenated_QTL, '<b>', ''), '</b>', ''), ... % save disease-specific data as .txt (tab-delimited, can be viewed in excel)
    file_name_to_txt(strrep(catalog_html_file, ...
    'table', 'all_loci_concatenated_QTL')), [], 1);



