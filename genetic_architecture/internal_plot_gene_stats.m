% Plot site-freq. specturm for a specific gene (This is temporariliy located in in it's own function)
%
% Input:
% gene_header - header representing master directory for dir
% gene_ind - index of gene in GeneStruct
% GeneStruct - Structure with basic genomic information on eacg gene
% SiteFreqSpecStruct - Structure with variants from sequencing information for gene
% MutationRateTable - table of genomic mutation rates for each gene and each mutation class
% MutationTypes - indices of mutation classes
% gene_dir - where to save output for gene
%
% Output: None. Information on gene is saved in file: [gene_dir]/[gene_name]'_Info.txt'
% 
function internal_plot_gene_stats(gene_header, gene_ind, GeneStruct, ExonsGeneStruct, SiteFreqSpecStruct, ...
    MutationRateTable, MutationTypes, gene_dir)

Assign24MammalsGlobalConstants;
alpha_s_fit = 'MLE'; % MLE (how to fit alpha and s) 
%legend_vec = vec2row(SiteFreqSpecStruct{1}.allele_types(SiteFreqSpecStruct{1}.good_allele_inds));
num_populations = length(SiteFreqSpecStruct);

I = gene_ind; gene_name = GeneStruct.gene_names{gene_ind}; % strfind_cell(GeneStruct.gene_names, gene_header); % get gene's index in genes struct
J = strmatch(gene_name, SiteFreqSpecStruct{1}.unique_genes, 'exact'); % get gene's index in SFS struct


if(isempty(I) || isempty(J)) % couldn't find gene in both files
    return;
end
%[sorted_f_vec sort_perm] = sort(SiteFreqSpecStruct{1}.f_vec{J}); % get allele frequencies for this gene

for figure_type = {'num_variants_cum_log', 'heterozygosity_cum_log', 'num_variants_cum_log_normalized', 'heterozygosity_cum_log_normalized'}
    figure; ctr=1; % plot cumulative allele frequency for each mutation type
    [~, allele_types_I, allele_types_J] = intersect(MutationTypes, SiteFreqSpecStruct{1}.allele_types_ind);
    legend_vec = strdiff_cell(SiteFreqSpecStruct{1}.allele_types(SiteFreqSpecStruct{1}.good_allele_inds), 'coding');
    %    for i=1:length(MutationTypes)
    %        cur_mutation_type_ind = allele_types_J(i); %    find(SiteFreqSpecStruct{1}.allele_types_ind == MutationTypes(i));
    for i=1:length(SiteFreqSpecStruct{1}.good_allele_inds) % take only synonymous, missense and stop codons 
        cur_mutation_type_ind = SiteFreqSpecStruct{1}.good_allele_inds(i); % index in list of alleles from ESP
        cur_mutation_type_rate_ind = ... % index in list of genomic mutation types
            find(SiteFreqSpecStruct{1}.allele_types_ind(cur_mutation_type_ind) == MutationTypes); % find index representing mutation
        switch figure_type{1} % Determine x-vec
            case {'num_variants_cum', 'heterozygosity_cum', 'num_variants_cum_log', 'heterozygosity_cum_log', ...
                    'num_variants_cum_log_normalized', 'heterozygosity_cum_log_normalized'}
                [x_vec, freq_sort_perm] = sort( SiteFreqSpecStruct{1}.gene_by_allele_type_freq_list{cur_mutation_type_ind, J} );
                x_label = 'Derived Allele Freq.';
        end
        legend_loc='southeast'; % default legend (for cumulative distributions)
        switch figure_type{1} % determine y_vec
            case {'num_variants_cum', 'num_variants_cum_log'}
                y_vec = (1:length(x_vec)) ./ length(x_vec);
                y_label = 'Cumulative # Variants';
            case {'num_variants_cum_normalized', 'num_variants_cum_log_normalized'}
                y_vec = ( (1:length(x_vec)) ./ length(x_vec) ) ./ MutationRateTable(gene_ind, cur_mutation_type_rate_ind); % normalize by mutation rate per gene
                y_label = 'Cumulative Heterozygosity (Normalized by Target Size)';
                
            case {'heterozygosity_cum', 'heterozygosity_cum_log'}
                y_vec = cumsum(SiteFreqSpecStruct{1}.gene_by_allele_type_het_list{cur_mutation_type_ind, J}(freq_sort_perm));
                y_label = 'Cumulative Heterozygosity';
            case {'heterozygosity_cum_normalized', 'heterozygosity_cum_log_normalized'}
                y_vec = ( cumsum(SiteFreqSpecStruct{1}.gene_by_allele_type_het_list{cur_mutation_type_ind, J}(freq_sort_perm)) ) ./ ...
                    MutationRateTable(gene_ind, cur_mutation_type_rate_ind);
                y_label = 'Cumulative Heterozygosity (Normalized by Target Size)';
        end
        switch figure_type{1} % determine log-scale
            case {'num_variants_cum', 'heterozygosity_cum', 'num_variants_cum_normalized', 'heterozygosity_cum_normalized'}
                log_x_flag = 0; log_y_flag = 0;
            case {'num_variants_cum_log', 'heterozygosity_cum_log'}
                log_x_flag = 1; log_y_flag = 0;
        end
        if(log_x_flag)
            if(log_y_flag)
                loglog(x_vec, y_vec, color_vec(ctr), 'linewidth', 2);
            else
                semilogx(x_vec, y_vec, color_vec(ctr), 'linewidth', 2);
            end
        else
            if(log_y_flag)
                semilogy(x_vec, y_vec, color_vec(ctr), 'linewidth', 2);
            else
                plot(x_vec, y_vec, color_vec(ctr), 'linewidth', 2);
            end
        end
        hold on; ctr=ctr+1;
    end % loop on mutations type
    
    legend(legend_vec, 'location', legend_loc); legend('boxoff'); 
    xlabel(x_label); ylabel(y_label);
    save_flag = 0;
    while(save_flag == 0)
        save_flag = 1;
        try
            my_saveas(gcf, fullfile(gene_dir, [gene_name '_' figure_type{1}]), 'pdf');
        catch exception % some problem with saving to pdf 
            save_flag = 0;
        end
    end
end % loop on figure type


% Save Gene's summary statistics in a text file 
num_allele_types = length(MutationTypes);
R = cell(10,num_allele_types+2); % New: Also save a 'card' with information for each gene

R{1,1} = ['Summary Statistics for Gene: ' gene_name]; % Save mutation target size
R{2,1} = '--------------------------------------------';
ctr=4;
R{ctr,1} = 'Chr:'; R{ctr,2} = GeneStruct.chr_vec(gene_ind);
R{ctr+1,1} = 'Strand:'; R{ctr+1,2} = num2strand(GeneStruct.strand(gene_ind));

% R{ctr+1,1} = 'Start:'; R{ctr+1,2} = num2str(GeneStruct.start_vec(gene_ind));
R{ctr+2,1} = 'Exonic-Length ';
R{ctr+2,2} = [num2str(length(GeneStruct.seqs{gene_ind})) ' (nt) ' ...
    num2str(length(GeneStruct.seqs{gene_ind})/3) ' (AA)']; ctr=ctr+1;
R{ctr+3,1} = 'Exons: '; R{ctr+3,2} = 'start'; R{ctr+3,3} = 'end'; R{ctr+3,4} = 'length (nt)';

gene_exon_inds = strmatch(gene_name, ExonsGeneStruct.gene_names, 'exact');  num_exons = length(gene_exon_inds);
for j=1:num_exons % Plot individual exons. Need the original GeneSturct (non-unique)!!!!
    R{ctr+j+3,1} = num2str(j);
    R{ctr+j+3,2} = num2str(ExonsGeneStruct.pos_start_vec(gene_exon_inds(j)));
    R{ctr+j+3,3} = num2str(ExonsGeneStruct.pos_end_vec(gene_exon_inds(j)));
    R{ctr+j+3,4} = num2str(abs(  ExonsGeneStruct.pos_end_vec(gene_exon_inds(j)) - ...
        ExonsGeneStruct.pos_start_vec(gene_exon_inds(j)) ) + 1);
end

ctr=ctr+6+num_exons; % separate mutation rates

R{ctr,1} = 'Mutation Rates and Target Size for Different Classes (ESP Data):'; ctr=ctr+1;
R{ctr,1} = '----------------------------------------------------------------'; ctr=ctr+1;
R{ctr,1} = 'Class:'; R{ctr+1,1} = 'Target Size:';
R{ctr+2,1} = 'Population:';

R{ctr+3,1} = '# Alleles:'; R{ctr+4,1} = 'Mean DAF:'; % R{ctr+5,1} = 'Cumulative freq.:';
R{ctr+6,1} = 'Heterozygosity:';
R{ctr+7+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),1} = 'Fitted s:';
R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),1} = 'Fitted \alpha (crude!):';

total_num_alleles_by_type = length_cell(SiteFreqSpecStruct{1}.gene_by_allele_type_freq_list(:,J));
[~, mutation_type_inds] = intersect(SiteFreqSpecStruct{1}.allele_types_ind, MutationTypes);
total_num_alleles = sum(total_num_alleles_by_type(mutation_type_inds));
allele_freq_vec = zeros(total_num_alleles,num_populations);
allele_n_vec = zeros(total_num_alleles,num_populations);
allele_position_vec = zeros(total_num_alleles,1);
allele_type_vec = cell(total_num_alleles,1); % same allele type for all populations
allele_reference_vec = cell(total_num_alleles,1);
allele_derived_vec = cell(total_num_alleles,1);
allele_amino_acid_change_vec = cell(total_num_alleles,1);

allele_ctr=1;
for i=1:length(MutationTypes)
    cur_mutation_type_ind = find(SiteFreqSpecStruct{1}.allele_types_ind == MutationTypes(i));  % index in list of alleles from ESP
    current_num_alleles = total_num_alleles_by_type(cur_mutation_type_ind);
    
    R{ctr,num_populations*i+1} = genome_types{MutationTypes(i)}; % string with mutation type
    R{ctr+1,num_populations*i+1} = MutationRateTable(gene_ind, i); % Total mutation rate (target size)
    for k=1:num_populations % here population-specific values
        non_zero_freq_alleles = find( SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{cur_mutation_type_ind, J} > 0); % New! take onle allele freq. > 0
        R{ctr+2,num_populations*i+k} = SiteFreqSpecStruct{k}.population_str; % set population
        R{ctr+3,num_populations*i+k} = length(non_zero_freq_alleles); %%% length(SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{cur_mutation_type_ind, J}(non_zero_freq_alleles)); % write # of alleles
        R{ctr+4,num_populations*i+k} = mean(SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{cur_mutation_type_ind, J}(non_zero_freq_alleles)); % write mean allele freq.
        %    R{ctr+4,i+1} = R{ctr+2,i+1}*R{ctr+3,i+1};
        %    R{ctr+5,i+1} = sum(SiteFreqSpecStruct{k}.gene_by_allele_type_het_list{cur_mutation_type_ind, J}(non_zero_freq_alleles)); % write total heterozygosity
        
        for j=1:length(SiteFreqSpecStruct{k}.upper_freq_vec) % loop on thresholds
            R{ctr+5+j,1} = ['Cumulative freq. % (<' num2str(100*SiteFreqSpecStruct{k}.upper_freq_vec(j), 3) '%):'];
            R{ctr+5+j,num_populations*i+k} =  100*SiteFreqSpecStruct{k}.total_freq_per_gene_mat{j}(cur_mutation_type_ind, J);
            R{ctr+5+length(SiteFreqSpecStruct{k}.upper_freq_vec)+j,1} = ['Heterozygisity %(<' num2str(100*SiteFreqSpecStruct{k}.upper_freq_vec(j), 3) '%):'];
            R{ctr+5+length(SiteFreqSpecStruct{k}.upper_freq_vec)+j,num_populations*i+k} =  100*SiteFreqSpecStruct{k}.total_heterozygosity_per_gene_mat{j}(cur_mutation_type_ind, J);
        end
        R{ctr+7+2*length(SiteFreqSpecStruct{k}.upper_freq_vec),num_populations*i+k} = '-'; % currently we don't have s
        R{ctr+8+2*length(SiteFreqSpecStruct{k}.upper_freq_vec),num_populations*i+k} = '-'; % no alpha yet either
        
        allele_freq_vec(allele_ctr:allele_ctr+current_num_alleles-1,k) = vec2column(SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{cur_mutation_type_ind, J}(:));
        allele_n_vec(allele_ctr:allele_ctr+current_num_alleles-1,k) = vec2column(SiteFreqSpecStruct{k}.gene_by_allele_type_n_list{cur_mutation_type_ind, J}(:));
        if(k==1)
            allele_type_vec(allele_ctr:allele_ctr+current_num_alleles-1) = ...
                repmat({genome_types{MutationTypes(i)}}, current_num_alleles, 1); % total_num_alleles_by_type(cur_mutation_type_ind))];
            allele_position_vec(allele_ctr:allele_ctr+current_num_alleles-1) = ...
                vec2column(SiteFreqSpecStruct{1}.gene_by_allele_type_pos_list{cur_mutation_type_ind, J}(:));
            allele_reference_vec(allele_ctr:allele_ctr+current_num_alleles-1) = ...
                SiteFreqSpecStruct{1}.REF(SiteFreqSpecStruct{1}.gene_by_allele_type_inds_list{cur_mutation_type_ind, J}); 
            allele_derived_vec(allele_ctr:allele_ctr+current_num_alleles-1) = ...
                SiteFreqSpecStruct{1}.ALT(SiteFreqSpecStruct{1}.gene_by_allele_type_inds_list{cur_mutation_type_ind, J}); 
            allele_amino_acid_change_vec(allele_ctr:allele_ctr+current_num_alleles-1) = ...
                SiteFreqSpecStruct{1}.aminoAcidChange(SiteFreqSpecStruct{1}.gene_by_allele_type_inds_list{cur_mutation_type_ind, J}); 
        end
    end % loop on populations
    allele_ctr=allele_ctr+current_num_alleles;    
end % loop on mutation types
allele_num_carriers_vec = round(allele_freq_vec .* allele_n_vec); 


for i=1:length(allele_amino_acid_change_vec) % Parse amino-acid change
    allele_amino_acid_change_vec{i} = cell2vec(unique(strsplit(allele_amino_acid_change_vec{i}, ',')),',');
end


switch alpha_s_fit % Temp: Fit alpha crudely for each gene 
    case 'crude'        
        x_0 = R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),2};
        x_stop = R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),4};
        y = R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),3}; % Set alpha (temp.)
        fitted_alpha = (y-x_0) / (x_stop-x_0); 
        fitted_s = 0; % no fit for s for now .. (why? can use crude fit from stop codons)  
    case 'MLE' % here perform maximum-likelihood reconstruction of alpha and s. Also compute confidence intervals 
        
        beta_vec = []; prevalence = []; trait_type = []; y=[]; full_flag = 0;% no phenotypes 
        rare_cumulative_per_gene = []; % currently this is not used (we use the target size variable)
        N=10000; D=[]; % human effective population size 
        maximize_parameters = [1 1 0]; % maximize over alpha and s (no beta) 
        s_null_vec = 0 -logspace(-6, -0.5, 23); %  take scalar! (initial guess)
        alpha_vec = [0:100] ./ 100; % possible mixture coefficients % take scalar! (initial guess)
        num_individuals = allele_n_vec; % different # of individuals for different alleles !
        X = [vec2row(allele_num_carriers_vec) vec2row(num_individuals)]'; % get allele frequency vector 
        null_w_vec = zeros(total_num_alleles,1); % get class vector
        null_w_vec(strmatch('stop', allele_type_vec)) = 1; 
        null_w_vec(strmatch('synonymous', allele_type_vec)) = 0; 
        null_w_vec(strmatch('missense', allele_type_vec)) = -1; 
        target_size_by_class_vec = MutationRateTable(gene_ind, 1:3); % get target size for each class 
        
        
        % get rid of alleles appearing zero times (they're here since they appear in another population)
        t_max_like = cputime; 
        [fitted_LL, fitted_s, fitted_alpha, fitted_beta] = ...
            maximize_two_class_likelihood(s_null_vec, alpha_vec, beta_vec, rare_cumulative_per_gene, target_size_by_class_vec, N, ...
            X, y, trait_type, prevalence, null_w_vec, maximize_parameters, full_flag, D, num_individuals, 'brute-force'); % get MLE estimator !!! 
        t_max_like = cputime - t_max_like 
end
R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),2+num_populations} = fitted_alpha;
R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),2+num_populations} = 0; % neutral alleles
R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),3+num_populations} = fitted_alpha * fitted_s; % missense (mixture)
R{ctr+8+2*length(SiteFreqSpecStruct{1}.upper_freq_vec),4+num_populations} = fitted_s; % stop-codons


allele_carriers_vec = round(allele_n_vec .* allele_freq_vec);
for i=1:length(R) % get rid of stop-lost
    for k=1:num_populations
        R{i,4*num_populations+k} = '';
    end
end

ctr=ctr+12+2*length(SiteFreqSpecStruct{1}.upper_freq_vec); % need to change
R{ctr,1} = ['List of Individual Alleles: (' SiteFreqSpecStruct{1}.population_str ')']; ctr=ctr+1; % Get individual alleles (need to combine populations)
R{ctr,1} = '-----------------------------------------------------'; ctr=ctr+1;
R{ctr,1} = 'Pos.'; R{ctr,2} = 'Allele-Type';
R{ctr,3} = 'Ref.-Allele'; R{ctr,4} = 'Der.-Allele'; R{ctr,5} = 'AA-Change';
R{ctr,6} = '#Profiled';
R{ctr,6+num_populations} = '#Carriers'; R{ctr,6+2*num_populations} = 'Freq.(%)'; ctr=ctr+1;
R{ctr,1} = 'Population:'; 
for i=1:3
    for k=1:num_populations % here population-specific values
        R{ctr,3+num_populations*i+k} = SiteFreqSpecStruct{k}.population_str; % set population
    end % loop on different populations
end % loop on different fields


%R_alleles = cell(total_num_alleles), num_allele_types+2);
R_alleles = [ allele_position_vec  allele_n_vec allele_carriers_vec ]; % length(allele_n_vec) may be smaller than total_num_alleles
R_alleles = num2str_cell([num2cell(R_alleles(:,1)) allele_type_vec ...
    allele_reference_vec  allele_derived_vec  allele_amino_acid_change_vec ...
    num2cell(R_alleles(:,2:end)) ...
    num2str_cell(num2cell(100*allele_freq_vec), 2) ]);
R = num2str_cell(R', 3);
if(~isempty(R_alleles))
    R_alleles{end,size(R, 1)} = []; % make lengths equal
    R = [R R_alleles']';
else
    R = R'; 
end

R = strrep_cell(R, 'NaN', '-'); % for values we don't have yet
savecellfile(R, fullfile(gene_dir, [gene_name '_Info.txt'])); % save 'report card' for each gene 


