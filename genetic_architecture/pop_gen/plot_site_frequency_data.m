% Plot heterozygosity and site frequency spectrum.
% Note: This function also computes new fields !!!
%
% Input:
% A - structure with allele counts information
% GeneStruct - structure with genomic properties for each gene
% population_str - which populations to plot (optional)
% mutation_rates_files - file with genomic mutation rates information
% n_vec - # of individuals for each allele (?? not used ??)
% count_vec - # of derived allele counts for each allele
% f_vec - frequency of derived alleles for each allele
% allele_types - class of each allele
% target_length - total length in nucleotides of target sequenced
% num_bins - histogram plot resulution
% output_file_name - where to save figures
%
function R = plot_site_frequency_data(A, GeneStruct, exome_struct, mutation_rates_files, n_vec, count_vec, f_vec, ...
    allele_types, target_length, num_bins, output_file_name)

AssignGeneralConstants; Assign24MammalsGlobalConstants; AssignRVASConstants;
bar_plot = 0; % plot bar or histogram
mean_gene_length_in_nt = 2500; % take average exonic length of human genes
if(ischar(A)) % load input data from file
    spectrum_data_file = A; % keep file name
    if(isfield(exome_struct, 'populations') && (~isempty(exome_struct.populations)))
        num_populations = length(exome_struct.populations);
    else
        exome_struct.populations = str2word('_', remove_suffix_from_file_name(spectrum_data_file), 'end');
        num_populations = 1;
    end
    A = cell(num_populations, 1); % 10); % SET DIMENSIONS LATER !
    unite_field_names = {'XXX_FEATURE_', 'GENE', 'XXX_CHROM', 'POS',  'ALLELE_FREQ',   'GENE_INDS', 'unique_genes', 'ProteinPos'}; % list of fields to take in union
    sfs_file_names =  GetFileNames(add_pop_to_file_name(spectrum_data_file, 'AllPop'), 1);
    num_variants_vec = zeros(num_populations, 1); %    spectrum_population_data_file = cell(length(sfs_file_names), 1); 
    
    if(~exist(fullfile(dir_from_file_name(sfs_file_names{1}), [exome_struct.prefix '_AllPop_union.mat']), 'file'))
        for i_c=1:length(sfs_file_names) % , 10) % 10 % TEMP!!! RUN ON FIRST 10 FILES FOR DEBUG.  % loop on all chunks (By chromosomes or otherwise)  % NEW! let many populations !!
            cur_B = load(sfs_file_names{i_c}, 'XXX_REF_ALLELE_COUNT_', 'XXX_VARIANT_COUNT_', 'num_genes', 'unique_genes', ... % 'GENE', ...
                'num_allele_types', 'num_alleles_per_gene_mat', 'total_heterozygosity_per_gene_mat', ...
                'upper_freq_vec', 'total_freq_per_gene_mat', 'num_genes', ...
                'allele_types', 'allele_types_ind', 'all_allele_types', 'num_all_allele_types', 'good_allele_inds', 'population', ...
                'count_vec', 'f_vec', 'n_vec', 'allele_inds_vec', 'allele_types', 'gene_by_allele_type_inds_list', 'GENE_INDS', 'ProteinPos');
            if(i_c==1) % first
                B = cur_B;
            else % next
                B = union_SFS_structs(B, cur_B, unite_field_names);
            end
            unite_file = i_c
        end % loop on different files in population
        B.good_allele_inds = get_good_allele_inds(B, exome_struct);
        B = internal_unite_by_class(B); % NEW! unite sub-classes into class
        save(fullfile(dir_from_file_name(sfs_file_names{1}), [exome_struct.prefix '_AllPop_union.mat']), '-struct', 'B'); % save union !!
    else
        B = load(fullfile(dir_from_file_name(sfs_file_names{1}), [exome_struct.prefix '_AllPop_union.mat']));
    end
else
    exome_struct.populations = '';     num_populations = 1;
end
load(mutation_rates_files); MutationTypes = MutationTypes(1:3); MutationRateTable = MutationRateTable(:,1:3);
if(exist('UniqueMutationRateTable', 'var'))
    UniqueMutationRateTable = UniqueMutationRateTable(:,1:3);
end
if(ischar(GeneStruct)) % load gene-struct file
    load(GeneStruct); % load gene-struct %     gene_struct_input_file = GeneStruct;
end
internal_plot_SFS_by_gene_position(B, GeneStruct);

% Compute theoretical constant population size distribution
N=10000; mu = mu_per_site; theta = 4*N*mu; % estimate for human effective population size and effective mutation rate
x_vec = (1:2*N-1)./(2*N); % allele frequency
constant_size_cum_phi_zero_vec = absorption_time_by_selection(0, 1, N, 1/(2*N), x_vec, 0);
constant_size_cum_phi_one_vec = absorption_time_by_selection(0, 1, N, 1/(2*N), x_vec, 'freq');
constant_size_cum_het_vec = absorption_time_by_selection(0, 1, N, 1/(2*N), x_vec, 'var');

total_num_alleles = zeros(num_populations, 1); fraction_allele_types = cell(num_populations, 1);
num_chr = max(B.XXX_REF_ALLELE_COUNT_ + B.XXX_VARIANT_COUNT_);
het_vec = cell(B.num_allele_types, num_populations); % vector of heterozygosities for each allele type
het_var_vec = cell(B.num_allele_types, num_populations); % used for confidence interval calculations
num_alleles_vec = cell(B.num_allele_types, num_populations);
for i=1:B.num_allele_types
    % % %     for j=1:num_populations
    % % %         f_vec{i,j} = A{j}.f_vec{i};
    % % %         count_vec{i,j} = A{j}.count_vec{i};
    % % %         het_vec{i,j} = 2 .* f_vec{i,j} .* (1-f_vec{i,j}); % multiply by two
    % % %         het_var_vec{i,j} = f_vec{i,j} .* (1-f_vec{i,j}) .* (1-2.*f_vec{i,j}).^2 ./ num_chr; % variance of heterozygosity (???)
    % % %         %    [h_freq bin_freq] = hist(f_vec{i}, num_bins);
    % % %         num_alleles_vec{i,j} = f_vec{i,j} .* f_vec{i,j}; % weight by allele frequencies
    % % %     end
    B.het_vec{i} = 2 .* B.f_vec{i} .* (1-B.f_vec{i}); % multiply by two
    B.het_var_vec{i} = B.f_vec{i} .* (1-B.f_vec{i}) .* (1-2.*B.f_vec{i}).^2 ./ repmat(num_chr, size(B.f_vec{i},1), 1); % variance of heterozygosity (???)
    B.num_alleles_vec{i} = B.f_vec{i} .* B.f_vec{i};
end
for j=1:num_populations
    total_num_alleles(j) = sum(cellfun('size', B.f_vec, 1));   % sum(length_cell(f_vec(:,j))); % determine the total number of alleles.
    fraction_allele_types{j} = cellfun('size', B.f_vec, 1) ./ total_num_alleles(j);   % length_cell(f_vec(:,j)) ./ total_num_alleles(j); % fraction of number of alleles for each classw
end

num_snps = length_cell(f_vec);
if(isfield(B, 'GENE')) % is this ESP data?
    B.num_genes = length(unique(B.GENE));
    %    for j=1:num_populations
    B.good_allele_inds{4} = B.good_allele_inds{4}(end:-1:1); % temp. hack - reverse NS vs. S
    %    end
else
    if(isfield(A, 'XXX_GENE_'))
        B.num_genes = length(unique(B.XXX_GENE_));
    end
end
good_allele_inds = B.good_allele_inds;

%new_A = cell(num_populations, 1);
alpha_fit = zeros(num_populations, 1); alpha_fit_by_freq = alpha_fit;
ratio_vec2 = zeros(num_populations,  B.num_allele_types);
[variants, carriers, singletons, heterozygosity] = ... % New: compute mutation rates per site
    compute_average_mutation_rates_per_class(B, B.f_vec, B.count_vec, ...
    B.het_vec, target_length, B.good_allele_inds);
%save('-append', spectrum_population_data_file{j}, 'het_vec', 'het_var_vec', 'good_allele_inds');
save('-append', sfs_file_names{1}, 'variants', 'carriers', 'singletons', 'heterozygosity');
new_B = struct('variants', variants, 'carriers', carriers, 'singletons', singletons, 'heterozygosity', heterozygosity);

for j=1:num_populations
    % % %     if(exist('spectrum_data_file', 'var')) % Save again to same file new fields
    % % %         tmp_het_vec = het_vec; het_vec = het_vec(:,j);
    % % %         tmp_het_var_vec = het_var_vec; het_var_vec = het_var_vec(:,j);
    % % %         save('-append', spectrum_population_data_file{j}, 'het_vec', 'het_var_vec', 'good_allele_inds');
    % % %         save('-append', spectrum_population_data_file{j}, 'variants', 'carriers', 'singletons', 'heterozygosity');
    % % %         new_A{j} = struct('variants', variants, 'carriers', carriers, ...
    % % %             'singletons', singletons, 'heterozygosity', heterozygosity);
    % % %         het_vec = tmp_het_vec; het_var_vec = tmp_het_var_vec;
    % % %     end
    % Fit alpha_0 (crude fit)
    alpha_fit(j) = ( heterozygosity.ratio_over_stop_gained_vec(j,B.good_allele_inds{5}(2)) - heterozygosity.ratio_over_stop_gained_vec(j,B.good_allele_inds{5}(3)) ) ./ ...
        ( heterozygosity.ratio_over_stop_gained_vec(j,B.good_allele_inds{5}(1)) - heterozygosity.ratio_over_stop_gained_vec(j,B.good_allele_inds{5}(3)) ); % fraction of missense which are roughly 'lethal'
    ratio_vec2(j,:) = (new_B.variants.per_gene ./ new_B.singletons.per_gene(j,:)) ./ ...
        (new_B.variants.per_gene(13) ./ new_B.singletons.per_gene(j,13)); % compute relative in a different way
    alpha_fit_by_freq(j) = ( ratio_vec2(j,B.good_allele_inds{5}(2)) - ratio_vec2(j,B.good_allele_inds{5}(3)) ) ./ ...
        ( ratio_vec2(j,B.good_allele_inds{5}(1)) - ratio_vec2(j,B.good_allele_inds{5}(3)) ); % fraction of missense which are roughly 'lethal'
    % alpha_fit = alpha_fit_by_freq;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Save also data as tab-delimited text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = compute_SFS_table(B, MutationRateTable, MutationTypes, exome_struct, output_file_name); % A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_type_vec = ...
    {'num_variants_cum', ... %
    'num_variants_cum_normalized', ... %
    'num_carriers_cum', ... %
    'heterozygosity_cum', ... %
    'heterozygosity_cum_normalized_log', ... %
    'singletons_hist', ... %
    'singletons_heterozygosity_hist', ... %
    'allele_freq_hist', ... %
    'heterozygosity_hist', ... %
    'heterozygosity_hist_zoom', ... %
    'num_alleles_per_gene_hist', ... %
    'heterozygosity_per_gene_hist', ... %
    'fraction_of_null_alleles', ... %
    'singletons_per_gene', ... %
    'heterozygosity_per_gene', ... %
    'heterozygosity_per_site', ... %
    'enrichment_missense_hist'};

bin_size = min(diff(num_bins));

new_fig_pop = 1; new_fig_allele_type = 1;

for figure_type = { ... % 'enrichment_missense_hist', ...  % figure_type_vec
        'num_variants_cum_log_x_normalized', 'num_variants_cum', ...
        'num_carriers_cum_normalized', ... % 'num_carriers_cum', 'heterozygosity_cum',
        'heterozygosity_cum_normalized'} % , 'fraction_of_null_alleles'} %   [0 1 1.5 1.8] %  3 7 9 10 11 12 13] % [2 3 7 8] % try different ways of plotting
    [normalized_flag, plot_x_vec, plot_y_vec, additional_plot_x_vec, additional_plot_y_vec, ...
        x_str, y_str, my_x_lim, my_y_lim, legend_vec, legend_loc, fig_str] = ...
        internal_plot_SFS_figure(figure_type, B, new_B, num_populations, ...
        x_vec, B.f_vec, B.het_vec, B.num_alleles_vec, ...
        constant_size_cum_phi_zero_vec, constant_size_cum_phi_one_vec, constant_size_cum_het_vec, ...
        GeneStruct, MutationTypes, UniqueMutationRateTable, theta, ...
        alpha_fit);
    
    if(normalized_flag) % normalize y figure
        for i=1:length(plot_y_vec)
            plot_y_vec{i} = bsxfun(@rdivide, plot_y_vec{i}, plot_y_vec{i}(end));
        end
        for i=1:length(additional_plot_y_vec)
            additional_plot_y_vec{i} = bsxfun(@rdivide, additional_plot_y_vec{i}, additional_plot_y_vec{i}(end));
        end
    end
    
    log_x_flag = strfind(figure_type{1}, 'log_x'); log_y_flag = strfind(figure_type{1}, 'log_y');
    if(log_x_flag)     % Determine figure type
        if(log_y_flag) % log-log
            plot_str = 'loglog';
        else  % log-linear
            plot_str = 'semilogx';
        end
    else
        if(log_y_flag) % linear-log
            plot_str = 'semilogy';
        else  % linear-linear
            plot_str = 'plot';
        end
    end
    if(exist('plot_x_vec', 'var') && (~isempty(plot_x_vec)))
        for j=1:num_populations % here just make sub-plots.
            if(new_fig_pop)
                subplot(2, 4, j); j_sim=1;  % figure;
            else
                j_sim = j;
            end
        end
        ctr=1;
        %             if(new_fig_allele_type)
        %                 figure;
        %             end
        for i=1:length(B.good_allele_inds{5})
            for j=1:num_populations % here plotting is done !!
                if(new_fig_pop)
                    subplot(2, 4, j); j_sim=1;  % figure;
                else
                    j_sim = j;
                end
                eval(['h(' num2str(ctr) ') = ' plot_str '(plot_x_vec{' num2str(ctr) '}, plot_y_vec{' num2str(ctr) ...
                    '}, ''' symbol_vec{j_sim} ''', ''linewidth'', 2, ''color'', ''' color_vec(i) ''');']); hold on; ctr=ctr+1;
                if(new_fig_pop)
                    title(exome_struct.populations{j});
                end
                if(i == length(B.good_allele_inds{5})) % last plot
                    add_faint_grid(0.5, 0); % add grid
                end
            end % loop on j (population)
        end
    end % loop on populations and ..
    if(exist('additional_plot_x_vec', 'var') && (~isempty(additional_plot_x_vec)))
        for i=1:length(additional_plot_x_vec)
            eval(['h(' num2str(ctr) ') = ' plot_str '(additional_plot_x_vec{' num2str(i) '}, additional_plot_y_vec{' num2str(i) ...
                '}, ''linewidth'', 2, ''color'', ''' color_vec(i+length(B.good_allele_inds{5})) ''');']); hold on;
            ctr=ctr+1;
        end
    end
    xlabel(x_str); ylabel(y_str);
    if(~isempty(my_x_lim))
        xlim(my_x_lim);
    end
    if(~isempty(my_y_lim))
        ylim(my_y_lim);
    end
    
    if(~new_fig_pop) % one title at end
        if(ismember(figure_type, 'heterozygosity_hist_zoom')) % [7])) % something weird is wrong with legend for 8 (zoom-in)
            legend(legend_vec, 'location', legend_loc); % h
        else
            %        legend(h([1:num_populations:num_populations*length(A{j}.good_allele_inds{5})  ctr-1]), legend_vec, ...
            legend(legend_vec, 'location', legend_loc, 'fontsize', 14, 'fontweight', 'bold'); % h
        end
    else % here legend for subplots
        legend({'Synonymous', 'Missense', 'Loss-of-Function'});
    end
    
    legend boxoff;
    if( (iscell(legend_vec)) && (~isempty(strfind(lower(legend_vec{1}), 'nons'))) )
        tmp_str = {'NS', 'S', 'STOP'};
    else
        tmp_str = {'S', 'NS', 'STOP'};
    end
    
    if(~new_fig_pop) % one title at end
        title(str2title([str2word('.', remove_dir_from_file_name(output_file_name), 'end') ...
            ' ' fig_str ', 2n\_s=' num2str(num_chr)  ', #GENES=' num2str(A{1}.num_genes) ', #SNPs=(' tmp_str{1} ' ' ...
            num2str(num_snps(B.good_allele_inds{5}(1))) ', ' tmp_str{2} ' ' num2str(num_snps(B.good_allele_inds{5}(2))) ...
            ', ' tmp_str{3} ' ' num2str(num_snps(B.good_allele_inds{5}(3))) ')'] ), 'fontsize', 8);
    end
    % ...
    %    '), het. per gene=(' tmp_str{1} ' ' num2str(heterozygosity.per_gene(A{j}.good_allele_inds{5}(1)),2) ', ' tmp_str{2} ' ' ...
    %    num2str(heterozygosity.per_gene(A{j}.good_allele_inds{5}(2)),2) ') ' ...
    %    ' het. per site=(' tmp_str{1} ' ' num2str(heterozygosity.per_site(A{j}.good_allele_inds{5}(1)),2) ', ' tmp_str{2} ' ' ...
    %    num2str(heterozygosity.per_site(A{j}.good_allele_inds{5}(2)),2) ')']), 'fontsize', 8);
    my_saveas(gcf, [output_file_name '_'  cell2vec(exome_struct.populations, '-') '_all_' fig_str], {'epsc', 'jpg'}); % {'epsc', 'pdf', 'jpg', 'fig'});
end % loop on figure types


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Start internal functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Internal function for plotting constraint vs. gene location
%
% Input:
% A - structure with inofrmation about alleles
% Output:
%
function internal_plot_SFS_by_gene_position(B, GeneStruct) % gene_position_vec, s_missense_vec)
AssignGeneralConstants; AssignRVASConstants;  exome_data_figs_dir;
have_protein_pos_inds = find(~isempty_cell(B.ProteinPos)); % find alleles with protein position

[gene_names, I, J] = intersect(upper(B.unique_genes), upper(GeneStruct.gene_names));
gene_lens = GeneStruct.gene_lens(J); gene_strand = GeneStruct.strand(J);
B.gene_lens = zeros(length(B.GENE_INDS), 1);
[aaa, bbb] = ismember(B.GENE_INDS, I);
B.gene_lens(aaa>0) = gene_lens(bbb(aaa>0));
B.gene_strand = zeros(length(B.GENE_INDS), 1); % take strand
B.gene_strand(aaa>0) = gene_strand(bbb(aaa>0));


smooth_window_vec = [7500 12000 2000];
figure;
for i=1:3 % loop on synonymous, missense, stop 
    B.allele_types(B.good_allele_inds{5})
    cur_f_vec = B.f_vec{B.good_allele_inds{5}(i)}(:,1);
    cur_pos_vec = B.ProteinPos(B.allele_inds_vec{B.good_allele_inds{5}(i)}(:,1)); % get positions
    gene_lens_vec = B.gene_lens(B.allele_inds_vec{B.good_allele_inds{5}(i)});
    gene_strand_vec = B.gene_strand(B.allele_inds_vec{B.good_allele_inds{5}(i)});
    have_protein_pos_inds = find(~isempty_cell(cur_pos_vec) & (gene_lens_vec>0)); % find alleles with protein position
    cur_pos_vec = 3*cell2mat(cur_pos_vec(have_protein_pos_inds)) ./ gene_lens_vec(have_protein_pos_inds); % move from nucleotides to amino acids and relative position
    cur_pos_vec(gene_strand_vec(have_protein_pos_inds) == 1) = 1-cur_pos_vec(gene_strand_vec(have_protein_pos_inds) == 1); % Flip strand !!! 
    [sorted_cur_pos_vec, sort_perm] = sort( cur_pos_vec );
    sorted_cur_f_vec = cur_f_vec(have_protein_pos_inds); sorted_cur_f_vec = sorted_cur_f_vec(sort_perm);
%    sorted_gene_lens_vec = gene_lens_vec(have_protein_pos_inds); sorted_gene_lens_vec = sorted_gene_lens_vec(sort_perm); 
    num_plot = sum(sorted_cur_f_vec>0)
    plot(sorted_cur_pos_vec(sorted_cur_f_vec>0), smooth(sorted_cur_f_vec(sorted_cur_f_vec>0), smooth_window_vec(i)), [color_vec(i)], 'linewidth', 2); hold on;
end
xlabel('Protein Pos.'); ylabel('Allele Freq.'); legend({'Synonymous', 'Missense', 'Loss-of-Function'}); legend('boxoff'); 
x_lim = xlim(gca); xlim([0 1]); %  x_lim(2)]);
%add_faint_grid(0,5, 0); 
my_saveas(gcf, fullfile(exome_data_figs_dir, 'allele_freq_by_protein_pos'), {'epsc', 'jpg', 'pdf'}); 

% take europe

% Extract gene position:
gene_position_vec = [];
s_missense_vec = [];

figure; hold on; plot(gene_position_vec, s_missense_vec);
xlabel('Gene Position'); ylabel('Missense s');



% Internal function for plotting different SFS figures
%
% Input:
% ... (lots of variables. Put in structrs)
% Output:
% ...
function [normalized_flag, plot_x_vec, plot_y_vec, additional_plot_x_vec, additional_plot_y_vec, ...
    x_str, y_str, my_x_lim, my_y_lim, legend_vec, legend_loc, fig_str] = ...
    internal_plot_SFS_figure(figure_type, A, new_A, num_populations, ...
    x_vec, f_vec, het_vec, num_alleles_vec, ...
    constant_size_cum_phi_zero_vec, constant_size_cum_phi_one_vec, constant_size_cum_het_vec, ...
    GeneStruct, MutationTypes, UniqueMutationRateTable, theta, ...
    alpha_fit)

new_fig_pop = 1; new_fig_allele_type = 0;
Assign24MammalsGlobalConstants;

normalized_flag = strfind(figure_type{1}, 'normalized'); % Get plots
% cumulative_flag = strfind(figure_type{1}, 'cum');
clean_figure_type = strdiff(strdiff(strdiff( figure_type{1}, '_log_x'), '_log_y'), '_normalized');
my_x_lim = []; my_y_lim = [];

figure; ctr=1; plot_x_vec = cell(num_populations*length(A.good_allele_inds{5}), 1); plot_y_vec = plot_x_vec; legend_vec = plot_x_vec;
for i=vec2row(A.good_allele_inds{5}) % loop on different allele types 1:min(6, A.num_allele_types)
    if(new_fig_allele_type)
        figure;
    end
    for j=1:num_populations % loop on different populations
        cur_allele_type = A.allele_types{i}; %j % good_allele_inds{4}(i)); % get string of allele type. Mis-match!!!
        [sorted_f_vec, sort_perm] = sort(f_vec{i}(:,j));
        sorted_het_vec = het_vec{i}(sort_perm,j);
        if(strncmp('all_stop', cur_allele_type, 8))  %        if(i == good_allele_inds{4}(1)) % stop codons
            [stop_f_vec, I] = unique(sorted_f_vec);
        end
        if(strncmp('all_synonymous', cur_allele_type, length('all_synonymous')))  % if(i == good_allele_inds{4}(3)) % synonymous
            [synonymous_f_vec, I] = unique(sorted_f_vec);
        end
        switch clean_figure_type
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'enrichment_missense_hist' % plot enrichment/depletion of missense vs. synonymous variants for each gene
                legend_loc = 'northeast';
                switch cur_allele_type
                    case 'all_stop' % A{j}.allele_types(A{j}.good_allele_inds{3})  % {'stop', 'stop-gained'}
                end
                %                 enrichment_vec = A{1}.total_heterozygosity_per_gene_mat{1}(i,:); % get total heterozygocity
                %                 [x_vec h_vec] = hist( enrichment_vec, num_bins); % histogram of enrichment vectors
                %
                %                 plot(x_vec, h_vec);
                x_str = 'log_2(h_{missense} / h_{synonymous})'; y_str = '# genes';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'num_variants_cum' % 0 % plot cumulative allele frequencies (per-site, un-normalized)
                plot_x_vec{ctr} = sorted_f_vec;
                plot_y_vec{ctr} = (1-new_A.variants.per_site(j,i)) + ...
                    new_A.variants.per_site(j,i) .* (1:length(f_vec{i})) ./ length(f_vec{i});
                if(ismember(i, A.good_allele_inds{5}(3))) % strncmp('stop', cur_allele_type, 4)) %% if(i == good_allele_inds{4}(1)) % stop codons
                    stop_hist = (1-new_A.variants.per_site(j,i)) + ...
                        new_A.variants.per_site(j,i) .* (1:length(f_vec{i})) ./ length(f_vec{i});
                end
                if(ismember(i, A.good_allele_inds{5}(1))) % strncmp('coding-synonymous', cur_allele_type, length('coding-synonymous')))  %% if(i == good_allele_inds{4}(3)) % synonymous
                    synonymous_hist = (1-new_A.variants.per_site(j,i)) + ...
                        new_A.variants.per_site(j,i) .* (1:length(f_vec{i})) ./ length(f_vec{i});
                end
                x_str = 'Derived Allele Freq.'; y_str = '# Alleles Per-Site (Cumulative)';
                legend_loc = 'southeast';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'num_carriers_cum' % 1.5 % plot cumulative allele frequency weighted by # carriers (frequency)
                sorted_num_alleles_vec = num_alleles_vec{i}(sort_perm,j);
                switch cur_allele_type % i % save distirbution
                    case A.allele_types(A.good_allele_inds{5}(3))  % {'stop', 'stop-gained'} % good_allele_inds{4}(1) % stop codons
                        stop_hist = new_A.carriers.per_site(j,i) .* ...
                            cumsum(sorted_num_alleles_vec) ./ sum(sorted_num_alleles_vec);
                    case A.allele_types(A.good_allele_inds{5}(2)) % {'missense'} % good_allele_inds{4}(2) % missense
                        missense_hist = new_A.carriers.per_site(j,i) .* ...
                            cumsum(sorted_num_alleles_vec) ./ sum(sorted_num_alleles_vec);
                    case A.allele_types(A.good_allele_inds{5}(1)) % 'coding-synonymous' % good_allele_inds{4}(3) % synonymous
                        synonymous_hist = new_A.carriers.per_site(j,i) .* ...
                            cumsum(sorted_num_alleles_vec) ./ sum(sorted_num_alleles_vec);
                end % switch cur_allele_type
                plot_x_vec{ctr} = sorted_num_alleles_vec;
                plot_y_vec{ctr} = new_A.carriers.per_site(j,i) .* ...
                    cumsum(sorted_num_alleles_vec) ./ sum(sorted_num_alleles_vec); %
                legend_loc = 'southeast';
                x_str = 'Derived Allele Freq.'; y_str = 'Num. Alleles Per-Site (Cumulative)';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'heterozygosity_cum' % 2 % plot cumulative heterozygosity per-site (un-normalized)
                plot_x_vec{ctr} = sorted_f_vec;
                plot_y_vec{ctr} = new_A.heterozygosity.per_site(j,i) .* cumsum(sorted_het_vec) ./ sum(sorted_het_vec);
                switch cur_allele_type              %   if(i == good_allele_inds{4}(1)) % stop codons
                    case A.allele_types(A.good_allele_inds{5}(3)) % {'stop', 'stop-gained'}
                        stop_hist = new_A.heterozygosity.per_site(j,i) .* cumsum(sorted_het_vec) ./ sum(sorted_het_vec);
                    case A.allele_types(A.good_allele_inds{5}(1)) % 'coding-synonymous'  %               if(i == good_allele_inds{4}(3)) % synonymous
                        synonymous_hist = new_A.heterozygosity.per_site(j,i) .* cumsum(sorted_het_vec) ./ sum(sorted_het_vec);
                end
                legend_loc = 'southeast';
                x_str = 'Derived Allele Freq.'; y_str = 'Heterozygosity Per-Site (Cumulative)';
                %               ylim([0 1]); % force limits
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'singletons_hist' % 4 % plot singletons
                tail_bins = 1:100;
                h_tail = hist(count_vec{i}, tail_bins);
                plot_x_vec{ctr} = tail_bins+0.2*(ctr-1);
                plot_y_vec{ctr} = h_tail ./ sum(h_tail); % bar(.., 0.2, color_vec(ctr)); % 1:length(f_vec(i))); % Plot the tail (singletons etc.)
                legend_loc = 'northeast';
                x_str = '# Carriers'; y_str = 'Frac. of # Alleles'; % title(allele_types{i});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'singletons_heterozygosity_hist' % 5 % plot histogram at singletons but with var. explained
                h_tail_het = zeros(length(tail_bins), 1);
                for k=1:length(tail_bins)
                    if(k < length(tail_bins))
                        cur_inds = find(count_vec{i} == tail_bins(k));
                    else
                        cur_inds = find(count_vec{i} >= tail_bins(k));
                    end
                    h_tail_het(k) = sum(het_vec{i}(cur_inds));
                end
                plot_x_vec{ctr} = tail_bins+0.2*(ctr-1);
                plot_y_vec{ctr} = h_tail_het ./ sum(h_tail_het); %               bar(tail_bins+0.2*(ctr-1), h_tail_het ./ sum(h_tail_het), 0.2, color_vec(ctr)); % 1:length(f_vec(i))); % Plot the tail (singletons etc.)
                legend_loc = 'northeast';
                x_str = '# Carriers'; y_str = 'Heterozygosity'; % title(allele_types{i});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'allele_freq_hist' % 6 % Plot histogram of # alleles
                [h_freq, bin_freq] = hist(f_vec{i}, num_bins);
                plot_x_vec{ctr} = bin_freq+0.2*(ctr-1)*bin_size; % the 0.2 is optional (for bar plots)
                plot_y_vec{ctr} = h_freq ./ sum(h_freq);
                % bar(bin_freq+0.2*(ctr-1)*bin_size, h_freq ./ sum(h_freq), color_vec(ctr)); plot(bin_freq, h_freq ./ sum(h_freq), color_vec(ctr), 'linewidth', 2);
                x_str = 'Derived Allele Freq.'; y_str =  'Frac. of # Alleles';
                legend_loc = 'northeast';
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'heterozygosity_hist' % 7 % plot histogram of heterozygosity
                [h_het, bin_het] = weighted_hist(f_vec{i, j}, het_vec{i, j}, num_bins); %  min(count_vec{i}));
                h_het = normalize_hist(bin_het, h_het) .* new_B.heterozygosity.per_site(j,i); % normalize and multiply by total area
                save_total_het = sum(h_het);                 save_total_het = 1;
                [h_freq, bin_freq] = hist(f_vec{i}, num_bins); % count # of alleles in each bin
                h_het_std = h_het ./ ( save_total_het .* sqrt(h_freq(1:length(h_het))) ); % take coefficient of variation ..
                
                plot_x_vec{ctr} =     bin_het+0.4*(ctr-1)*bin_size;
                plot_y_vec{ctr} = h_het ./ save_total_het;
                % h(ctr) = bar(bin_het+0.4*(ctr-1)*bin_size, h_het ./ save_total_het, 0.2, color_vec(ctr)); hold on; %                   h(ctr) = plot(bin_het, h_het ./ save_total_het, color_vec(ctr), 'linewidth', 2);
                h2 = plot(bin_het, 2 .* new_B.heterozygosity.per_site(j,i) .* (1-bin_het), [color_vec(ctr) '--'], 'linewidth', 2); % plot fitted line
                if(i==B.good_allele_inds{5}(1)) % plot once, at the first allele (synonymous)
                    h2 = plot(bin_het, 2 .* theta .* (1-bin_het), 'k--', 'linewidth', 2); %
                end
                x_str = 'Derived Allele Freq.'; y_str =  'Heterozygosity per site';
                legend_loc = 'northeast';
                my_x_lim = [0, 1];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'heterozygosity_hist_zoom' % 8 % plot histogram of heterozygosity bins at the very low end
                eps = 0.1;  % contraction factor (we look at DAF < eps)
                %                    rare_inds = find(f_vec{i} < max(num_bins .* eps));
                frac_rare = sum(het_vec{i}(f_vec{i} < max(num_bins .* eps))) / sum(het_vec{i});
                [h_het, bin_het] = weighted_hist(min(f_vec{i}, max(num_bins .* eps)), het_vec{i}, num_bins .* eps); %  min(count_vec{i}));
                h_het = h_het(1:end-1); bin_het = bin_het(1:end-1);
                h_het = normalize_hist(bin_het, h_het) .* new_A.heterozygosity.per_site(j,i) .* frac_rare; % normalize and multiply by total area
                [h_freq, bin_freq] = hist(f_vec{i}, num_bins .* num_bins); % count # of alleles in each bin
                h_het_std = h_het ./ ( save_total_het .* sqrt(h_freq(1:length(h_het))) ); % take coefficient of variation ..
                save_total_het = 1;
                %              h2 = errorbar(bin_het+0.4*(ctr-1)*bin_size.*eps, h_het ./ save_total_het, h_het_std, color_vec(ctr), 'linestyle', 'none'); hold on;
                h(ctr) = bar(bin_het+0.3*(ctr-1)*bin_size.*eps, h_het ./ save_total_het, 0.3, color_vec(ctr));
                h2 = plot(bin_het, 2 .* new_A.heterozygosity.per_site(j,i) .* (1-bin_het), [color_vec(ctr) '--'], 'linewidth', 2); % plot fitted line
                if(i==B.good_allele_inds{5}(1)) % plot once (synonymous)
                    N=10000; mu = mu_per_site;
                    theta = 4*N*mu;
                    h2 = plot(bin_het, 2 .* theta .* (1-bin_het), 'k--', 'linewidth', 2); %
                end
                x_str = 'Derived Allele Freq.'; y_str =  'Heterozygosity per site';
                legend_loc = 'northeast';
                my_x_lim = [0, eps];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'num_alleles_per_gene_hist' % 9 % plot histogram of number of alleles per gene
                hist_density(A{1}.num_alleles_per_gene_mat(i,:), 100, color_vec(ctr));
                x_str = 'Num. Alleles'; y_str =  '# Genes';
                fig_str = '_num_allele_genes_hist'; legend_loc = 'northeast';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'heterozygosity_per_gene_hist' % 10 % plot histogram of total heterozygosity per gene
                % [tmp_h tmp_bins] = hist(A{1}.total_heterozygosity_per_gene_mat(i,:), 100);
                hist_density(A{1}.total_heterozygosity_per_gene_mat(i,:), 100, color_vec(ctr));
                x_str = 'Total Heterozygosity'; y_str =  '# Genes';
                legend_loc = 'northeast';
        end % switch figure type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fig_str = figure_type{1};
        switch clean_figure_type
            case {'num_variants_cum', 'heterozygosity_cum'} % {0,2}
                switch cur_allele_type
                    case A.allele_types(A.good_allele_inds{5}(3)) % {'stop', 'stop-gained'} %                 if(i == good_allele_inds{4}(1)) % stop codons
                        stop_hist = vec2column(stop_hist(I)); % get unique
                    case A.allele_types(A.good_allele_inds{5}(1)) %   {'coding-synonymous', ??} %                if(i == good_allele_inds{4}(3)) % synonymous
                        synonymous_hist = vec2column(synonymous_hist(I)); % get unique
                end % switch cur_allele_type
                if(i == A.good_allele_inds{5}(end))  % finished last allele type
                    [missense_fit_bins, missense_fit_hist, stop_hist_interp, synonymous_hist_interp] = ...
                        sum_hist(stop_f_vec, alpha_fit(j).*stop_hist, ...
                        synonymous_f_vec, (1-alpha_fit(j)).*synonymous_hist, 10, 0);
                    fit_inds = intersect( find(missense_fit_bins < min( max(stop_f_vec), max(synonymous_f_vec) )), ...
                        find(missense_fit_bins > max( min(stop_f_vec), min(synonymous_f_vec) )) );
                    missense_fit_bins = missense_fit_bins(fit_inds); missense_fit_hist = missense_fit_hist(fit_inds);
                    stop_hist_interp = stop_hist_interp(fit_inds); synonymous_hist_interp = synonymous_hist_interp(fit_inds);
                end % if last allele type
        end % switch figure type again ..
        legend_vec{ctr} = [strdiff(A.allele_types{i}, 'coding-') ', ' A.population ...
            ', het.=' num2str(new_A.heterozygosity.per_gene(j,i),2)];
        %         num_points = length(plot_x_vec); num_pos_points = sum(plot_x_vec > 0);
        %         legend_vec{ctr} = [legend_vec{ctr} '(' num2str(num_pos_points) ', ' num2str(num_points) ')']; % add #points to legend
        ctr=ctr+1;
    end % loop on different populations
end % loop on different allele types
%    legend_vec = vec2row(allele_types(A{j}.good_allele_inds{5})); % why only last population ??
%    for i=1:length(A{j}.good_allele_inds{5})       % add total heterozygosity
%        legend_vec{i} = [legend_vec{i} ', het.=' num2str(new_A{j}.heterozygosity.per_gene(A{j}.good_allele_inds{5}(i)),2)];
%    end

%stop_ind = find(A{1}.allele_types_ind == STOP); %target_stop_ind = find(MutationTypes == STOP);
missense_ind = find(A.allele_types_ind == MISSENSE); %target_missense_ind = (MutationTypes == MISSENSE);
synonymous_ind = find(A.allele_types_ind == SYNONYMOUS); %target_synonymous_ind = (MutationTypes == SYNONYMOUS); % Get different mutation types
[~, I_genes, J_genes] = intersect(A.unique_genes, GeneStruct.gene_names);

additional_plot_x_vec = []; additional_plot_y_vec = [];
switch clean_figure_type % additional plots specific to each type of plot
    case 'enrichment_missense_hist' % plot enrichment/depletion of missense vs. synonymous variants for each gene
        num_het_bins = 100;
        heterozygosity_h_vec = zeros(length(A.upper_freq_vec), num_het_bins); heterozygosity_x_vec = heterozygosity_h_vec;
        freq_h_vec = heterozygosity_h_vec; freq_x_vec = heterozygosity_h_vec; num_alleles_h_vec = heterozygosity_h_vec; num_alleles_x_vec = heterozygosity_h_vec;
        for k=1:length(A.upper_freq_vec)
            target_enrichment_vec = log2( UniqueMutationRateTable(:, (MutationTypes == MISSENSE)) ./ ...
                UniqueMutationRateTable(:, (MutationTypes == SYNONYMOUS)) );
            heterozygosity_enrichment_vec = log2( A.total_heterozygosity_per_gene_mat{k}(missense_ind,:) ./ ...
                A.total_heterozygosity_per_gene_mat{k}(synonymous_ind,:) ); % get total heterozygocity
            freq_enrichment_vec = log2( A.total_freq_per_gene_mat{k}(missense_ind,:) ./ ...
                A.total_freq_per_gene_mat{k}(synonymous_ind,:) ); % get total heterozygocity
            num_alleles_enrichment_vec = log2( A.num_alleles_per_gene_mat{k}(missense_ind,:) ./ ...
                A.num_alleles_per_gene_mat{k}(synonymous_ind,:) ); % get total heterozygocity
            
            heterozygosity_ratio_enrichment_vec = heterozygosity_enrichment_vec(I_genes)' - target_enrichment_vec(J_genes);
            freq_ratio_enrichment_vec = freq_enrichment_vec(I_genes)' - target_enrichment_vec(J_genes);
            num_alleles_ratio_enrichment_vec = num_alleles_enrichment_vec(I_genes)' - target_enrichment_vec(J_genes);
            [heterozygosity_h_vec(k,:), heterozygosity_x_vec(k,:) ] = ...
                hist_density( heterozygosity_ratio_enrichment_vec, num_het_bins, [], 0, A.num_genes); % histogram of enrichment vectors
            [freq_h_vec(k,:), freq_x_vec(k,:) ] = ...
                hist_density( freq_ratio_enrichment_vec, num_het_bins, [], 0, A.num_genes); % histogram of enrichment vectors
            [num_alleles_h_vec(k,:), num_alleles_x_vec(k,:) ] = ...
                hist_density( num_alleles_ratio_enrichment_vec, num_het_bins, [], 0, A.num_genes); % histogram of enrichment vectors
        end
        plot(heterozygosity_x_vec(end,:)', heterozygosity_h_vec(end,:)', 'linewidth', 2); hold on;
        plot(freq_x_vec(end,:)', freq_h_vec(end,:)', 'linewidth', 2, 'linestyle', '--');
        plot(num_alleles_x_vec(end,:)', num_alleles_h_vec(end,:)', 'linewidth', 2, 'linestyle', ':');
        line([0 0], [0 max(num_alleles_x_vec(:))], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
        legend_vec = num2str(vec2column(A.upper_freq_vec)); %            legend_vec = {'log2(missense/synonymous)'};
        legend_vec = [repmat('DAF<', length(A.upper_freq_vec), 1) legend_vec];
        
    case 'num_variants_cum' % 0 % plot cumulative of fitted distribution for missense (why needed??)
        % % %         additional_plot_x_vec = {x_vec}; % missense_fit_bins,
        % % %         additional_plot_y_vec = {constant_size_cum_phi_zero_vec ./ constant_size_cum_phi_zero_vec(end)}; % missense_fit_hist,
        % % %         legend_vec{end+1} = 'neutral, const. N'; %  'missense-fit',
        % % %         if(normalized_flag)
        % % %             my_y_lim = [0.98 1];
        % % %         end
    case 'fraction_of_null_alleles' % 1.8 % 0.5 plot alpha (fraction of nulls ??) as function of allele frequency
        stop_hist_interp = stop_hist_interp ./ alpha_fit;
        synonymous_hist_interp = synonymous_hist_interp ./ (1-alpha_fit);
        %             alpha_vec = alpha_fit .* ( stop_hist_interp - (1-variants.per_site(1)) ) ./ ...
        %                 ( alpha_fit .* ( stop_hist_interp - (1-variants.per_site(1)) )  + ...
        %                 (1-alpha_fit) .* ( synonymous_hist_interp - (1-variants.per_site(3)) ) );
        
        %             alpha_vec = alpha_fit .* ( stop_hist_interp - stop_hist_interp(1) ) ./ ...
        %                 ( missense_fit_hist -  missense_fit_hist(1) ); % use alpha0 best fit
        alpha_vec = ( (missense_fit_hist - synonymous_hist_interp) ./ (stop_hist_interp - synonymous_hist_interp) ) .* ...
            ( stop_hist_interp - stop_hist_interp(1) ) ./ ( missense_fit_hist -  missense_fit_hist(1) ); % use fit for each x
        
        bayes_alpha_vec = ...  % compute bayes-factor. Doesn't depend on differences in birth rates
            ( (stop_hist_interp - (1-new_A.variants.per_site(j,A.good_allele_inds{5}(1)))) ./ new_A.variants.per_site(j,A.good_allele_inds{5}(1)) ) ./ ...
            ( ( (stop_hist_interp - (1-new_A.variants.per_site(j,A.good_allele_inds{5}(1)))) ./ new_A.variants.per_site(j,A.good_allele_inds{5}(1)) )  + ...
            ( synonymous_hist_interp - (1-new_A.variants.per_site(j,A.good_allele_inds{5}(3))) ) ./ new_A.variants.per_site(j,A.good_allele_inds{5}(3)) );
        additional_plot_x_vec{1} = missense_fit_bins;
        additional_plot_y_vec{1} = alpha_vec;
        %            plot(missense_fit_bins, alpha_vec, 'linewidth', 2);
        legend_vec = '\alpha_s(f)'; legend_loc = 'northeast';  % just one legend
    case 'num_carriers_cum' % 1.5 not normalized
        additional_plot_x_vec{1} = x_vec;
        additional_plot_y_vec{1} = theta.*constant_size_cum_phi_one_vec ./ constant_size_cum_phi_one_vec(end);
        %            plot(x_vec, theta.*constant_size_cum_phi_one_vec ./ constant_size_cum_phi_one_vec(end), 'k--', 'linewidth', 2); % plot the theoretical constant population size distribution
        legend_vec{end+1} = 'neutral, const. N';
        my_x_lim = [10^(-5) 1];
    case 'heterozygosity_cum' % 2 % plot cumulative of fitted distribution
        additional_plot_x_vec = {x_vec}; % missense_fit_bins,
        additional_plot_y_vec = {constant_size_cum_het_vec ./ constant_size_cum_het_vec(end)}; % missense_fit_hist,
        legend_vec{end+1} = 'neutral, const. N'; % 'missense-fit',
    case 'singletons_per_gene' % 11 % ??? Plot # singletons per gene
        ratio_vec = singletons.per_gene ./ singletons.per_gene(13);
        normalized_ratio_vec = ratio_vec ./ sum(ratio_vec(A.good_allele_inds{5}));
        additional_plot_x_vec{1} = 1:length(A.good_allele_inds{5});
        additional_plot_y_vec{1} = singletons.per_gene(A.good_allele_inds{5});
        %            bar(singletons.per_gene(A{j}.good_allele_inds{5}));
        set(gca, 'xtick', [1:length(A.good_allele_inds{5})]);  x_str = '';
        set(gca, 'XTicklabel', allele_types(A.good_allele_inds{5})); y_str = '# singletons per gene';
        title(['# singletons per gene. Ratio: (stop,missense,synom.) ' num2str(ratio_vec(A.good_allele_inds{5}), 3)]);
    case 'heterozygosity_per_gene' % 12 % plot just a bar showing average heterozygosity per gene
        ratio_vec = heterozygosity.per_gene ./ heterozygosity.per_gene(13);
        additional_plot_x_vec{1} = 1:length(A.good_allele_inds{5});
        additional_plot_y_vec{1} = heterozygosity.per_gene(A.good_allele_inds{5});
        %            bar(heterozygosity.per_gene(A{j}.good_allele_inds{5}));
        set(gca, 'xtick', [1:length(A.good_allele_inds{5})]);  x_str = '';
        set(gca, 'XTicklabel', allele_types(A.good_allele_inds{5})); y_str = 'Heterozygosity per gene';
        title(['Heterozygosity per gene. Ratio: (stop,missense,synom.) ' num2str(ratio_vec(A.good_allele_inds{5}), 3)]);
    case 'heterozygosity_per_site' % 13 % plot heterozygosity per-site (in the Target!)
        ratio_vec = (heterozygosity.per_gene ./ singletons.per_gene) ./ ...
            (heterozygosity.per_gene(13) ./ singletons.per_gene(13));
        additional_plot_x_vec{1} = 1:length(A.good_allele_inds{5});
        additional_plot_y_vec{1} = heterozygosity.per_gene(A.good_allele_inds{5}) ./ singletons.per_gene(A.good_allele_inds{5});
        %            bar(heterozygosity.per_gene(A{j}.good_allele_inds{5}) ./ singletons.per_gene(A{j}.good_allele_inds{5}));
        set(gca, 'xtick', [1:length(A.good_allele_inds{5})]); x_str = '';
        set(gca, 'XTicklabel', allele_types(A.good_allele_inds{5})); y_str = 'Heterozygosity per site (in target)';
        title(['Heterozygosity per site (in target). Ratio: (stop,missense,synom.) ' num2str(ratio_vec(A.good_allele_inds{5}), 3)]);
end % switch figure type
add_faint_grid(0.5, 0); % add grid







% Internal function for computing rates for different mutations
%
% Input:
% A - structure with allele counts information
% f_vec - frequency of derived alleles for each allele
% count_vec - # of derived allele carriers for each allele
% het_vec - # heterozygisity for each allele
% num_genes - total number of genes with data
% target_length - total length in nucleotides of target sequenced
% good_allele_inds - which allele types do we want to plot
%
% Output:
% variants - structure indicating the number of UNIQUE variants per gene and per site
% carriers - structure indicating the number of alleles (non-UNIQUE, in the POPULATION) per gene and per site
% singletons - structure indicating the number of singletons per gene and per site
% heterozygosity - structure indicating the heterozygosity per gene and per site
%
function [variants, carriers, singletons, heterozygosity] = ... % Compute mutation rates per site
    compute_average_mutation_rates_per_class(A, f_vec, count_vec, het_vec, target_length, good_allele_inds)

num_pop = size(f_vec{1}, 2);

% Compute averages per-gene and per-site in target
variants.per_gene = cellfun('size', f_vec, 1) ./ A.num_genes; % here count each variant only once
singletons.per_gene = zeros(num_pop, A.num_allele_types);
for i=1:A.num_allele_types
    for j=1:num_pop
        singletons.per_gene(j,i) = sum(count_vec{i}(:,j) == 1);
    end
end
singletons.per_gene = singletons.per_gene ./ repmat(A.num_genes, num_pop, 1);
fraction_singleton_allele_types = singletons.per_gene ./ sum(singletons.per_gene);  % fraction of number of alleles for each class at birth
singletons.per_site = ( singletons.per_gene .* A.num_genes ./ (target_length .* fraction_singleton_allele_types) ) .* sum(fraction_singleton_allele_types(good_allele_inds{4}));
variants.per_site = ( cellfun('size', f_vec, 1) ./ (target_length .* fraction_singleton_allele_types) ) .* sum(fraction_singleton_allele_types(good_allele_inds{4}));
sum_f = sum_cell(f_vec, [], 1); sum_het = sum_cell(het_vec, [], 1);
carriers.per_gene = sum_f ./ A.num_genes; % here count alleles with multiplicty (this is per-one individual)
carriers.per_site = ( sum_f' ./ (target_length .* fraction_singleton_allele_types) ) .* ...
    repmat(sum(fraction_singleton_allele_types(:,good_allele_inds{4}), 2), 1, length(f_vec));
heterozygosity.per_gene = sum_het'./ A.num_genes; % compute heterozygosity per gene
heterozygosity.per_site = sum_het' ./ (target_length .* fraction_singleton_allele_types); % fraction_allele_types);
heterozygosity.ratio_over_stop_gained_vec = (heterozygosity.per_gene ./ singletons.per_gene) ./ ...
    repmat(heterozygosity.per_gene(:,13) ./ singletons.per_gene(:,13), 1, length(f_vec)); % compute relative heterozygosity for each allele type compared to stop-gained alleles




% Unite different sub-class annotations into three major classes: stop, missense and synonymous
% Input:
% A - strcuture with alleles information
% Output:
% A_new - new structure, with three major classes:  stop, missense and synonymous
function A_new = internal_unite_by_class(A)

A.good_allele_inds{5} = (A.num_allele_types+1):(A.num_allele_types+3);
A.num_allele_types = A.num_allele_types+3;

A.allele_types = [A.allele_types' {'all_synonymous', 'all_missense', 'all_stop'}]';
for i=1:3 % loop on class types
    A.f_vec{A.good_allele_inds{5}(i)} = cell2vec(A.f_vec(A.good_allele_inds{i}));
    A.n_vec{A.good_allele_inds{5}(i)} = cell2vec(A.n_vec(A.good_allele_inds{i}));
    A.count_vec{A.good_allele_inds{5}(i)} = cell2vec(A.count_vec(A.good_allele_inds{i}));
    A.allele_inds_vec{A.good_allele_inds{5}(i)} = cell2vec(A.allele_inds_vec(A.good_allele_inds{i}));
end
A_new=A; % copy output
