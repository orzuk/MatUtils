% Script reading rare-alleles table and computing variance explained for each study
% function parse_rare_alleles()
run_rare_allele_computation = 1;
rare_alleles_file_name = '../../common_disease_model/data/rare/RareVariants_2011_03_07.txt';
%rare_alleles_results_file_name = [remove_suffix_from_file_name(rare_alleles_file_name) '_results.txt'];
rare_alleles_results_file_name = fullfile(dir_from_file_name(rare_alleles_file_name), ...
    ['RareVariants_' datestr(now, 'yyyy-mm-dd') '_results.txt']);

alpha_vec = [0.3 0.7 1]; % assumed fraction of rare variants which are functional (1 is good for the case where one already filtered for functional variants)

if(run_rare_allele_computation) % here run all the stuff on rare alleles
    if(~exist(file_name_to_mat(rare_alleles_file_name), 'file'))
        rare = ReadDataFile(rare_alleles_file_name, file_name_to_mat(rare_alleles_file_name), 1, [], 9); % don't convert cell to mat
    else
        rare = load(file_name_to_mat(rare_alleles_file_name));
    end
end
full_inds = find(~isempty_cell(rare.Reference));
rare = struct_by_inds(rare, full_inds);

rare.num_genes(rare.num_genes == 21) = 3; % Ahituv et al. (localized effect)
rare.Cases = cell2mat(empty_cell_to_numeric_val(str2num_cell(rare.Cases),0));
rare.Controls = cell2mat(empty_cell_to_numeric_val(str2num_cell(rare.Controls),0));
rare.num_rare_variants_cases = cell2mat(empty_cell_to_numeric_val(str2num_cell(rare.num_rare_variants_cases),0));
rare.num_rare_variants_controls = cell2mat(empty_cell_to_numeric_val(str2num_cell(rare.num_rare_variants_controls),0));

f_range = 0.5 .* [rare.num_rare_variants_controls ./ rare.Controls ... % consider two chromosomes
    rare.num_rare_variants_cases ./ rare.Cases]; % determine minimum and maximum possible values for f

rare.min_f_per_gene = min(f_range,[],2) ./ rare.num_genes;
rare.max_f_per_gene = max(f_range,[],2) ./ rare.num_genes; % get a range for f
num_traits = length(rare.Trait);
f = zeros(num_traits,length(alpha_vec)); beta = zeros(num_traits,length(alpha_vec)); alpha_lowerbound = zeros(num_traits,1);
V = zeros(num_traits,length(alpha_vec));
for i=1:num_traits % 16:16 % Loop over all variants and compute var. explained 9:9 %
    if( any(i == [7 8 9 10 11] ) )
        continue % bug in some of the trait
    end
    if(i == 11) % temp stop for debugging
        xxxxx = 13242134
    end
    run_trait = i
    if( isempty(rare.R(i)) || (rare.R(i) == 0))
        rare.R(i) = (rare.num_rare_variants_cases(i) / rare.Cases(i)) ./ ...
            (rare.num_rare_variants_controls(i) / rare.Controls(i));
    end
    if(~isfield(rare, 'beta') || rare.beta(i) == 0)
        for j=1:length(alpha_vec) % loop on different alphas
            switch rare.Type{i}
                case {'QTL', 'Quantitative'} % Use Shamil's formula
                    [beta(i,j) f(i,j) V(i,j) alpha_lowerbound(i)] = ... %                         [B F VV] = ...
                        ratio_QTL_to_var_explained(rare.Controls(i)*2, rare.Cases(i)*2, ...
                        rare.num_rare_variants_controls(i), rare.num_rare_variants_cases(i), ...
                        rare.t1(i), rare.t2(i), alpha_vec(j), 0.09); % USe shamil's formula
                    if(alpha_lowerbound(i) > alpha_vec(j)) % mark as impossible
                        V(i,j) = 0/0;
                    end
                    % Define the prevalence based on quantile
                    rare.prevalence(i) = 1-rare.t2(i);
                    V_binary(i,j) = ratio_to_var_explained(rare.R(i), rare.prevalence(i)); % simple grr computation
                case {'Binary', 'binary'}
                    V(i,j) = ratio_to_var_explained(rare.R(i), rare.prevalence(i)); % simple grr computation
            end
            
        end
    else % here we already know beta
        beta(i,:) = rare.beta(i); % the 1:2 refers to two values of alpha?
        f(i,:) = rare.f(i);
        V(i,:) = beta_to_variance_explained(beta(i,:), f(i,:), 1, 'diploid');
    end
end
rare.prevalence = num2str_cell(num2cell(rare.prevalence));
rare.prevalence = strrep_cell(rare.prevalence, '0', '-');
rare = rmfield(rare, 'prevalence'); % don't display for now

rare.Comments = empty_cell_to_empty_str(rare.Comments);
rare.Beta = num2cell(rare.Beta);
rare.CombinedAlleleFrequency = num2cell(rare.CombinedAlleleFrequency);

for i=1:length(alpha_vec) % new: loop on different alpha
    eval_str = ['rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = V(:,' num2str(i) ');'];
    eval(eval_str);
    eval_str = ['rare.Var_Explained_per_gene_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = ' ...
        'rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' ./ rare.num_genes;' ];
    eval(eval_str);
end
%rare.Var__explained___lower_bound_ = V(:,2);
rare.alpha_lowerbound = alpha_lowerbound; % get minimal possible value of alpha
rare.Var__explained__Case_Control_Modeling = V_binary(:,2); % alternative computation: on disease scale

for i=1:num_traits
    rare.Beta{i} = ['(' cell2vec(num2str_cell(num2cell(beta(i,:)), 2), ' , ') ')'];
    rare.CombinedAlleleFrequency{i} = ['(' cell2vec(num2str_cell(num2cell(f(i,:))), ' , ') ')'];
end
good_inds = intersect(find(isfinite(rare.Var__explained__upper_bound_)), ...
    find(isfinite(rare.R)));
rare = struct_by_inds(rare, good_inds);
rare.R = max(rare.R, 1./rare.R); % take always the high effect size

%rare.Var__explained___lower_bound_alpha_0_7 = rare.Var__explained___lower_bound_;
%rare.Var__explained___upper_bound_alpha_0_3 = rare.Var__explained__upper_bound_;
rare = rmfield(rare, {'Var__explained__upper_bound_', 'Var__explained___lower_bound_'});
rare = rmfield(rare, {'Var__Explained_by_Common', 'Heritability', 'Comments'})

rare.R = num2str_cell(rare.R, 2);
rare.t1 = num2str_cell(rare.t1.*100, 2,[],1);
rare.t2 = num2str_cell(rare.t2.*100, 2,[],1);
rare.alpha_lowerbound = num2str_cell(rare.alpha_lowerbound.*100, 2,[],1);
rare.min_f_per_gene = num2str_cell(rare.min_f_per_gene.*100, 2,[],1);
rare.max_f_per_gene = num2str_cell(rare.max_f_per_gene.*100, 2,[],1);
%rare.Var__explained__Case_Control_Modeling = num2str_cell(rare.Var__explained__Case_Control_Modeling.*100, 2,[],1);


for i=1:length(alpha_vec)
    eval_str = ['rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = ' ...
        'num2str_cell(rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') '.*100, 2,[],1)'];
    eval(eval_str);
    eval_str = ['rare.Var_Explained_per_gene_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = ' ...
        'num2str_cell(rare.Var_Explained_per_gene_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') '.*100, 2,[],1)'];
    eval(eval_str);
    
    % replace 0 with '-'
    eval_str = ['rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = ' ...
        'strrep_cell(rare.Var_Explained_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_'), ...
        ', ''NaN%'', ''-'');'];
    eval(eval_str);
    eval_str = ['rare.Var_Explained_per_gene_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_') ' = ' ...
        'strrep_cell(rare.Var_Explained_per_gene_Alpha_' strrep(num2str(alpha_vec(i),2), '.', '_'), ...
        ', ''NaN%'', ''-'');'];
    eval(eval_str);
end
%rare.Var__explained___lower_bound_alpha_0_7 = num2str_cell(rare.Var__explained___lower_bound_alpha_0_7.*100, 2,[],1);
%rare.Var__explained___upper_bound_alpha_0_3 = num2str_cell(rare.Var__explained___upper_bound_alpha_0_3.*100, 2,[],1);

%rare.f = f;  %just copy input
rare = rmfield(rare, 'Var__explained__Case_Control_Modeling');
WriteDataFile(rare, rare_alleles_results_file_name,0); % Save results







%%% New: parse the new T2D paper - here we've got functional data 
%%  parse_T2D_MTNR1B()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



