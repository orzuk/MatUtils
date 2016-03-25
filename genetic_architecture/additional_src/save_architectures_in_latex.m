% Save good architectures in latex format
% 
% Input: 
% good_architectures - structure with good architectures parameters 
% good_architectures_latex_file - output latex file 
% good_architectures_file - output .mat file 
% good_architectures_plot_file - output figures file for  
% disease_type_str - type of model 
% 
function save_architectures_in_latex(good_architectures, ...
    good_architectures_latex_file, ...
    good_architectures_file, good_architectures_plot_file, disease_type_str) % save to nice latex tables format

precision = 3;
if(~exist('disease_type_str', 'var') || isempty(disease_type_str))
    disease_type_str = ''; % 'ALL-POSSIBLE-DISEASES'; % need something that will not appear in any file string
end

good_architectures_latex_file = ... % Use only '/' in latex (the backslash '\' doesn't work)
    strrep(good_architectures_latex_file, '\', '/');
good_architectures_file = ...
    strrep(good_architectures_file, '\', '/');
good_architectures_plot_file = ...
    strrep(good_architectures_plot_file, '\', '/');
good_architectures_plot_file = ...
    ['../../figures/' remove_dir_from_file_name(good_architectures_plot_file)]; % why force this ??? 
good_architectures_dir = '../../';

num_architectures = length(good_architectures);
if(num_architectures > 0) % not empty
    for i=1:num_architectures % shorten arch names
        good_architectures(i).plot_name_arch = arch_name_to_plot_name(good_architectures(i).arch);
%         fi = strfind(good_architectures(i).plot_name_arch, '(');
%         if(~isempty(fi))
%             good_architectures(i).plot_name_arch = good_architectures(i).plot_name_arch(1:fi-1);
%         end
        good_architectures(i).arch = strdiff(good_architectures(i).arch, 'moids');
    end
    
    num_fields_in_table = good_architectures(1).num_fields_in_table;
    % r_add = mat_into_vec(r_add); % output everything to a latex file
    num_good_architectures = length(good_architectures); num_fields = length(fieldnames(good_architectures));
    latex_tab = reshape(struct2cell(good_architectures), num_fields, num_good_architectures)';
    latex_tab = latex_tab(:,1:num_fields_in_table); % remove some fields which do not appear in the table
    
    latex_str_title = fieldnames(good_architectures); latex_str_title = latex_str_title(1:num_fields_in_table);
    for i=1:length(latex_str_title)
        if(~isempty(strfind(latex_str_title{i}, '_')))
            f = strfind(latex_str_title{i}, '_');
            latex_str_title{i} = ['$' latex_str_title{i}(1:f) '{' ...
                latex_str_title{i}(f+1:end) '}$'];
        end
    end
    latex_str_title = ['$\#$' latex_str_title']'; % first index is arch #
    interesting_vec = zeros(size(latex_tab, 1),1);
    for i=1:size(latex_tab, 1) % highlight significant ratios/ interesting architectures
        interesting_vec(i) = good_architectures(i).interesting_flag;
        if(good_architectures(i).interesting_flag)
            for j=1:size(latex_tab, 2)
                latex_tab{i,j} = ['{ \bf ' num2str(latex_tab{i,j}, precision) ' }'];
            end
        end
    end
    latex_tab = latex_tab(find(interesting_vec),:); %  = % Display only interesting ones
    latex_tab = [num2str_cell(num2cell(1:size(latex_tab, 1)))' latex_tab];
    for i=1:size(latex_tab,1)
        latex_tab{i,1} = ['\htmlref{' latex_tab{i,1} '}{link:arch_' num2str(i) '}'];
    end
    
    latex_tab = [latex_str_title latex_tab']';
    latex_str = latex(latex_tab, 2, precision); % convert to latex format
    latex_str = mat2cell(latex_str, ones(size(latex_str,1),1));
    latex_str{1} = strrep(latex_str{1}, 'MAF', 'M.A.F.'); % not sure if its the first .. need to check
    for i=1:size(latex_str, 1)
        latex_str{i} = strrep(latex_str{i}, '\begin{bmatrix}', '');
        latex_str{i} = strrep(latex_str{i}, '\end{bmatrix}', '');
    end
    %    latex_str = [ latex_str' ]';
    tab_header = {'\begin{table}[h!] % put tables summarizing all architectures', ...
        '\begin{center}', '{\scriptsize', latex_str{1}};
    tab_footer = {'\end{tabular}', '}', '\caption{Architectures Summary Statistics}', ...
        '\label{tab:additive_vs_interaction_effects}', ...
        '\end{center}', '\end{table}', '', '\newpage', ''};
    latex_str = split_latex_table(latex_str(2:end-1), tab_header, tab_footer, 60); % remove end{tabular} (moved to header)
    latex_str = [latex_str' {'\clearpage'}]';
    
    % latex_str = [latex_str' {'', '', '', '', 'Architectures:'} mat_into_vec(best_architecture_formula)']'; % add the formulas
    my_mkdir(dir_from_file_name(good_architectures_latex_file)); 
    savecellfile(latex_str, good_architectures_latex_file, [], 1);
    
    
    minimal_fig = {'\begin{figure}[thb!]', '\begin{center}'}; % Add a figure on minimal lods-ratio
    minimal_fig{end+1} = ...
        ['\psfig{file=' good_architectures_plot_file ...
        '_minimal_lods_ratio.eps,width=10cm}'];
    minimal_fig{end+1} = ...
        ['\caption{ Minimal lods-ratio achieved ' ...
        '(while satisfying other constraint) for various architectures'];
    minimal_fig{end+1} = ...
        '\label{fig:common_disease_minimal_lods_ratio}}';
    minimal_fig{end+1} = '\end{center}';
    minimal_fig{end+1} = '\end{figure}';
    minimal_fig{end+1} = '\newpage';
    minimal_fig{end+1} = '';
    latex_str = [latex_str' {''} minimal_fig]'; % concatenate minimal figure
    
    
    power_fig = {'\begin{figure}[thb!]', '\begin{center}'}; % Add a figure on minimal lods-ratio
    power_fig{end+1} = ...
        ['\psfig{file=' good_architectures_plot_file ...
        '_sample_size_power_vs_lods_ratio_power.eps,width=12cm}'];
    power_fig{end+1} = ...
        ['\caption{ Arch. power - num samples needed for 0.1, 0.5, 0.9 ' ...
        'power. MAF=0.1 is comparable to Altchuler et al.'];
    power_fig{end+1} = ...
        '\label{fig:common_disease_all_arch_power}}';
    power_fig{end+1} = '\end{center}';
    power_fig{end+1} = '\end{figure}';
    %    power_fig{end+1} = '\newpage';
    power_fig{end+1} = '';
    latex_str = [latex_str' {''} power_fig]'; % concatenate power figure
    
    power_fig = {'\begin{figure}[thb!]', '\begin{center}'}; % Add a figure on minimal lods-ratio
    power_fig{end+1} = ...
        ['\psfig{file=' dir_from_file_name(good_architectures_plot_file) ...
        'altschuler_et_al_fig_2_power.eps,width=16cm,height=9cm}'];
    power_fig{end+1} = ...
        ['\caption{ Power from Altchuler et al.'];
    power_fig{end+1} = ...
        '\label{fig:common_disease_altchuler_power}}';
    power_fig{end+1} = '\end{center}';
    power_fig{end+1} = '\end{figure}';
    power_fig{end+1} = '\newpage';
    power_fig{end+1} = '';
    latex_str = [latex_str' {''} power_fig]'; % concatenate power figure
    
    max_N_display = 10; % maximal value that fits into screen without latex errors
    good_ctr=1;
    for i=1:length(good_architectures) % Prepare a 'card' for each good architecture
        if(good_architectures(i).interesting_flag) % make a card only for interesting ones
            good_architectures(i).card = [];
            good_architectures(i).card{1} = ...
                ['\subsubsection*{Architecture ' num2str(i) ...
                ':} \label{link:arch_' num2str(i) '}']; % table with parameters
            good_architectures(i).card{2} = ['{\scriptsize ' good_architectures(i).formula  ' } '];
            good_architectures(i).card(3:5) = {'\begin{table}[!h]', ...
                '\begin{center}', '{\scriptsize'};
            good_architectures(i).card(6:9) = mat2cell(latex(latex_tab([1 good_ctr+1],:), ...
                2, precision), ones(4,1)); % table parameter names
            good_architectures(i).card{8} = strrep(good_architectures(i).card{8}, ...
                ['\htmlref{' num2str(i) '}{link:arch_' num2str(i) '}'], ...
                ['\htmlref{' num2str(i) '}{split_1' ...
                '_tab:additive_vs_interaction_effects}']); % send hyperling back to big table
            
            good_architectures(i).card(10:13) = {'}', '\end{center}', '\end{table}', ''};
            good_architectures(i).card{14} = 'Inheritance interaction matrix $h_{ij}$ (Marginals $h_i$ on diagonal).';
            N = length(good_architectures(i).h_ij_full);
            N = min(N, max_N_display);
            tmp_c = mat2cell(good_architectures(i).h_ij_full(1:N,1:N), ones(N,1), ones(1,N));
            for j=1:N
                tmp_c{j,j} = ['{\bf ' num2str(tmp_c{j,j}) '}'];
            end
            if(N < length(good_architectures(i).h_ij_full))
                for j=1:N
                    tmp_c{N,j} = '...'; tmp_c{j,N} = '...';
                end
            end
            
            good_architectures(i).card{15} = '{\tiny';
            good_architectures(i).card(16:16+N-1) = ...
                strrep(mat2cell(latex(tmp_c), ones(N,1)), ...
                'NaN', '-'); % matrix of pairwise interactions     % good_architectures(i).h_ij_full
            good_architectures(i).card{16+N} = '}';
            good_architectures(i).card{17+N} = 'Log-ratio interaction matrix $L_{ij}$ (Marginals $L_i$ on diagonal).';
            tmp_c = mat2cell(good_architectures(i).L_ij_full(1:N,1:N), ones(N,1), ones(1,N));
            for j=1:N
                tmp_c{j,j} = ['{\bf ' num2str(tmp_c{j,j}) '}'];
            end
            if(N < length(good_architectures(i).h_ij_full))
                for j=1:N
                    tmp_c{N,j} = '...'; tmp_c{j,N} = '...';
                end
            end
            good_architectures(i).card{17+N+1} = '{\tiny';
            good_architectures(i).card(17+N+2:17+2*N+1) = ...
                strrep(mat2cell(latex(tmp_c), ones(N,1)), ...
                'NaN', '-'); % matrix of pairwise interactions     % good_architectures(i).L_ij_full
            good_architectures(i).card{end+1} = '}';
            
            %         latex_z_risk_expected_tab = ...
            %             vec_into_mat(good_architectures(i).mu_pairwise_expected, 2);
            
            latex_z_risk_marginal_tab = cell(2,3);
            latex_z_risk_marginal_tab{1,1} = '$x_i$';
            latex_z_risk_marginal_tab{1,2} = '0';
            latex_z_risk_marginal_tab{1,3} = '1';
            latex_z_risk_marginal_tab{2,2} = ...
                ['{\bf ' num2str(good_architectures(i).p_z_x_marginal(2) / ...
                sum(good_architectures(i).p_z_x_marginal(1:2)), precision) '} (' ...
                num2str(good_architectures(i).freq, precision) ')'];
            latex_z_risk_marginal_tab{2,3} = ...
                ['{\bf ' num2str(good_architectures(i).p_z_x_marginal(4) / ...
                sum(good_architectures(i).p_z_x_marginal(3:4)), precision) '} (' ...
                num2str(good_architectures(i).freq, precision) ')'];
            good_architectures(i).card{end+1} = '{\scriptsize';
            p_vec = pop_prob_to_case_control_prob(good_architectures(i).p_z_x_marginal);
            good_architectures(i).card{end+1} = ...
                ['Disease prob., 1st locus: $Pr(Z=1 | x_i)$. {\bf Observed} ' ...
                '(Expected based on no effect). Effect size: $\Delta=' ...
                num2str(abs(p_vec(1) - p_vec(2)), precision) '$']; % correct effect for case-control studies
            %                 num2str(abs(good_architectures(i).p_z_x_marginal(4) - ...
            %                 good_architectures(i).freq * good_architectures(i).f_vec(1)), precision) '$'];
            good_architectures(i).card(end+1:end+2) = ...
                mat2cell(latex(latex_z_risk_marginal_tab), ones(2,1));
            good_architectures(i).card{end+1} = '}';
            
            latex_z_risk_pairwise_tab = cell(3);  % an extra small table
            latex_z_risk_pairwise_tab{1,1} = '$x_i \backslash x_j$';
            latex_z_risk_pairwise_tab{1,2} = '0';
            latex_z_risk_pairwise_tab{1,3} = '1';
            latex_z_risk_pairwise_tab{2,1} = '0';
            latex_z_risk_pairwise_tab{3,1} = '1';
            ctr=1;
            for j=1:2
                for k=1:2
                    latex_z_risk_pairwise_tab{j+1,k+1} = ...
                        ['{\bf ' num2str(good_architectures(i).mu_pairwise(ctr), precision) ...
                        '} (' ...
                        num2str(good_architectures(i).mu_pairwise_expected(ctr), precision) ...
                        ')'];
                    ctr=ctr+1;
                end
            end
            %         latex_z_risk_pairwise_tab(2:3,2:3) = ...
            %             num2cell(vec_into_mat(good_architectures(i).mu_pairwise,2));
            good_architectures(i).card{end+1} = '{\scriptsize';
            
            % Fitting formula is: Pr(z = 1 | x_1, x_2) / Pr(z = 0 | x_1, x_2) =
            % log(alpha + alpha_1*x_1 + alpha_2*x_2.
            % This is the linear logistic model
            x_12_logit_matrix = [0 0; 0 1; 1 0; 1 1]; % const, x_1, x_2
            z_logit_vec = good_architectures(i).mu_pairwise; % Prob. of Z equal one
%             N_logit_vec = good_architectures(i).mu_pairwise ...
%                 + good_architectures(i).mu_pairwise; % counts of Z equals one and zero
            
            [b dev stats] =  glmfit(x_12_logit_matrix, ... % Model: log(Pr(z=1) / Pr(z=0)) = b_0 + b_1*x_1 + b_2*x_2
                [z_logit_vec ], 'binomial', 'link', 'logit');
            delta = abs(good_architectures(i).mu_pairwise(4) - ...
                exp(sum(b)) / (1 + exp(sum(b)))); % subtract fitted value from actual value
            
            
            good_architectures(i).card{end+1} = ...
                ['Disease prob., 1st interaction: $Pr(Z=1 | x_i, x_j)$. ' ...
                '{\bf Observed} (Expected based on no-epistasis). Effect size: ' ...
                '$\Delta=' ...
                num2str(delta, precision) '$'];
%                 num2str(abs((good_architectures(i).mu_pairwise(4) - ...
%                 good_architectures(i).mu_pairwise_expected(4)) * ...
%                 good_architectures(i).f_vec(1)^2), precision) '$'];
            
            good_architectures(i).card(end+1:end+3) = ...
                mat2cell(latex(latex_z_risk_pairwise_tab), ones(3,1));
            good_architectures(i).card{end+1} = '}';
            
            %        good_architecture(i).card{3} = latex_tab(i+1,:); % table with parameters
            latex_str = [latex_str' {''} good_architectures(i).card '\newpage']'; % concatenate all cards
            
            fig_str = {'_all_four_figs_', '_risk_given_k_ones_of_M_known_loci_', '_family_risk_'}; % '_prob_disease_given_k_ones_arch_', '_power_arch_', 
            caption_str = { ... % {'Prob. disease given k risk loci', 'Power to detect best locus', ...
                ['Various statistics for architecture. In top left we display the power to ' ...
                'detect the strongest loci in a GWAS - the first marginal effect and the ' ...
                'strongest pairwise effect beyond marginal. n is the total number of samples ' ...
                'in the study where we assume $n/2$ cases and $n/2$ controls. ' ...
                'In the top right we display the prob. of having the disease when exactly k ' ...
                'out of the N loci are present (set to one). In the bottom left we show ' ...
                'the prob. of having exactly k loci present (set to one). In the bottom ' ...
                'right we show the risk of disease for a relative of a person with the ' ...
                'disease (penetrance). The x-axis denotes the relative degree, where 0 ' ...
                'means twin (penet.), 1 means son, 2 grand-son and so on for further ' ...
                'generations (relative risk was computed with sampling, so a small error may ' ...
                'be expected - the exact value at zero should match the penetrance value and ' ...
                'at the far right should match the population frequency value).'], ...
                ['Disease risk for carriers of k loci out of M known risk loci, where M is assumed to ' ...
                'be either $N/2$ or $N/3$ '], ...
                ['Segregation of disease risk in families. Shown is relative risk (compared to population ' ...
                'frequency) of a relative of diseased individual (gray). Based on simulating $5000$ families. ' ...
                'Number in gray refers to Monozygous twins relative risk']};
            fig_cm_vec = [18 18 18]; 
            [cur_good_architectures_file cur_good_architectures_plot_file] = ...
                get_architecture_file_names(good_architectures_dir, ...
                [good_architectures(i).plot_name_arch '_' disease_type_str], good_architectures(i).N);
            cur_good_architectures_plot_file = ...
                strrep(cur_good_architectures_plot_file, '\', '/');
            if(~isfield(good_architectures(i), 'index'))
                good_architectures(i).index = i;
            end
            %             fig_name_vec = {good_architectures_plot_file, ...
%                 good_architectures_plot_file, good_architectures_plot_file};
            for j=1:3 % loop over three figures. First is prob(disease|k ones) and second is family risk for observed loci, 
                % and third is family seggregation (used to be power)
                good_architectures(i).fig = {'\begin{figure}[thb!]', '\begin{center}'}; % Prepare a 'figure' for each architecture
                good_architectures(i).fig{end+1} = ...
                    ['\psfig{file=' cur_good_architectures_plot_file ...
                    fig_str{j} num2str(good_architectures(i).index) '.eps,width=' num2str(fig_cm_vec(j)) 'cm}']; % change to 6-7 for single figures
                good_architectures(i).fig{end+1} = ...
                    ['\caption{ ' caption_str{j} ' (all in arch. ' num2str(i) ')'];
                good_architectures(i).fig{end+1} = ...
                    ['\label{fig:common_disease_'  fig_str{j} num2str(i) '}}'];
                good_architectures(i).fig{end+1} = '\end{center}';
                good_architectures(i).fig{end+1} = '\end{figure}';
                good_architectures(i).fig{end+1} = '';
                latex_str = [latex_str' {''} good_architectures(i).fig]'; % concatenate all cards
                latex_str = [latex_str' {'\clearpage', '\newpage'}]'; % clear page between figures
            end
            
            
            family_size = length(good_architectures(i).family_tree);
            gender_vec = repmat([0 1 1 0], 1, ceil(family_size/4)); % get female/male shapes (box/circle)
%            node_shapes = node_shapes(1:family_size);
            pedigree_dir = remove_suffix_from_file_name(good_architectures_file)
            my_mkdir(pedigree_dir); 
            traits_vec = zeros(family_size,1); traits_vec(end) = 1;
            pedigree_names_vec = num2str_cell(num2cell(good_architectures(i).family_risk ./ good_architectures(i).freq), 2); % use precision 2
            for j=1:length(pedigree_names_vec)
                pedigree_names_vec{j} = ['A'+j-1 ':' pedigree_names_vec{j}];
            end
            save_pedigree(good_architectures(i).family_tree, pedigree_names_vec, gender_vec, traits_vec, ...
                fullfile(pedigree_dir, ['pedigree_arch_' num2str(i) '.txt']), 'pedigraph')
            
            % New! save architecture's relative-risks as pedigrees            
            
            % % %
            % % %         \begin{figure}[thb!]
            % % % \begin{center}
            % % % \psfig{file=../figures/common_disease_prob_relatives.eps,width=6cm}
            % % % \caption{Disease probability for relatives of an effected person for a {\it common} disease.
            % % % All models were normalized such that the twin probability of disease $h$ is the same and was set to $0.7$.
            % % % The overall population disease probability $p$ is also the same for all models and stands at $0.1$.
            % % % Similar to the rare-disease case, the Mendelian model actually shows the slowest decay (here the decay is not
            % % % by a factor of two) but now there is less difference between models. The common-weak and rare-strong models still
            % % % look very similar.The protective model cannot fit the data here as well.
            % % % \label{fig:common_disease_prob_relatives}}
            % % % \end{center}
            % % % \end{figure}
            good_ctr = good_ctr+1;
        end % if good architecture
    end  % loop on architectures 
    
    latex_str = [latex_str'  {'', '\end{document}'}]';
    
    R = loadcellfile(fullfile(dir_from_file_name(good_architectures_latex_file), 'architectures_header.tex'));
    latex_str = [R' latex_str']';
    saving_latex_file = good_architectures_latex_file
    my_mkdir(dir_from_file_name(good_architectures_latex_file));
    savecellfile(latex_str, good_architectures_latex_file, [], 1); % save architectures in latex format
    saving_mat_file = file_name_to_mat(good_architectures_file)
    my_mkdir(dir_from_file_name(good_architectures_file));
    save(file_name_to_mat(good_architectures_file), 'good_architectures'); % save good architectures in .mat file
end % if not empty

