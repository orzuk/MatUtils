% Plot and save heritability deviation figs for QTL
%
% Input:
% h_pop - epidemiological estimates of h (from twin studies)
% h_pop_str - types of estimates
% h_z - genetic estimate (sum of additive contribution of all loci)
% qtl_R - correlation of relatives in trait
% operation - type of architecture
% same_inds_prob - probability that the same liability was chosen as maximal for two relatives
% N_vec - vector of different number of liabilities
% main_ind - index of ...
% h_x - ??
% h_x_vec - vector measuring heritability in each pathway 
% h_x_shared_env_vec - vector measuring the amount of shared environment for relatives in EACH PATHWAY  
% h_shared_env_vec - vector measuring the amount of shared environment for relatives in the trait  
% isoheritability_flag - should all models have the same heritability
% (on phenotype scale), or same 'raw' parameters
%
function  plot_heritability_deviation_figures_QTL(h_pop, h_pop_str, h_z, qtl_R, operation, ...
    same_inds_prob, N_vec, main_ind, h_x, h_x_vec, h_x_shared_env_vec, h_shared_env_vec, isoheritability_flag)

AssignGeneralConstants;
format_fig_vec = 'epsc';
%fig_outdir = '../../common_disease_model/docs/figs/new';
fig_outdir = '../../common_disease_model/docs/pnas/genetic_interactions/figs';
data_outdir = '../../common_disease_model/docs/pnas/genetic_interactions/tables';
visscher_table_file = '../../common_disease_model/data/epidemiology/visscher_QTL_twin_correlation_statistics.txt';

num_h_pop = length(h_pop);

shared_inv_inds = union(find(h_x_shared_env_vec == 0), find(h_x_shared_env_vec == 0.5)); %[1:length(h_x_shared_env_vec)]; % [1:2]; % drop the 100%
save_h_x_shared_env_vec = h_x_shared_env_vec;
save_h_z = h_z; h_z = h_z(:,shared_inv_inds,:);
save_qtl_R = qtl_R; qtl_R = qtl_R(:,shared_inv_inds);
h_x_shared_env_vec = h_x_shared_env_vec(shared_inv_inds);
save_h_pop = h_pop; 
for i_pop = 1:num_h_pop % loop on different estimators
    save_h_phantom_vec{i_pop} = (h_pop{i_pop} - save_h_z) ./ h_pop{i_pop};
    h_pop{i_pop} = h_pop{i_pop}(:,shared_inv_inds,:);
    h_phantom_vec{i_pop} = (h_pop{i_pop} - h_z) ./ h_pop{i_pop};
end
%h_phantom_ADE_vec = (h_pop_ADE - h_z) ./ h_pop_ADE;

num_h = length(h_x_vec); num_shared_env = length(h_x_shared_env_vec); num_N = length(N_vec);
qtl_r_overestimation_vec = zeros(num_h,num_shared_env,size(qtl_R{1,1},2));
r_mz_vec = zeros(num_h,num_shared_env,size(qtl_R{1,1},2));
r_dz_vec = zeros(num_h,num_shared_env,size(qtl_R{1,1},2));
h_legend_vec = num2str_cell(num2cell(h_x_vec'*100),3,[],1);

for i=1:num_h
    if(isoheritability_flag)
        h_legend_vec{i} = ['h_x^2=' h_legend_vec{i}];
    else
        h_legend_vec{i} = [h_legend_vec{i} '%'];
%        h_legend_vec{i} = ['h_{pathway}^2=' h_legend_vec{i}];
    end
    for j=1:num_shared_env
        qtl_r_overestimation_vec(i,j,:) = ...
            (qtl_R{i,j}(1,:)./2 - qtl_R{i,j}(2,:)) ./ ...
            (qtl_R{i,j}(1,:)./2); % how much is the sib-correlation underestimating heritability    %    h_phantom_vec; % need to change that!!!
        r_mz_vec(i,j,:) = qtl_R{i,j}(1,:);
        r_dz_vec(i,j,:) = qtl_R{i,j}(2,:);
    end
end
h_pop_MZ = r_mz_vec;
h_pop_DZ = 2*r_dz_vec;
h_phantom_MZ_vec = (h_pop_MZ - h_z) ./ h_pop_MZ;
h_phantom_DZ_vec = (h_pop_DZ - h_z) ./ h_pop_DZ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('num_shared_envs', 'var') ||isempty(num_shared_envs))
    num_shared_envs=length(h_x_shared_env_vec);
end
for i_pop = 1:num_h_pop % loop on different heritability esitmators
    full_figure; hold on; % plot over-estimation of heritability and sib-correlation
    shared_env_legend_vec = num2str_cell(num2cell(h_x_shared_env_vec'*100),3,[],1);
    for j=1:num_shared_envs
        shared_env_legend_vec{j} = ['V_c/V_e=' shared_env_legend_vec{j}];
        plot(100*reshape(h_phantom_vec{i_pop}(:,j,:), num_h, num_N)', ...
            100*reshape(qtl_r_overestimation_vec(:,j,:), num_h, num_N)', ...
            symbol_vec{1+0*j}, 'linewidth', 2);
        plot_colors = get(gca, 'colororder');
        
        for i=1:num_h
            hh = arrow([100*h_phantom_vec{i_pop}(i,j,end-1), 100*qtl_r_overestimation_vec(i,j,end-1)], ...
                [100*h_phantom_vec{i_pop}(i,j,end), 100*qtl_r_overestimation_vec(i,j,end)], ...
                'length', 30); % Need to color this!!!, ...
            set(hh, 'FaceColor', plot_colors(i,:));
            set(hh, 'EdgeColor', plot_colors(i,:));
        end
        
        if(j == num_shared_env)
            plot(100*mat2vec(h_phantom_vec{i_pop}(main_ind,j,:)), ...
                mat2vec(100*qtl_r_overestimation_vec(main_ind,j,:)), ...
                'ko', 'linewidth', 2);
            plot(100*mat2vec(h_phantom_vec{i_pop}(main_ind,j,:)), ...
                mat2vec(100*qtl_r_overestimation_vec(main_ind,j,:)), ...
                'k.', 'linewidth', 2);
            plot(100*mat2vec(h_phantom_vec{i_pop}(main_ind,j,:)), ...
                mat2vec(100*qtl_r_overestimation_vec(main_ind,j,:)), ...
                'k*', 'linewidth', 2);
        end
    end % loop on different environments
    for i=1:length(N_vec)
        for j=num_shared_env:num_shared_env % choose lowest h_x, highest h_shared
            text(2+100*h_phantom_vec{i_pop}(1,j,i) + 0.6*(sign(1.5-i)+1), ...
                max(-100, 100*qtl_r_overestimation_vec(1,j,i)), ...
                [num2str(N_vec(i))], 'fontsize', 10, 'fontweight', 'bold');
        end
        % %     text(2+1.5*length(num2str(N_vec(i)))+100*h_phantom_vec(main_ind,i) + 0.6*(sign(1.5-i)+1), ...
        % %         max(1, 100*qtl_r_overestimation_vec(main_ind,i)), ...
        % %         ['(h_x^2=' num2str(100*h_z(main_ind,i),3) '%, r_{MZ}=' ...
        % %         num2str(qtl_R{main_ind}(1,i)*100,3) '%, r_{DZ}=', ...
        % %         num2str(qtl_R{main_ind}(2,i)*100,3) '%)'], 'fontsize', 8, 'fontweight', 'light', ...
        % %         'color', [0.25 0.25 0.25]); % gray
        % %     %            , h_x^2=' num2str(100*h_x_output_mat(main_ind,i),3) '%)'], 'fontsize', 12, 'fontweight', 'bold');
    end
    title(['h^2 overestimation for ' operation ' architecture. ' ...
        'Common environment V_c/V_e=0% (solid), 50% (dotted), 100% (dashed)']);
    xlabel('heritability overestimation (%): \pi');
    ylabel('r_{DZ} underestimation  (%): 1- 2r_{DZ}/r_{MZ}');
    legend(h_legend_vec); % shared_env_legend_vec);
    
    xlim([-0.1 100]); ylim([-100 80]); %    xlim([0 100]); ylim([0 11]);
    my_saveas(gcf, ['../../common_disease_model/figs/qtl_heritability_deviations_curve_' operation '_main'], ...
        format_fig_vec); % plot main figure
end % loop on different heritability estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shared_env_legend_vec = num2str_cell(num2cell(h_x_shared_env_vec'*100),3,[],1);

R = ReadDataFile(visscher_table_file,[],0,1); % don't convert yet cell to mat !! 
R = struct_by_inds(R, 2:87); % take 86 traits
R.all_r_DZ = cell2mat(empty_cell_to_numeric_val(R.all_r_DZ,0));
R.all_r_MZ = cell2mat(empty_cell_to_numeric_val(R.all_r_MZ,0));

marker_size = [13 4 24 8];
for i_pop = 1:num_h_pop
    switch h_pop_str{i_pop}
        case 'ACE'
            plot_ind_vec = 0:3;
        otherwise
            plot_ind_vec = 4;
    end
    for plot_ind = plot_ind_vec; % 0:6 % different plots
        full_figure; hold on; % plot over-estimation of heritability and sib-correlation
        for j=1:num_shared_env
            shared_env_legend_vec{j} = ['V_c/V_e=' shared_env_legend_vec{j}];
            
            switch plot_ind
                case 0  %plot r_MZ vs. r_DZ
                    fig_outfile = fullfile(fig_outdir, 'fig1b');  title_str = 'fig1b';
                    xlabel('$r_{MZ} (%)$', 'Interpreter','latex'); 
                    ylabel('$r_{DZ} (%)$', 'Interpreter','latex');
                    x_vec = 100*reshape(r_mz_vec(:,j,:), num_h, num_N)';
                    y_vec = 100*reshape(r_dz_vec(:,j,:), num_h, num_N)';
                case 1 % plot h_pop vs. h_phantom
                    fig_outfile = fullfile(fig_outdir, 'fig1a'); title_str = 'fig1a';
                    xlabel('$h_{pop}^2 (\%)$', 'Interpreter','latex');
                    ylabel('$\pi_{phantom} (\%)$', 'Interpreter','latex');
                    
                    x_vec = 100*reshape(min(h_pop{i_pop}(:,j,:),1), num_h, num_N)';
                    y_vec = 100*reshape(h_phantom_vec{i_pop}(:,j,:), num_h, num_N)';
                case 2 % plot h_pop vs. h_loci
                    x_vec  = 100*reshape(min(h_pop{i_pop}(:,j,:),1), num_h, num_N)';
                    y_vec  = 100*reshape(h_z(:,j,:), num_h, num_N)';
                    fig_outfile = fullfile(fig_outdir, 'supp_fig2b');  title_str = 'supp_fig2b';
                    xlabel('$h_{pop}^2 (\%)$', 'Interpreter','latex');
                    ylabel('$h_{all}^2 (\%)$', 'Interpreter','latex'); % how much explained by loci
                case 3 % plot r_MZ vs. r_DZ overestimation
                    fig_outfile = fullfile(fig_outdir, 'supp_fig2a');  title_str = 'supp_fig2a';
                    x_vec = 100*reshape(r_mz_vec(:,j,:), num_h, num_N)';
                    y_vec = 100*reshape(2.*r_dz_vec(:,j,:)./r_mz_vec(:,j,:) -1 , num_h, num_N)';
                    xlabel('$r_{MZ} (%)$', 'Interpreter','latex'); 
                    ylabel('$2r_{DZ}/r_{MZ}-1 (%)$', 'Interpreter','latex');
                    xlim([0, 100]); ylim([-50, 100]); % temp to make all arrows on same scale
                case 4 % just like case 1, but with ADE heritability (instead of ACE)
                    fig_outfile = fullfile(fig_outdir, ['supp_fig2' h_pop_str{i_pop}]);
                    title_str = ['supp_fig2'  h_pop_str{i_pop}];
                    xlabel(['${h_{pop}}^2 ' h_pop_str{i_pop} '(\%)$'], 'Interpreter','latex');
                    ylabel(['phantom heritability ' h_pop_str{i_pop} ' (\%): $1- \frac{h^2}{{h_{pop}}^2}$'], 'Interpreter','latex');
                    
                    x_vec = 100*reshape(min(h_pop{i_pop}(:,j,:),1), num_h, num_N)';
                    y_vec = 100*reshape(h_phantom_vec{i_pop}(:,j,:), num_h, num_N)';
                case 5 % just like case 1, but with r_MZ heritability (instead of ACE)
                    fig_outfile = fullfile(fig_outdir, 'supp_fig2MZ'); title_str = 'supp_fig2MZ';
                    xlabel('${h_{pop}}^2 r_{MZ} (\%)$', 'Interpreter','latex');
                    ylabel('phantom heritability r_{MZ} (\%): $1- \frac{h^2}{{h_{pop}}^2}$', 'Interpreter','latex');
                    
                    x_vec = 100*reshape(min(h_pop_MZ(:,j,:),1), num_h, num_N)';
                    y_vec = 100*reshape(h_phantom_MZ_vec(:,j,:), num_h, num_N)';
                case 6 % just like case 1, but with r_DZ heritability (instead of ACE)
                    fig_outfile = fullfile(fig_outdir, 'supp_fig2DZ'); title_str = 'supp_fig2DZ';
                    xlabel('${h_{pop}}^2 r_{DZ} (\%)$', 'Interpreter','latex');
                    ylabel('phantom heritability r_{DZ} (\%): $1- \frac{h^2}{{h_{pop}}^2}$', 'Interpreter','latex');
                    
                    x_vec = 100*reshape(min(h_pop_DZ(:,j,:),1), num_h, num_N)';
                    y_vec = 100*reshape(h_phantom_DZ_vec(:,j,:), num_h, num_N)';
            end
            plot(x_vec, y_vec, symbol_vec{1+0*j}, 'linewidth', 1); % plot line
            plot(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(j-1)
                symbol_vec{12-5*(j-1)}, 'linewidth', 1, 'markersize', marker_size(j)+5); % plot all dots except first
            
            plot_colors = get(gca, 'colororder');
            for k=1:num_h % plot first dot (bigger)
                plot(x_vec(1,k), y_vec(1,k), ... % 12-7*(j-1)
                    symbol_vec{12-5*(j-1)}, 'color', plot_colors(k,:), 'linewidth', 1, 'markersize', marker_size(2+j));
            end            
            
            for i=1:num_h
                hh = arrow([x_vec(end-1,i), y_vec(end-1,i)], ...
                    [x_vec(end,i), y_vec(end,i)], ...
                    'length', 45); %  'BaseAngle', 0); % Need to color this!!!, ...
                %        hh = arrow([100*h_pop(i,j,end-1), 100*h_phantom_vec(i,j,end-1)], ...
                %            [100*h_pop(i,j,end), 100*h_phantom_vec(i,j,end)]); % Need to color this!!!, ...
                set(hh, 'EdgeColor', plot_colors(i,:));
                if(j == 1)  %no shared environment
                    set(hh, 'FaceColor', plot_colors(i,:));
                else
                    set(hh, 'FaceColor', 'w');
                end
                
            end
            if(j==2) % text exampled
                plot(x_vec(4,3), y_vec(4,3), '*r', 'markersize', 16);
            end
            
            
            %     if(j == num_shared_env)
            %         plot(100*mat2vec(h_pop(main_ind,j,:)), ...
            %             mat2vec(100*h_phantom_vec(main_ind,j,:)), ...
            %             'ko', 'linewidth', 2);
            %         plot(100*mat2vec(h_pop(main_ind,j,:)), ...
            %             mat2vec(100*h_phantom_vec(main_ind,j,:)), ...
            %             'k.', 'linewidth', 2);
            %         plot(100*mat2vec(h_pop(main_ind,j,:)), ...
            %             mat2vec(100*h_phantom_vec(main_ind,j,:)), ...
            %             'k*', 'linewidth', 2);
            %     end
        end % loop on different environments
        switch plot_ind
            case 0
                plot(100*R.all_r_MZ, 100*R.all_r_DZ, '.k', 'linewidth', 1, 'MarkerSize', 13);
                xlim([-0.1 100]); ylim([0 75]); %    xlim([0 100]); ylim([0 11]);
            case 1
                xlim([-0.1 100]); ylim([-1 100]); %    xlim([0 100]); ylim([0 11]);
            case 2
                plot(0:100, 0:100, '--k', 'linewidth', 1);
                xlim([-0.1 100]); ylim([0 91]); %    xlim([0 100]); ylim([0 11]);
            case 3
                plot(100*R.all_r_MZ, 100*(2.*R.all_r_DZ./R.all_r_MZ-1), '.k', 'linewidth', 1, 'MarkerSize', 13);
                xlim([-0.1 100]); % ylim([0 100]); %    xlim([0 100]); ylim([0 11]);
        end
        
        h_leg = legend(h_legend_vec, 2); % shared_env_legend_vec);
        legend boxoff;
        leg_pos = get(h_leg, 'position');
        leg_pos(1:2) = leg_pos(1:2) + [0.07 -0.03];
        set(h_leg, 'position', leg_pos);
        title([repmat(' ', 1, 100) str2title(title_str)], 'fontsize', 16, 'fontweight', 'bold');
        my_saveas(gcf, fig_outfile,     format_fig_vec); % plot main figure
    end % loop on plot ind
    close all;
end % loop on i_pop

h_x_shared_env_vec = save_h_x_shared_env_vec;
num_shared_env = length(h_x_shared_env_vec);
num_vals = num_h*num_shared_env*num_N;
h_z = save_h_z;
qtl_R = save_qtl_R;
r_mz_vec = zeros(num_h,num_shared_env,size(qtl_R{1,1},2));
r_dz_vec = zeros(num_h,num_shared_env,size(qtl_R{1,1},2));
for i=1:num_h
    for j=1:num_shared_env
        r_mz_vec(i,j,:) = qtl_R{i,j}(1,:);
        r_dz_vec(i,j,:) = qtl_R{i,j}(2,:);
    end
end

if(isempty(h_shared_env_vec)) % fill shared environment. Always assume that the first index has 0 shared environment
    h_shared_env_vec = zeros(num_h, num_shared_env, num_N); 
    for i=1:num_h
        for j=1:num_shared_env
            h_shared_env_vec(i,j,:) = r_mz_vec(i,j,:) - r_mz_vec(i,1,:); % subtract the case with no shared environment
        end
    end
end


S = 100*[0.01*vec2column(mat2vec(repmat(N_vec, num_h*num_shared_env, 1))) ...
    vec2column(repmat(h_x_vec, 1, num_shared_env*num_N)) ...
    vec2column(repmat(mat2vec(repmat(h_x_shared_env_vec, num_h , 1)), num_N, 1)) ... % c_R
    vec2column( repmat(1-h_x_vec, 1, num_shared_env*num_N)' .* ...
    repmat(mat2vec(repmat(h_x_shared_env_vec, num_h , 1)), num_N, 1) ) ... % h_{path,c}
    mat2vec(r_mz_vec) mat2vec(r_dz_vec) mat2vec(h_z) ...
    mat2vec(h_shared_env_vec)]; % New! add shared environment
S_header = {'$k$', '$h^2_{path} (\%)$', '$c_R (\%)$', '$V_{path,c} (\%)$',  ...
    '$r_{MZ} (\%)$', '$r_{DZ} (\%)$', '$h^2_{all} (\%)$', '$V_c (\%)$'}; % new! add overall shared environment !!! 
h_pop = save_h_pop; h_phantom_vec = save_h_phantom_vec;
for i_pop = [1 2 5] % 1:num_h_pop
    S = [S ...
        100*mat2vec(min(1,h_pop{i_pop})) 100*mat2vec(max(0,h_phantom_vec{i_pop}))]; % concatenate and set bounds [0,100] to all heritabilities
    S_header = [S_header ['$h^{2 (' h_pop_str{i_pop} ')} _{pop} (\%)$'] ...
        ['$\pi_{phan}^{(' h_pop_str{i_pop} ')} (\%)$']];
end

S = [S_header' num2cell(S)']';

savecellfile(S, fullfile(data_outdir, 'fig1_data_quantitative_traits_tab.txt'),[],1);  % Save figures data
S = S(:, [1:3 5:end]); % remove shared environment within the pathway 
for i=2:size(S,1)
    for j=4:size(S,2)
        S{i,j} = sprintf('%.1f', S{i,j});
    end
end
S_latex = latex(S, 2, precision);
S_latex = mat2cell(S_latex, ones(size(S_latex,1),1));
for i=2+(num_h*num_shared_env):(num_h*num_shared_env):length(S_latex)-num_N % add lines when k is changed
    S_latex{i} = [S_latex{i} '  \hdashline'];
end

S_latex{1} = strrep(S_latex{1},  '|c', '|r'); % align to right
tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', ...
    ['\caption[Model parameters for quantitative traits]{Model parameters for ' ...
    ' the LP model $\LP(k, \hx, c_R)$ for quantitative traits. }'], ...
    '{\tiny', S_latex{1}};
%    '{\scriptsize', S_latex{1}};
tab_footer = {'}', ...
    '\label{table:fig1_data_quantitative_traits}', ...
    '\end{center}', '\end{table}', '', '\newpage', ''};
% S_latex = [tab_header  S_latex(2:end)' tab_footer]';
S_latex = split_latex_table(S_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 60, -1); % 50 for scriptsize, 60 for tiny font
savecellfile(S_latex, fullfile(data_outdir, 'fig1_data_quantitative_traits.txt'),[],1);  % Save figures data

return;

% % % for i= [1:length(N_vec)-3 length(N_vec)]
% % %     for j=num_shared_env:num_shared_env % choose lowest h_x, highest h_shared
% % %         text(2+100*h_pop(end,j,i) + 0.6*(sign(1.5-i)+1), ...
% % %             max(-100, 100*h_phantom_vec(end,j,i)), ...
% % %             [num2str(N_vec(i))], 'fontsize', 10, 'fontweight', 'bold');
% % %     end
% % % % %     text(2+1.5*length(num2str(N_vec(i)))+100*h_phantom_vec(main_ind,i) + 0.6*(sign(1.5-i)+1), ...
% % % % %         max(1, 100*qtl_r_overestimation_vec(main_ind,i)), ...
% % % % %         ['(h_x^2=' num2str(100*h_z(main_ind,i),3) '%, r_{MZ}=' ...
% % % % %         num2str(qtl_R{main_ind}(1,i)*100,3) '%, r_{DZ}=', ...
% % % % %         num2str(qtl_R{main_ind}(2,i)*100,3) '%)'], 'fontsize', 8, 'fontweight', 'light', ...
% % % % %         'color', [0.25 0.25 0.25]); % gray
% % % % %     %            , h_x^2=' num2str(100*h_x_output_mat(main_ind,i),3) '%)'], 'fontsize', 12, 'fontweight', 'bold');
% % % end
% % % title(['Phantom heritability for ' operation ' architecture. ' ...
% % %     'Common environment V_c/V_e=0% (circle), 50% (star)']);
%['../../common_disease_model/figs/qtl_phantom_heritability_' operation '_main'], ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % full_figure; % plot over-estimation of heritability and of sib-correlation (why plot again?? supp. info.)
% % num_h = length(h_x_vec);
% % plot(100*h_phantom_vec(1:min(num_h,7),:)', 100*qtl_r_overestimation_vec(1:min(num_h,7),:)', 'linewidth', 2);
% % if(num_h >= 8)
% %     plot(100*h_phantom_vec(8,:)', 100*qtl_r_overestimation_vec(8,:)', 'g', 'linewidth', 2);
% % end
% % if(num_h >= 9)
% %     plot(100*h_phantom_vec(9,:)', 100*qtl_r_overestimation_vec(9,:)', 'm', 'linewidth', 2);
% % end
% % plot(100*h_phantom_vec, 100*qtl_r_overestimation_vec, 'k.', 'linewidth', 2);
% % % %     for i=1:length(h_x_vec)
% % % %         for j=1:length(N_vec)
% % % %             text(0.5+100*h_phantom_vec(i,j), 100*qtl_r_overestimation_vec(i,j), ...
% % % %                 [num2str(N_vec(j))], 'fontsize', 12, 'fontweight', 'bold');
% % % %             %' (h_{pop}^2=' num2str(100*h_pop(i,j),3) ...
% % % %             %    '%, h_{loci}^2=' num2str(100*h_z(i,j),3) '%, h_x^2=' num2str(100*h_x_output_mat(i,j),3) '%)'], 'fontsize', 12, 'fontweight', 'bold');
% % % %         end
% % % %     end
% % title(['Overestimation of heritability for non-linear ' operation ' genetic architecture']);
% % xlabel('h^2 overestimation (%): (h_{pop}^2 - h_{loci}^2) / h_{pop}^2');
% % ylabel('r_s underestimation  (%): (r_s^{add} - r_s^{obs})/r_s^{add}');
% % %    xlim([0 100]); ylim([0 20]);
% % xlim([0 100]); ylim([0 90]);
% % legend_vec = num2str_cell(num2cell(100*h_x_vec),[],[], 1);
% % for i=1:length(legend_vec)
% %     legend_vec{i} = ['h_{pop}^2=' legend_vec{i}];
% % end
% % legend(legend_vec, 2);
% % my_saveas(gcf, ['../../common_disease_model/figs/qtl_heritability_deviations_curve_' operation], ...
% %     format_fig_vec); % save main figure
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % full_figure; % plot mz-correlation, sib-correlation and heritability explained by all loci
% % % % % for j=1:num_shared_env
% % % % %     plot(N_vec, qtl_R{main_ind,j}(1,:).^2, color_vec(j)); % plot r_MZ^2
% % % % %     plot(N_vec, qtl_R{main_ind,j}(2,:).^2, [color_vec(j) symbol_vec{1+0*(j+1)}]); % plot r_DZ^2
% % % % %     plot(N_vec, mat2vec(h_pop(main_ind,j,:)), 'g'); % plot h_pop^2
% % % % %     plot(N_vec, mat2vec(h_z(main_ind,j,:)), 'm'); % plot h_loci^2
% % % % %     %        plot(N_vec, h_z2, 'm--'); % plot h_loci^2 computed differently
% % % % %     %       plot(N_vec, mu, 'c');
% % % % %     %       plot(N_vec, sigma, 'k');
% % % % %     plot(N_vec, mat2vec(same_inds_prob{main_ind,j}(1,:)), 'm.'); % MZ had the same gaussian
% % % % %     plot(N_vec, mat2vec(same_inds_prob{main_ind,j}(2,:)), 'c.'); % DZ had the same gaussian
% % % % % end % loop on shared environment variable
% % % % % legend('r_{MZ}^2', 'r_{DZ}^2', 'h_{pop}^2', 'h_{loci}^2', 'h_{loci2}^2', ...
% % % % %     'P_{DZ}', 'P_{MZ}'); % '\mu', '\sigma',
% % % % % %%title(['Familial correlations for ' operation ' architecture. Single Gaussian heritability h_x^2=' num2str(h_x)]);
% % % % % xlabel('N'); ylabel('Heritabilities');
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_figure; hold on; % Plot just the r_DZ vs. r_MZ for different models
%for k=1:num_h % plot for different values of the heritabity
for j=1:num_shared_env
    plot(100*reshape(r_mz_vec(:,j,:), num_h, num_N)', ...
        100*reshape(r_dz_vec(:,j,:), num_h, num_N)', ...
        symbol_vec{1+0*j},  'linewidth', 1);
    plot(100*reshape(r_mz_vec(:,j,2:end-3), num_h, num_N-4)', ...
        100*reshape(r_dz_vec(:,j,2:end-3), num_h, num_N-4)', ...
        symbol_vec{(3-j)*5},  'linewidth', 1, 'markersize', 4);
    
    plot_colors = get(gca, 'colororder');
    for k=1:num_h
        plot(100*r_mz_vec(k,j,1), 100*r_dz_vec(k,j,1), ...
            symbol_vec{(3-j)*5}, 'color', plot_colors(k,:),  'linewidth', 1);
    end
    
    %     plot(100*reshape(r_mz_vec(:,j,:), num_h, num_N)', ...
    %         100*reshape(r_dz_vec(:,j,:), num_h, num_N)', ...
    %         '.', 'linewidth', 1);
    
    
    
    %     plot(100*reshape(r_mz_vec(:,j,:), num_h, num_N)', ...
    %         100*reshape(r_dz_vec(:,j,:), num_h, num_N)', ...
    %         'o', 'linewidth', 1);
    for i=1:num_h
        %        plot(100*r_mz_vec(i,j,end), 100*r_dz_vec(i,j,end), 'm*');
        hh = arrow([100*r_mz_vec(i,j,end-1), 100*r_dz_vec(i,j,end-1)], ...
            [100*r_mz_vec(i,j,end), 100*r_dz_vec(i,j,end)], ...
            'length', 30); % Need to color this!!!, ...
        %            'EdgeColor','r');
        set(hh, 'FaceColor', plot_colors(i,:));
        set(hh, 'EdgeColor', plot_colors(i,:));
    end
    
    
    % plot(100*reshape(r_mz_vec(main_ind,:,:), num_shared_env, num_N)', ...
    %     2*100*reshape(r_dz_vec(main_ind,:,:), num_shared_env, num_N)', ...
    %     'ko', 'linewidth', 2);
    % plot(100*reshape(r_mz_vec(main_ind,:,:), num_shared_env, num_N)', ...
    %     2*100*reshape(r_dz_vec(main_ind,:,:), num_shared_env, num_N)', ...
    %     'k*', 'linewidth', 2);
    
end % loop on shared environment
% % % No need for text in the fig
% % for i=[1:length(N_vec)-3 length(N_vec)] % plot also numbers
% %     for k=num_h:num_h % num_shared_env % plot just for one shared environment
% %         text(2+100*r_mz_vec(k,num_shared_env,i) + 0.6*(sign(1.5-i)+1), ...
% %             max(1, 100*r_dz_vec(k,num_shared_env,i)), ...
% %             num2str(N_vec(i)), 'fontsize', 10, 'fontweight', 'bold');
% %     end
% % end
plot(0:100, 0.5.*(0:100), 'k--');

% Plot empirical data from Visscher
visscher_table_file = '../../common_disease_model/data/epidemiology/visscher_QTL_twin_correlation_statistics.txt';
R = ReadDataFile(visscher_table_file,[],[],1);
R = struct_by_inds(R, 2:87); % take 86 traits
R.all_r_DZ = cell2mat(R.all_r_DZ);
R.all_r_MZ = cell2mat(R.all_r_MZ);
plot(100*R.all_r_MZ, 100*R.all_r_DZ, '.k', 'linewidth', 1, 'MarkerSize', 11);
xlabel('r_{MZ} (%)'); ylabel('r_{DZ} (%)'); legend(h_legend_vec,2); % shared_env_legend_vec);

% % % Don't put a title
% % if(isoheritability_flag) % adjust title to parameters
% %     title(['Familial correlations for ' operation ' architecture. Common environment V_c/V_e=' ...
% %         '0% (solid), 50% (dotted), 100% (dashed)']); %num2str(100*h_x_vec(1),2) '%']);
% % else
% %     title(['Familial correlations for ' operation ' architecture. Common environment V_c (1-liab) =' ...
% %         '0% (circle), 50% (star)']); %num2str(100*h_x_vec(1),2) '%']);
% % end

xlim([0 105]);
legend boxoff;

my_saveas(gcf, fig_outfile, format_fig_vec);

%'../../common_disease_model/figs/qtl_r_MZ_vs_r_DZ_' operation ...
%    '_iso_' num2str(isoheritability_flag)], format_fig_vec); % save shared-environment curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
