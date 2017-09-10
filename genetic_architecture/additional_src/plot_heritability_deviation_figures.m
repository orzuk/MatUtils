% Plot heritability estimation for different lambda_R values
function plot_heritability_deviation_figures(heritability_deviation_file)

if(~exist('h_pop_str', 'var') || isempty(h_pop_str))
    h_pop_str = {'ACE', 'ADE', 'MZ', 'DZ', 'PO'}; num_h_pop = length(h_pop_str);
end

initial_plots = 0; % additional optional plots
AssignGeneralConstants;
format_fig_vec = {'epsc'};
load(heritability_deviation_file);

h_x_shared_env_vec = h_shared_env_vec;  h_shared_env_vec = []; 

if(~exist('h_phantom_vec', 'var'))
    h_phantom_vec = heritability_overestimation_vec;
end
if(~exist('h_liab_MZ_vec', 'var'))
    h_liab_MZ_vec = zeros(size(lambda_mz_vec)); h_liab_DZ_vec = h_liab_MZ_vec ;
    for i1=1:size(lambda_mz_vec,1)
        for i2=1:size(lambda_mz_vec,2)
            loop_on_i1 = i1, loop_on_i2 = i2
            for i3=1:size(lambda_mz_vec,3)
                for i4=1:size(lambda_mz_vec,4)
                    h_liab_MZ_vec(i1,i2,i3,i4) = ...
                        familial_risk_to_heritability(lambda_mz_vec(i1,i2,i3,i4), 'liability', ...
                        mu_vec(i1), 1);
                    h_liab_DZ_vec(i1,i2,i3,i4) = ...
                        familial_risk_to_heritability(lambda_s_vec(i1,i2,i3,i4), 'liability', ...
                        mu_vec(i1), 0.5);
                end
            end
        end
    end
    h_liab_ADE_vec = 2.*h_liab_DZ_vec-h_liab_MZ_vec; % 4*r_DZ - r_MZ
    
    save('-append', heritability_deviation_file, ...
        'h_liab_MZ_vec', 'h_liab_DZ_vec', 'h_liab_ADE_vec');
end
if(~exist('h_phantom_ADE_vec', 'var'))
    h_phantom_ADE_vec = 1- h_liab_loci_vec./ h_liab_ADE_vec;
    h_phantom_MZ_vec = 1- h_liab_loci_vec./ h_liab_MZ_vec;
    h_phantom_DZ_vec = 1- h_liab_loci_vec./ h_liab_DZ_vec;
    save('-append', heritability_deviation_file, ...
        'h_phantom_MZ_vec', 'h_phantom_DZ_vec', 'h_phantom_ADE_vec');
end
if(~exist('H01_phantom_vec', 'var'))
    H01_phantom_vec = 1 - h01_loci ./ H01;
end

% Test heritabilities relations
should_be_zero = max(abs(h_liab_pop_vec{1}(:) - (2*h_liab_MZ_vec(:)-h_liab_DZ_vec(:))))


%fig_outdir = '../../common_disease_model/docs/figs/new/';
fig_outdir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/';
data_outdir = '../../common_disease_model/docs/pnas/genetic_interactions/tables/';
fig_outfile = fullfile(fig_outdir, 'fig');  % New! set fig numbering like in the paper !


% h_phantom_vec, ...
%     lambda_s_overestimation_vec, good_inds, lambda_mz_vec, lambda_s_vec, mu_multidim_vec, ...
%     mu_vec, h_x_one_liab_vec, disease_mu_vec, main_mu_ind, max_N, ...
%     MAX_LAMBDA_MZ, MIN_LAMBDA_MZ) % set independent parameters)
lambda_mz_vec(lambda_mz_vec > 98) = max(lambda_mz_vec(lambda_mz_vec > 98), 100); % round lambdas (fix numerical errors)

num_shared_env = length(h_x_shared_env_vec);
num_h = length(h_x_one_liab_vec);
num_mu = length(mu_vec);
if(~exist('N_vec', 'var'))
    N_vec = 1:max_N;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New figure: plot lambda's as function of family relation
figure; hold on; legend_vec = [];
kinship_vec = 2.^(0:-1:-5); h_phantom_vec(:,:,:,1)=0; 
for i=1:6
    plot(kinship_vec, log2(lambda_R{2,5,1,i}), color_vec(i));
    legend_vec{i} = ['N=' num2str(i) ', \pi_{phantom} = ' num2str(h_phantom_vec(2,5,1,i)*100,3) '%'];
end
for i=1:6
    plot(kinship_vec, log2(lambda_R{2,5,1,i}), [color_vec(i) '*']);
end

xlabel('Kinship k_R'); ylabel('$\log(\lambda_R)$', 'interpreter', 'latex');
title(['Familial risk for $LP(N, \mu= 1\%, h_{path}^2=90\%)$'], 'interpreter', 'latex');  
legend(legend_vec,2); 
my_saveas(gcf, fullfile(fig_outdir, 'familial_risk', 'familial_risk_LP_model'), {'epsc', 'pdf', 'jpg', 'fig'}); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[wray_disease_vec, wray_disease_short_vec, ... % Plot real disease data from Wray's paper
    wray_lambda_s_vec, wray_lambda_mz_vec, wray_prevalence_vec, ...
    ~, ~, ~, ~, wray_lambda_s_overestimation_vec] = wray_lambda_table(0);
k=1; alpha=0.1; num_N = max_N;

for num_h_shared = [1 length(h_x_shared_env_vec)] % loop on main and supp. info figs
    if(num_h_shared == 1) % main figure
        add_zoomed = 1;
        run_shared_env_vec_inds = 1;
    else
        add_zoomed = 0; % supplumentary
        run_shared_env_vec_inds = union(find(h_x_shared_env_vec == 0), ...
            find(h_x_shared_env_vec == 0.5)); 
    end
    marker_size = [11 3 18 6]; % dot gets the larger font 
    for i_pop = [1 2 5]; % 1:length(h_pop_str)
        switch i_pop
            case 1 % ACE model. Plot once all plots
                plot_type_vec = 0:3
            otherwise
                plot_type_vec = 2;
        end
        for plot_type=plot_type_vec %0:8 % loop on different plots % 0 - lambda_mz vs. lambda_s, 1 - lambda_mz vs. difference, 2 - lambda_mz vs. h_phantom
            if((plot_type < 2) && add_zoomed)
                cur_add_zoomed = 1;
            else
                cur_add_zoomed = 0;
            end
            for zoom_flag = 0:cur_add_zoomed % First is figure, second is inset
                if(zoom_flag == 0)
                    h_fig = full_figure(0); % plot lambda_s vs. lambda_mz for real disease and for models
                else
                    h_axes = get(h_fig,'CurrentAxes');
                    switch plot_type
                        case 0
                            h_inset = axes('pos',[0.16 0.51 .35 .4]);
                        case 1
                            h_inset = axes('pos',[0.17 0.51 .42 .4]);
                    end
                end
                range_vec = 1:length(mu_vec);
                for i=range_vec % loop on mu
                    for j_shared=1:length(run_shared_env_vec_inds) %  1:num_h_shared % loop on shared environment
                        i_shared = run_shared_env_vec_inds(j_shared);
                        cur_good_inds = intersect(find(lambda_mz_vec(i,:,i_shared,max_N,k) < MAX_LAMBDA_MZ), ...
                            find(lambda_mz_vec(i,:,i_shared,max_N,k) > MIN_LAMBDA_MZ));
                        cur_good_inds = intersect(cur_good_inds, [1:20]); % f 4 7 10 13 16 19]);
                        cur_good_inds = 1:length(lambda_mz_vec(i,:,i_shared,max_N,k)); %TEMP: plot every lambda
                        num_lambda_mz = length(cur_good_inds);
                        
                        %    for j=1:max(1,length(h_x_one_liab_vec)*min(1,max_N-1)) % loop on heritability/lambda_mz
                        switch plot_type
                            case 0 % lambda_mz vs lambda_s
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k),  num_lambda_mz, num_N)';
                                y_vec = reshape(lambda_s_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                loglog(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i},  'linewidth', 1); hold on;
                                loglog(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(j_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared)); hold on; % marker_size(i_shared)
                                %                    '.',  'linewidth', 1, 'markersize', 12); hold on;
                                plot_colors = get(gca, 'colororder');
                                if(num_h_shared == 1)
                                    cur_marker_size = marker_size(j_shared+2);
                                else
                                    cur_marker_size = marker_size(j_shared+2);
                                end

                                for j=1:num_h % plot the first big point
                                    loglog(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), ...
                                        'markersize', marker_size(j_shared+2));
                                end
                                
                                
                                
                                %                     reshape(lambda_mz_vec(i,cur_good_inds,i_shared,1:end-1,k), num_lambda_mz, max_N-1)', ...
                                %                         reshape(lambda_s_vec(i,cur_good_inds,i_shared,1:end-1,k), num_lambda_mz,  max_N-1)', ...
                                
                            case 1 % lambda_mz vs lambda_s over-estimation
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*lambda_s_overestimation_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(i_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                                
                                
                            case 2 % phantom heritability vs. lambda_mz
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
%                                y_vec = reshape(100*h_phantom_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*pi_liab_phantom_vec{i_pop}(i,:,i_shared,:,k), num_lambda_mz, num_N)';

                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, ...
                                    'markersize', marker_size(j_shared)+3); % increase font size
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), ...
                                        'markersize', marker_size(j_shared+2)+3);
                                end
                                
                            case 3 %Here plot h_epi vs. h_phantom or h_loci
                                x_vec = reshape(100*min(1,h_liab_pop_vec{i_pop}(i,:,i_shared,:,k)), num_lambda_mz, num_N)';
                                y_vec = reshape(100*h_liab_loci_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                plot(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                plot(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    plot(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                                
                            case 4 % phantom heritability (like 2) with ADE
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*h_phantom_ADE_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                            case 5 % phantom heritability (like 2) with MZ
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*h_phantom_MZ_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                            case 6 % phantom heritability (like 2) with DZ
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*h_phantom_DZ_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(i_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(i_shared+2));
                                end
                            case 7 % phantom heritability (like 2) with binary h01
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*h01_phantom_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                            case 8 % phantom heritability (like 2) with binary h01. Use true H01 as epi
                                x_vec = reshape(lambda_mz_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                y_vec = reshape(100*H01_phantom_vec(i,:,i_shared,:,k), num_lambda_mz, num_N)';
                                semilogx(x_vec(1:end,:), y_vec(1:end,:), ...
                                    symbol_vec{i}, 'linewidth', 1); hold on;
                                semilogx(x_vec(2:end-3,:), y_vec(2:end-3,:), ... % 12-7*(i_shared-1)
                                    symbol_vec{12-5*(j_shared-1)}, 'linewidth', 1, 'markersize', marker_size(j_shared));
                                hold on;
                                for j=1:num_h
                                    semilogx(x_vec(1,j), y_vec(1,j), ... % 12-7*(i_shared-1)
                                        symbol_vec{12-5*(j_shared-1)}, 'color', plot_colors(j,:), 'markersize', marker_size(j_shared+2));
                                end
                                
                        end % switch plot type
                        if(zoom_flag == 0)
                            arrow_vec = 1:num_lambda_mz;
                        else % draw arrows only in zoomed area
                            arrow_vec = find(x_vec(end-1,:) < 5);
                        end
                        if(zoom_flag == 0)
                            arrow_len = 25;
                        else % make arrowhead smaller (since it's zoomed
                            arrow_len = 7; %                        set(hh, 'length', 5); % make smaller arrowhead
                        end
                        
                        for j=arrow_vec
                            hh = arrow([alpha*x_vec(end-1,j)+(1-alpha)*x_vec(end,j), ...
                                alpha*y_vec(end-1,j)+(1-alpha)*y_vec(end,j)], ...
                                [x_vec(end,j), y_vec(end,j)], 'length', arrow_len); % Need to color this!!!, ...
                            if(i_shared==1)
                                set(hh, 'FaceColor', plot_colors(j,:));
                            else % shared environment
                                set(hh, 'FaceColor', 'w');
                            end
                            set(hh, 'EdgeColor', plot_colors(j,:));
                            set(hh, 'lineStyle', symbol_vec{i});
                            if(i > 1) % arrow not solid line. Draw another one at the tip
                                alpha2=0.001;
                                hh = arrow([alpha*x_vec(end-1,j)+(1-alpha)*x_vec(end,j), ...
                                    alpha*y_vec(end-1,j)+(1-alpha)*y_vec(end,j)], ...
                                    [x_vec(end,j), y_vec(end,j)], 'length', arrow_len); % Need to color this!!!, ...
                                if(i_shared==1)
                                    set(hh, 'FaceColor', plot_colors(j,:));
                                else % shared environment
                                    set(hh, 'FaceColor', 'w');
                                end
                                set(hh, 'EdgeColor', plot_colors(j,:));
                            end
                            
                        end
                        if((i_shared==1) && (mu_vec(i) == 0.01)) % text exampled
                            plot(x_vec(3,3), y_vec(3,3), '*r', 'markersize', 12);
                        end
                        
                    end % loop on shared environment
                    %    end % loop on heritaility
                end % loop on mu
                
                
                xlabel('$\lambda_{MZ}$', 'Interpreter','latex');  % 1 - lambda_s/lambda_mz
                legend_place = [4 1 4 2 4 4 4 4 4];
                switch  plot_type
                    case 0 % lambda_mz vs. lambda_s 
                        if(num_h_shared == 1)
                            title_str = '(b.)'; % 'supp_fig2b'; % ['\lambda_{mz} vs. \lambda_s for'];
                            if(zoom_flag == 0)
                                set(gca,'Box','off'); % remove frame 
                            end
                        else % more than one shared environment
                            title_str = 'supp_fig4a';
                        end
                        ylabel('$\lambda_s$', 'Interpreter','latex');
                        loglog(wray_lambda_mz_vec, wray_lambda_s_vec , '.k', 'markersize', 20);
                        fig_outfile = fullfile(fig_outdir, title_str); % '_underestimation'];
                    case 1
                        if(num_h_shared == 1)
                            title_str = 'supp_fig3a'; % ['\lambda_{mz} vs. \lambda_s under-estimation for'];
                        else
                            title_str = 'not_used_supp_fig';
                        end
                        
                        fig_outfile = fullfile(fig_outdir, title_str); % '_underestimation'];
                        ylabel('$\lambda_s$ overestimation: $\frac{\lambda_{s,obs}}{\lambda_{s,exp}} - 1$', 'interpreter', 'latex');
                        semilogx(wray_lambda_mz_vec, 100*wray_lambda_s_overestimation_vec , '.k', 'markersize', 20);
                    case 2 % lambda_z vs. phantom (no real data)
                        if(num_h_shared == 1)
                            switch h_pop_str{i_pop}
                                case 'ACE'
                                    title_str = 'fig1b'; % ['\lambda_{mz} vs. phantom heritability for'];
                                otherwise
                                    title_str = ['not_used_fig2a_' h_pop_str{i_pop}];
                            end
                        else % multiple shared environments
                            switch h_pop_str{i_pop} % this is with shared environment
                                case 'ACE'
                                    title_str = 'supp_fig4b';
                                otherwise
                                    title_str = ['not_used_supp_fig4b_' h_pop_str{i_pop}];
                            end
                        end
                        
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('$\pi_{phantom} (\%)$', 'interpreter', 'latex');
                        ylim([0 100]); % max(100*h_phantom_vec(:))]); % set y limit
                        set(gca,'box','off')
%                        xlim([0.9 600]);
                    case 3 %Here plot h_epi vs. h_phantom or h_loci
                        if(num_h_shared == 1)
                            title_str = 'supp_fig3b'; % ['h_{pop}^2 vs. h^2 for'];
                        else
                            title_str = 'not_used_supp_fig';
                        end
                        
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('$h_{all}^2 (\%)$', 'interpreter', 'latex');
                        ylim([0 max(100*h_liab_loci_vec(:))]); % set y limit
                        xlim([0 101]);
                        plot(0:100, 0:100, 'k--', 'linewidth', 1);
                        xlabel('$h_{pop}^2 (\%)$', 'interpreter', 'latex');  % 1 - lambda_s/lambda_mz
                        
                    case 4 % phantom ADE (no real data)
                        if(num_h_shared == 1)
                            title_str = 'not_used_supp_fig_ADE'; % ['\lambda_{mz} vs. phantom heritability for'];
                        else
                            title_str = 'supp_fig4ADE';
                        end
                        
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('phantom heritability ADE (\%) : $1- h^2 / h_{pop}^2$', 'interpreter', 'latex');
                        ylim([0 max(100*h_phantom_ADE_vec(:))]); % set y limit
                    case 5 % phantom MZ (no real data)
                        if(num_h_shared == 1)
                            title_str = 'not_used_supp_fig_MZ'; % ['\lambda_{mz} vs. phantom heritability for'];
                        else
                            title_str = 'supp_fig4MZ';
                        end
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('phantom heritability MZ (\%) : $1- h^2 / h_{pop}^2$', 'interpreter', 'latex');
                        ylim([0 max(100*h_phantom_MZ_vec(:))]); % set y limit
                    case 6 % phantom DZ (no real data)
                        if(num_h_shared == 1)
                            title_str = 'not_used_supp_fig_DZ'; % ['\lambda_{mz} vs. phantom heritability for'];
                        else
                            title_str = 'supp_fig4DZ';
                        end
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('phantom heritability DZ (\%) : $1- h^2 / h_{pop}^2$', 'interpreter', 'latex');
                        ylim([0 max(100*h_phantom_DZ_vec(:))]); % set y limit
                    case 7 % phantom 01 (no real data)
                        if(num_h_shared == 1)
                            title_str = 'not_used_supp_fig_h01'; % ['\lambda_{mz} vs. phantom heritability for'];
                        else
                            title_str = 'supp_fig4h01';
                        end
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('phantom heritability h_{01} (\%) : $1- h_{01}^2 / h_{01,pop}^2$', 'interpreter', 'latex');
                        ylim([0 max(100*h01_phantom_vec(:))]); % set y limit
                    case 8 % phantom 01 (no real data)
                        if(num_h_shared == 1)
                            title_str = 'not_used_supp_fig_H01_MZ'; % ['\lambda_{mz} vs. phantom heritability for'];
                        else
                            title_str = 'supp_fig4H01_MZ';
                        end
                        fig_outfile = fullfile(fig_outdir, title_str); % '_phantom'];
                        ylabel('phantom heritability h_{01} (\%) : $1- h_{01}^2 / H_{01,pop}^2$', 'interpreter', 'latex');
                        ylim([0 max(100*H01_phantom_vec(:))]); % set y limit
                        
                end % switch plot_type
                
                if(any(plot_type == [0 1 2 4 5 6 7 8])) % set x ticks
                    x_ticks = get(gca, 'xtick');
                    new_x_ticks = []; % x_ticks;
                    for j=1:length(x_ticks)-1
                        new_x_ticks = [new_x_ticks [1 2 5].*x_ticks(j)];
                    end
                    new_x_ticks = [new_x_ticks x_ticks(end)];
                    set(gca, 'xtick', new_x_ticks);
                    set(gca, 'xticklabel', num2str(new_x_ticks'))
                end
                if(any(plot_type == [0])) % set x ticks
                    y_ticks = get(gca, 'ytick');
                    new_y_ticks = []; % y_ticks;
                    for j=1:length(y_ticks)-1
                        new_y_ticks = [new_y_ticks [1 2 3 4 5].*y_ticks(j)];
                    end
                    new_y_ticks = [new_y_ticks y_ticks(end)];
                    set(gca, 'ytick', new_y_ticks);
                    set(gca, 'yticklabel', num2str(new_y_ticks'))
                end
                
                
                
                %    title_str = [title_str ' Common environment V_c/V_e=0% (circle), 50% (star)']; % No shared environment!!!
                if(isoheritability_flag)
                    lambda_legend_vec = num2str_cell(num2cell(input_lambda_mz_vec));
                    for i=1:length(lambda_legend_vec)
                        lambda_legend_vec{i} = ['\lambda_{MZ}=' lambda_legend_vec{i}];
                    end
                else % here legend is h_x
                    lambda_legend_vec = num2str_cell(num2cell(100*h_x_one_liab_vec), 2, [], 1);
                    for i=1:length(lambda_legend_vec)
                        lambda_legend_vec{i} = ['h_{pathway}^2=' lambda_legend_vec{i}];
                    end
                end
                if((zoom_flag==0) & (plot_type ~= 2))
                    leg_hand = legend(lambda_legend_vec(cur_good_inds),legend_place(plot_type+1));
                    %        leg_hand = legend(lambda_legend_vec(cur_good_inds),'location', 'best');
                    leg_pos = get(leg_hand, 'position');
                    delta_pos = 0.1;
                    switch legend_place(plot_type+1)
                        case 1
                            leg_pos(1:2) = leg_pos(1:2) + delta_pos.*[-1 -1];
                        case 2
                            leg_pos(1:2) = leg_pos(1:2) +  delta_pos.*[1 -1];
                        case 3
                            leg_pos(1:2) = leg_pos(1:2) +  delta_pos.*[1 1];
                        case 4
                            leg_pos(1:2) = leg_pos(1:2) +  delta_pos.*[-1 1];
                    end
                    if(plot_type == 1)
                        legend('boxoff');
                        set(leg_hand, 'position', leg_pos);
                    else
                        set(leg_hand,'box','on', ...
                            'position', leg_pos, ...
                            'color','w', ...
                            'ycolor', [0.8 0.8 0.8], ...
                            'xcolor', [0.8 0.8 0.8], ...
                            'visible','on'); % this should show as white in eps !!!
                    end
                end
                %    set(leg_hand, 'edgecolor', 'w');
                %    legend(leg_hand, 'boxoff');
                
                % % %     for i=1:length(mu_vec)
                % % %         title_str = [title_str ' \mu=' num2str(mu_vec(i)) ' (' symbol_str_vec{i} ')'];
                % % %     end
                title([repmat(' ', 1, 100) str2title(title_str)], 'fontsize', 16, 'fontweight', 'bold');
                
                %     \mu=' num2str(mu_vec(1)) ' (solid), \mu=' ...
                %         num2str(mu_vec(2)) ' (dotted)' num2str(mu_vec(3)) ' (dashed)']);
                if(zoom_flag) %  add_zoomed)
                    switch  plot_type
                        case {0,1}
                            xlim([1 5]); % zoom-in
                            legend('hide');
                            title(''); xlabel(''); ylabel([]);
                            switch plot_type
                                case 0
                                    ylim([1 3.3]);
                            end
                            %                title([repmat(' ', 1, 100) str2title(title_str) '-inset'], 'fontsize', 12);
                            %                my_saveas(gcf, [fig_outfile '_inset'], format_fig_vec)
                    end
                    
                end
            end % loop on zoomed
            my_saveas(gcf, fig_outfile, format_fig_vec)
        end  %loop on plot type
    end % loop on i_pop
end % loop on main and supp info figs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(initial_plots)
    full_figure; %hold on; % plot two types of deviations
    plot(mat2vec(h_phantom_vec(good_inds)), ...
        mat2vec(-lambda_s_overestimation_vec(good_inds)), 'k*');
    %    plot(heritability_overestimation_MZ_vec(:), -lambda_s_overestimation_vec(:), 'r.');
    for i=1:length(mu_vec) % 0.04 % [0.001 0.01 0.05 0.1]
        for j=1:length(h_x_one_liab_vec)
            for N=1:max_N %
                for k=[1]
                    if( (lambda_mz_vec(i,j,N,k) < MAX_LAMBDA_MZ) && (lambda_mz_vec(i,j,N,k) > MIN_LAMBDA_MZ) )
                        text(h_phantom_vec(i,j,N,k), -lambda_s_overestimation_vec(i,j,N,k), ...
                            ['MLT(' num2str(N) ',' num2str(k) ',' ...
                            num2str(disease_mu_vec(i,j,N,k),2) ',' num2str(lambda_mz_vec(i,j,N,k),3) ')'], ...
                            'fontsize', 8);
                    end
                end
            end
        end
    end
    xlabel('h^2 overestimation: (h_{pop}^2 - h_{loci}^2) / h_{pop}^2');
    ylabel('\lambda_s underestimation: (\lambda_s^{LT} - \lambda_s^{obs})/\lambda_s^{LT}');
    xlim([0 1]); ylim([0 0.2]); % max(-lambda_s_overestimation_vec(:))*1.05]);
    my_saveas(gcf, '../../common_disease_model/figs/lambda_s_and_heritability_deviations2', format_fig_vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for plot_type = 0:1
        for main_fig_ind = 0:1 % plot main and supp. info figs
            if(main_fig_ind)
                range_vec = main_mu_ind;
            else
                range_vec = 1:length(mu_vec);
            end
            
            full_figure; % plot two types of deviations with a curve while varying heritability
            for i=range_vec
                %        for j=1:length(h_x_one_liab_vec)
                for N=1:max_N %
                    for k=[1]
                        cur_good_inds = intersect(find(lambda_mz_vec(i,:,N,k) < MAX_LAMBDA_MZ), ...
                            find(lambda_mz_vec(i,:,N,k) > MIN_LAMBDA_MZ));
                        
                        cur_good_inds = intersect(cur_good_inds, [1 3 5 7 9 11 13 15 17 19]);
                        if(plot_type)
                            if(main_fig_ind)
                                plot(100*h_phantom_vec(i,cur_good_inds,N,k), ...
                                    -100*lambda_s_overestimation_vec(i,cur_good_inds,N,k), color_vec(N), ...
                                    'linewidth', 2);
                            else
                                h_leg(i,N,k) = plot(100*h_phantom_vec(i,cur_good_inds,N,k), ...
                                    -100*lambda_s_overestimation_vec(i,cur_good_inds,N,k), ...
                                    [color_vec(N) symbol_vec{i+4}],  'linewidth', 2);
                            end
                        else % Here plot h_epi vs. overestimation
                            if(exist('h_liab_pop_vec', 'var'))
                                plot(100*h_liab_pop_vec(i,cur_good_inds,N,k), ...
                                    -100*h_phantom_vec(i,cur_good_inds,N,k), color_vec(N), ...
                                    'linewidth', 2); % flip over-estimation sign
                            end
                        end
                        cur_text_inds = find(-lambda_s_overestimation_vec(i,cur_good_inds,N,k) < 0.2, 1, 'first');
                    end
                end
                %        end
            end
            for i=range_vec % loop on mu (only for supp. info.)
                for N=1:max_N % loop on N
                    for j=1:max(1,length(h_x_one_liab_vec)*min(1,N-1)) % loop on heritability/lambda_mz
                        for k=[1]
                            %                     if(~main_fig_ind)
                            %                        mu_multidim_vec(i,j,N,k) = mu_vec(i);
                            %                     end
                            if(N == 1)
                                range_str = [num2str(lambda_mz_vec(i,j,N,k),3) '-' ...
                                    num2str(round(min(1./mu_vec(main_mu_ind), max(lambda_mz_vec(good_inds)))))];
                            else
                                range_str = num2str(lambda_mz_vec(i,j,N,k),3);
                            end
                            if(main_fig_ind)
                                if( (lambda_mz_vec(i,j,N,k) < MAX_LAMBDA_MZ) && (lambda_mz_vec(i,j,N,k) > MIN_LAMBDA_MZ) )
                                    text(100*h_phantom_vec(i,j,N,k)+1, ...
                                        -100*lambda_s_overestimation_vec(i,j,N,k), ...
                                        range_str, ...
                                        'fontsize', 12, 'color', color_vec(N), 'fontweight', 'bold');
                                end
                            end
                        end
                    end
                end
            end
            if(main_fig_ind)
                main_mu_good_inds = find(mu_multidim_vec(:) == mu_vec(main_mu_ind));
                good_inds = intersect(good_inds, main_mu_good_inds);
                plot(100*h_phantom_vec(good_inds), -100*lambda_s_overestimation_vec(good_inds), 'k*');
                plot(100*h_phantom_vec(1,1,1,1)+0.01, -100*lambda_s_overestimation_vec(1,1,1,1)+0.01, 'k*'); % This is for zero
            end
            legend_vec = cell(1,max_N);
            for i=1:max_N
                legend_vec{i} = ['N=' num2str(i)];
            end
            if(main_fig_ind || (~plot_type))
                legend(legend_vec, 2);
            else
                legend(h_leg(1,:,1), legend_vec, 2);
            end
            if(plot_type)
                xlabel('h^2 overestimation (%): (h_{pop}^2 - h_{loci}^2) / h_{pop}^2');
                ylabel('\lambda_s underestimation  (%): (\lambda_s^{LT} - \lambda_s^{obs})/\lambda_s^{LT}');
                xlim([0 80]); ylim([0 40]); % (percents) get rid of too large underestimations % max(-lambda_s_overestimation_vec(:))*1.05]);
                fig_outfile = '../../common_disease_model/figs/lambda_s_and_heritability_deviations_curve_lambda_mz';
            else
                xlabel('$h_{pop}^2$', 'Interpreter','latex');
                ylabel('phantom heritability (\%): $(h_{pop}^2 - h_{loci}^2) / h_{pop}^2$', ...
                    'interpreter', 'latex');
                %    xlim([0 80]); ylim([0 40]); % (percents) get rid of too large underestimations % max(-lambda_s_overestimation_vec(:))*1.05]);
                fig_outfile = '../../common_disease_model/figs/h_pop_vs_phantom_heritability';
                
            end
            if(main_fig_ind)
                fig_outfile = [fig_outfile '_main'];
            end
            my_saveas(gcf, fig_outfile, format_fig_vec)
        end  % loop on main/supp. info fig
    end % loop on plot type
    
end % initial plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save stuff in latex !!! 
h_shared_env_vec = zeros(num_mu, num_h, num_shared_env, num_N); 

for i=1:num_mu % loop on mu
    for j_h = 1:num_h % loop on heritability 
        for j_shared=1:num_shared_env %  1:num_h_shared % loop on shared environment
            h_shared_env_vec(i,j_h,j_shared,:) = h_liab_MZ_vec(i,j_h,j_shared,:) - h_liab_MZ_vec(i,j_h,1,:);
        end
    end
end

num_vals = num_h*num_shared_env*num_N*num_mu;
h_phantom_vec(h_phantom_vec<0.001) = 0; % rounding down !!!!!
S = [vec2column(mat2vec(repmat(N_vec, num_h*num_shared_env*num_mu, 1))) ...
    vec2column(repmat(mat2vec(repmat(100*h_x_one_liab_vec, num_mu, 1)),  num_shared_env*num_N, 1)) ...
    vec2column(repmat(mat2vec(repmat(100*h_x_shared_env_vec, num_h*num_mu , 1)), num_N, 1)) ... % c_R
    vec2column( repmat(mat2vec(repmat(1-h_x_one_liab_vec, num_mu, 1)),  num_shared_env*num_N, 1) .* ...
    repmat(mat2vec(repmat(100*h_x_shared_env_vec, num_h*num_mu , 1)), num_N, 1) ) ... % h_{c,env}
    vec2column(repmat(100*mu_vec, 1, num_N*num_shared_env*num_h)) ...
    mat2vec(lambda_mz_vec) mat2vec(lambda_s_vec) ...
    mat2vec(100*h_liab_loci_vec) mat2vec(100*h_shared_env_vec)]; % new! add shared environment vec 

S_header = {'$k$', '$h^2_{path} (\%)$', '$c_R (\%)$', '$V_{path,c} (\%)$', ...
    '$\mu (\%)$', '$\lambda_{MZ}$', '$\lambda_{s}$', ...
    '$h_{all}^2 (\%)$', '$V_c (\%)$'}; % New: add shared environment for the trait
for i_pop = [1 2 5] % 1:length(h_pop_str)
    S = [S mat2vec(100*min(1,h_liab_pop_vec{i_pop})) ...
        mat2vec(100*max(0,pi_liab_phantom_vec{i_pop}))];   % save figures data in latex format
    S_header = [S_header ['$h^{2 (' h_pop_str{i_pop} ')}_{pop} (\%)$'] ...
        ['$\pi_{phan}^{(' h_pop_str{i_pop} ')} (\%)$']];
end
S = [S_header' num2cell(S)']';

for i=2:size(S,1)
    for j=4:size(S,2)
        S{i,j} = sprintf('%.1f', S{i,j});
    end
end

savecellfile(S, fullfile(data_outdir, 'fig2_data_disease_traits_tab.txt'),[],1);  % Save figures data
S = S(:, [1:3 5:end]); % get rid of shared environment within the pathway 
S_latex = latex(S, 2, precision);
S_latex = mat2cell(S_latex, ones(size(S_latex,1),1));
S_latex{1} = strrep(S_latex{1},  '|c', '|r'); % align to right
%for i=2+num_N*num_mu:num_N*num_mu:length(S_latex)-num_N*num_mu % add lines when k is changed
for i=2+(num_h*num_shared_env*num_mu):(num_h*num_shared_env*num_mu):length(S_latex)-(num_h*num_shared_env*num_mu)
    S_latex{i} = [S_latex{i} '  \hdashline'];
end


tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
    '\begin{center}', '{\tiny', ... % '{\scriptsize', ...
    ['\caption[Model parameters for disease traits]{Model parameters for ' ...
    ' the LP model $\archd(k,\hx,c_R,\mu)$ for disease traits. \label{table:fig2_data_disease_traits}}'], ...
    S_latex{1}};
tab_footer = {'}',  ...
    '\end{center}', '\end{table}', '', '\newpage', ''};
%S_latex = [tab_header  S_latex(2:end)' tab_footer]';
%S_latex = strrep_cell(S_latex, 'tabular', 'longtable');
S_latex = split_latex_table(S_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 60, -1); % 50 for scriptsize, 60 for tiny font 
savecellfile(S_latex, fullfile(data_outdir, 'fig2_data_disease_traits.txt'),[],1);  % Save figures data


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_figure; % plot two types of deviations with a curve while varying mu
%    for i=1:length(mu_vec) % 0.04 % [0.001 0.01 0.05 0.1]
h_leg = zeros(length(h_x_one_liab_vec), max_N, 2);
for j=1:length(h_x_one_liab_vec)
    for N=1:max_N %
        for k=[1]
            cur_good_inds = intersect(find(lambda_mz_vec(:,j,N,k) < MAX_LAMBDA_MZ), ...
                find(lambda_mz_vec(:,j,N,k) > MIN_LAMBDA_MZ));
            
            cur_good_inds = intersect(cur_good_inds, [1 3 5 7 9 11 13 15 17 19]);
            if(isempty(cur_good_inds))
                cur_good_inds = 1; % just temp
                %                continue;
            end
            h_leg(j,N,k) = plot(h_phantom_vec(cur_good_inds,j,N,k), ...
                -lambda_s_overestimation_vec(cur_good_inds,j,N,k), color_vec(N));
            text(h_phantom_vec(cur_good_inds(end),j,N,k), ...
                -lambda_s_overestimation_vec(cur_good_inds(end),j,N,k), ...
                ['MLT(' num2str(N) ',' num2str(k) ',\mu' ...
                ',' num2str(lambda_mz_vec(1,j,N,k),2) ')'], ...
                'fontsize', 8); % one text for each curve
        end
    end
end
%    end
plot(h_phantom_vec(good_inds), -lambda_s_overestimation_vec(good_inds), 'k.');
legend(h_leg(1,:,1), 'N=1', 'N=2', 'N=3', 'N=4', 'N=5');
xlabel('h^2 overestimation: (h_{pop}^2 - h_{loci}^2) / h_{pop}^2');
ylabel('\lambda_s underestimation: (\lambda_s^{LT} - \lambda_s^{obs})/\lambda_s^{LT}');
xlim([0 1]); ylim([0 0.2]); % max(-lambda_s_overestimation_vec(:))*1.05]);
my_saveas(gcf, '../../common_disease_model/figs/lambda_s_and_heritability_deviations_curve_mu', format_fig_vec)
