% Plot how much of genetic variation resides at each allele frequency
% according to pop. gen. theory
%
% Input:
% N - effective population size
% s_vec - selection coefficent (abs. value)
%
function plot_var_explained_by_allele_freq(N, s_vec)

AssignGeneralConstants();
[machine machine_delim html_outdir] = get_machine_type();
spectrum_figs_dir = 'c:/research/common_disease_model/docs/allele_freq_spectrum/figs';
% fullfile(html_outdir, '/genetic_architecture/allele_freq_spectrum/figs');
res = 1/N; x_vec = res:res:1-res; % res = 0.001;
N2 = floor(length(x_vec)/2);


allele_freq_vec = zeros(length(s_vec),length(x_vec)); allele_freq_vec2 = zeros(length(s_vec),N2);
cum_allele_freq_vec = allele_freq_vec; cum_allele_freq_vec2 = allele_freq_vec2;
var_freq_vec = allele_freq_vec; var_freq_vec2 = allele_freq_vec2;
cum_var_vec = allele_freq_vec; cum_var_vec2 = allele_freq_vec2;
for i=1:length(s_vec)
    s = s_vec(i);

    % We use 4 (not 2)!. Population size is in individual (not chromosomes)
    allele_freq_vec(i,:) = (exp(-4.*N.*s.*(1-x_vec)) - 1) ./ (x_vec.*(1-x_vec).*(exp(-4.*N.*s) - 1)); % distribution of allele frequencies
    allele_freq_vec(i,:) = normalize_hist(x_vec, allele_freq_vec(i,:));
    allele_freq_vec2(i,:) = allele_freq_vec(i,1:N2) + allele_freq_vec(i,end:-1:end-N2+1); % collapse at zero
    allele_freq_vec2(i,:) = normalize_hist(x_vec(1:N2), allele_freq_vec2(i,:));
    cum_allele_freq_vec(i,:) = cumsum(allele_freq_vec(i,:)); cum_allele_freq_vec(i,:) = cum_allele_freq_vec(i,:) ./ cum_allele_freq_vec(i,end);
    cum_allele_freq_vec2(i,:) = cum_allele_freq_vec(i,1:N2) + 1 - cum_allele_freq_vec(i,end:-1:end-N2+1);
    cum_allele_freq_vec2(i,:) = cum_allele_freq_vec2(i,:) ./ cum_allele_freq_vec2(i,end);
    var_freq_vec(i,:) = allele_freq_vec(i,:) .* x_vec .* (1-x_vec); var_freq_vec(i,:) = normalize_hist(x_vec, var_freq_vec(i,:));
    var_freq_vec2(i,:) = var_freq_vec(i,1:N2) + var_freq_vec(i,end:-1:end-N2+1); var_freq_vec2(i,:) = normalize_hist(x_vec(1:N2), var_freq_vec2(i,:));
    
%    cum_var_vec(i,:) = x_vec - (exp(-4.*N*s.*(1-x_vec)) - exp(-4.*N.*s)) ./ (4.*N.*s); % compute analytically. Wrong formula!!!  
%     cum_var_vec(i,:) = (exp(-4.*N*s) .* (exp(4.*N*s.*x_vec) - 1) - 0.*4.*N.*s.*x_vec) ./ ...
%         ((4.*N.*s) .* (exp(-4.*N.*s) - 1)) - x_vec ./ (exp(-4.*N.*s) - 1); % compute analytically (new formula - numerical problems!!!)
    cum_var_vec(i,:) = cumsum(var_freq_vec(i,:)); % compute empirically
    cum_var_vec(i,:) = cum_var_vec(i,:) ./ cum_var_vec(i,end);
    %    cum_var_vec2 = cumsum(var_freq_vec); cum_var_vec2 = cum_var_vec2 ./ cum_var_vec2(end);
    %    figure; hold on; plot(cum_var_vec, cum_var_vec2, '.'); title('2 ways computing cumulative');
    %    plot(0:res:1, 0:res:1, 'r');
    
    cum_var_vec2(i,:) = cum_var_vec(i,1:N2) + (1 - cum_var_vec(i,end:-1:end-N2+1)); cum_var_vec2(i,:) = cum_var_vec2(i,:) ./ cum_var_vec2(i,end);
    %(2.*z_vec - (exp(-2.*N*s.*(1-z_vec)) - exp(-2.*N.*s.*z_vec)) ./ (2.*N.*s)) ./ f_vec;
end % loop on s values

if(N == 10000)
    param_str = ' N=10^4'; %  num2str(N)];
else
    param_str = num2str(N);
end
legend_vec = num2str_cell(num2cell(s_vec));
for i=1:length(s_vec)
    legend_vec{i} = ['s=' legend_vec{i} ', S=' num2str(4*N*s_vec(i))];
end
z_vec = x_vec;
for two_side_flag = [0 1] % 0 - use dervied allele frequency [0,1]. 1 - use minor allele frequenc [0,0.5]
    switch two_side_flag
        case 0
            x_label = 'DAF'; sup_title_str = 'Derived';
        case 1
            x_label = 'MAF'; sup_title_str = 'Minor';
    end
    for scale = {'linear', 'log'}
        full_figure;
        switch scale{1}
            case 'linear'
                suptitle([sup_title_str ' allele, linear scale ' param_str]);
                if(two_side_flag)
                    x_vec = z_vec(1:N2);
                else
                    x_vec = z_vec;
                end
            case 'log'
                suptitle([sup_title_str ' allele, log scale ' param_str]);
                if(two_side_flag)
                    x_vec = log10(z_vec(1:N2));
                else
                    x_vec = log10(z_vec);
                end
                
                x_label = [x_label ' (log_{10})'];
        end
        for plot_ind = 1:4 % four different plots 
            subplot(2,2,plot_ind); hold on;
            for i=1:length(s_vec) % one plot for each selection coefficient
                switch plot_ind
                    case 1 % allele frequency density
                        if(two_side_flag)
                            plot(x_vec, allele_freq_vec2(i,:), color_vec(i), 'linewidth', 2); % Plot allele freq. vec
                        else
                            plot(x_vec, allele_freq_vec(i,:), color_vec(i), 'linewidth', 2); % Plot allele freq. vec
                        end
                        title_str = 'Allele freq. dist.';
                        if(strcmp(scale{1}, 'linear'))
                            ylim([0 10]);
                        end
                         ylabel('Density');
                    case 2 % allele frequency cumulative
                        if(two_side_flag)
                            plot(x_vec, cum_allele_freq_vec2(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        else
                            plot(x_vec, cum_allele_freq_vec(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        end
                        title_str = 'Cum. allele freq. dist.';
                         ylabel('Cumulative');
                    case 3 % variance explained density
                        if(two_side_flag)
                            plot(x_vec, var_freq_vec2(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        else
                            plot(x_vec, var_freq_vec(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        end
                        title_str = 'Var. explained distribution';
                        if(strcmp(scale{1}, 'linear'))
                            ylim([0 10]);
                        end
                         ylabel('Density');
                    case 4 % variance explained cumulative
                        if(two_side_flag)
                            plot(x_vec, cum_var_vec2(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        else
                            plot(x_vec, cum_var_vec(i,:), color_vec(i), 'linewidth', 2); % Plot contribution to variance by freq. vec
                        end
                        title_str = 'Cum. Var. explained distribution';
                        ylim([0 1]);
                         ylabel('Cumulative');
                        %                plot(x_vec, cum_var_2side_vec, '.'); % Plot contribution to variance by freq. vec
                        %                title_str = 'Cum. Var. explained (two-sided) distribution';
                end
            end % loop on s
            if(plot_ind == 2)
                legend(legend_vec, 4);
            end
            title(title_str); xlabel(x_label);
        end % loop on sub-plot
        my_saveas(gcf, fullfile(spectrum_figs_dir, ...
            ['allele_freq_spectrum_' sup_title_str '_' scale{1} '_scale']), format_fig_vec);
    end % loop on scale
end % loop on two side flag
