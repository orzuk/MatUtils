% Plot results of FisherWright simulation
%
% Input: 
% freq_struct, absorption_struct, simulation_struct, ...
%    N_vec, expansion_factor, s, mu, num_bins, init_str, iters, ...
%    fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, all_s_output_file, mathematica_flag
%
% Output: 
% plot_x_vec 
% plot_y_vec 
% legend_vec 
% title_str
%
function  [plot_x_vec plot_y_vec legend_vec title_str] = ...
    FisherWrightPlotResults(freq_struct, absorption_struct, simulation_struct, ...
    N_vec, expansion_factor, s, mu, num_bins, init_str, iters, ...
    fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, all_s_output_file, mathematica_flag)


AssignGeneralConstants;
if(~exist('num_generations', 'var') || isempty(num_generations))
    num_generations = length(freq_struct{1}.x_vec)-1;
end
if(~exist('N', 'var') || isempty(N))
    N = N_vec{1}(1);
end


format_fig_vec = {}; % {'pdf'};
new_sim = 1;
num_plots = 6; plot_vec = max(1, floor(linspace(1, num_generations, num_plots)))
legend_vec = [repmat('K=', num_plots, 1) num2str(plot_vec')]; %  repmat(', H=', num_plots, 1) num2str(total_het_vec(plot_vec))];

switch init_str
    case 'newly_born'
        print_str = 'newly_born';
    case 'equilibrium'
        print_str = 'old'
end

for normalized_flag = 1:1 % 0:1
    figure; hold on; ctr=1; max_quantile=0;
    for i=plot_vec  % plot how heterozygosity distribution evolves
        if(new_sim)
            %            freq_struct{1}.x_vec = (0:2*N_vec{1}(i)) ./ (2*N_vec{1}(i)); % assume we get freq_struct{1}.x_vec as input
        else
            freq_struct{1}.x_vec = (1:2*N_vec{1}(i)-1) ./ (2*N_vec{1}(i));
        end
        if(normalized_flag)
            norm_str = '_normalized';
        else
            norm_str = ' (not normalized)';
        end
        
        if(iscell(freq_struct{1}.het_vec))
            max_quantile_ind = find(cumsum(freq_struct{1}.het_vec{i}) >= sum(freq_struct{1}.het_vec{i}) * 0.9999, 1);
            max_quantile = max(max_quantile, freq_struct{1}.x_vec{i}(max_quantile_ind));
            
            if(normalized_flag)
                freq_struct{1}.het_vec{i} = normalize_hist(reshape(freq_struct{1}.x_vec{i}, size(freq_struct{1}.het_vec{i})), freq_struct{1}.het_vec{i}); % normalization (ruins un-normalized!)
            end
            plot(freq_struct{1}.x_vec{i}(2:end-1), freq_struct{1}.het_vec{i}(2:end-1), color_vec(ctr), 'linewidth', 2);
        else
            %        freq_struct{1}.x_vec = (1:N_vec{1}(i+1)-1) ./ N_vec{1}(i+1);
            if(ctr == length(plot_vec)) % plot only once (last one)
                plot(freq_struct{1}.x_vec, normalize_hist(freq_struct{1}.x_vec, freq_struct{1}.het_vec), color_vec(ctr), 'linewidth', 2);
            end
        end
        ctr=ctr+1;
    end % loop on indices to plot
    ylim([0 max(10^(-20), min(10, 1.1*max(max(cell2vec(freq_struct{1}.het_vec)))))]);
    xlim([0, min(1, max_quantile)]); % don't show stuff that's too high
    xlabel('Derived Allele Freq.'); ylabel(str2title(['Heterozygosity ' norm_str])); legend(legend_vec); % num2str(plot_vec'));
    title_str = ['N_0=' num2str(N*1) str2title([', expan. = 1+' num2str(unique(vec2row(expansion_factor))-1, 3), ...
        ', K=' num2str(num_generations*1) ', s=' num2str(s,3) ', \mu=' num2str(mu,2) ...
        ', iters (old=' num2str(simulation_struct{1}.old_num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(1)) ...
        ', new=' num2str(simulation_struct{1}.new_num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(2)) '), ' print_str])];
    title(title_str); % str2title(title_str) legend('200');
    my_saveas(gcf, [fisher_wright_output_file '_heterozygosity' norm_str], format_fig_vec);
end % loop on normalized flag

% % % figure; hold on; ctr=1;
% % % for i=plot_vec  %num_generations  % plot how distribution of #alleles evolves
% % %     if(~new_sim)
% % %         freq_struct{1}.x_vec = (1:2*N_vec{1}(i)-1) ./ (2*N_vec{1}(i));
% % %     end
% % %     if(iscell(freq_struct{1}.p_vec))
% % %         plot(freq_struct{1}.x_vec{i}(2:end-1), freq_struct{1}.p_vec{i}(2:end-1), color_vec(ctr), 'linewidth', 2);
% % %     else
% % %         %        freq_struct{1}.x_vec = (1:N_vec{1}(i+1)-1) ./ N_vec{1}(i+1);
% % %         if(ctr == length(plot_vec)) % plot only once
% % %             plot(freq_struct{1}.x_vec, freq_struct{1}.p_vec(2:end-1), color_vec(ctr), 'linewidth', 2);
% % %         end
% % %     end
% % %     ctr=ctr+1;
% % % end
% % % xlabel('Derived Allele Freq.'); ylabel('Frac. Alleles'); legend(legend_vec); % num2str(plot_vec'));
% % % title(str2title(title_str));
% % % my_saveas(gcf, [fisher_wright_output_file '_frac_alleles'], format_fig_vec);


figure;  ctr=1;
for i=plot_vec  %num_generations  % plot how distribution of #alleles evolves
    if(~new_sim)
        freq_struct{1}.x_vec = (1:2*N_vec{1}(i)-1) ./ (2*N_vec{1}(i));
    end
    if(iscell(freq_struct{1}.p_vec))
        semilogy(freq_struct{1}.x_vec{i}(2:end-1), normalize_hist(vec2row(freq_struct{1}.x_vec{i}(2:end-1)), vec2row(freq_struct{1}.p_vec{i}(2:end-1))), ...
            color_vec(ctr), 'linewidth', 2);
    else
        %        freq_struct{1}.x_vec = (1:N_vec{1}(i+1)-1) ./ N_vec{1}(i+1);
        if(ctr == length(plot_vec)) % plot only once
            semilogy(freq_struct{1}.x_vec, normalize_hist(freq_struct{1}.x_vec, freq_struct{1}.p_vec(2:end-1)), color_vec(ctr), 'linewidth', 2);
        end
    end
    hold on; ctr=ctr+1;
end
xlabel('Derived Allele Freq.'); ylabel('Frac. Alleles Density'); legend(legend_vec); % num2str(plot_vec'));
title(str2title(title_str));
my_saveas(gcf, [fisher_wright_output_file '_frac_alleles_normalized'], format_fig_vec);


figure; hold on; ctr=1;
for i=plot_vec
    q_sorted = sort(simulation_struct{1}.q(:,i));
    het_sorted = cumsum(2 .* q_sorted .* (1-q_sorted)); het_sorted = het_sorted ./ het_sorted(end);
    plot(q_sorted  , het_sorted, color_vec(ctr)); ctr=ctr+1;
end
title(str2title(title_str));
xlabel('Derived Allele Freq.'); ylabel('Heterozygosity (cumulative)');  legend(legend_vec);
my_saveas(gcf, [fisher_wright_output_file '_heterozygosity_cumulative'], format_fig_vec);





% New: plot both old distribution + new distribution
if(~mathematica_flag)
    display_x_vec = (0:num_bins) ./ num_bins;
else
    display_x_vec = num_bins;
end
display_equilibrium_x_vec = (0:20*N_vec{1}(1)) ./ (20*N_vec{1}(1));


old_theta = 4.*N_vec{1}(1).*mu; % get theta at start
new_theta = 4.*N_vec{1}(num_generations).*mu; % get theta at end

old_equilibrium_het_vec_density = old_theta.*2.*exp(allele_freq_spectrum(display_equilibrium_x_vec, s, N_vec{1}(1), 0, 'log', 1)); % Add factor two!!!
new_equilibrium_het_vec_density = new_theta.*2.*exp(allele_freq_spectrum(display_equilibrium_x_vec, s, N_vec{1}(num_generations), 0, 'log', 1));
[~, old_equilibrium_het_vec_cumulative] = integral_hist(display_equilibrium_x_vec, old_equilibrium_het_vec_density, 1);
[~, new_equilibrium_het_vec_cumulative] = integral_hist(display_equilibrium_x_vec, new_equilibrium_het_vec_density, 1);

old_equilibrium_absolute_het = absorption_time_by_selection(abs(s), old_theta, N_vec{1}(1), 0, 1, 'var');
new_equilibrium_absolute_het = absorption_time_by_selection(abs(s), new_theta, N_vec{1}(num_generations), 0, 1, 'var');


switch init_str % Use these plots
    case 'newly_born' % we currently assume only newly_born flag is available 
        [plot_x_vec plot_y_vec legend_vec dir_str] = ...
            fisher_right_heterozygosity_plot_internal( ...
            freq_struct, absorption_struct, simulation_struct, ...
            N_vec, s, mu, num_bins, title_str, display_x_vec, display_equilibrium_x_vec, ...
            old_equilibrium_het_vec_density, new_equilibrium_het_vec_density,  ...
            old_equilibrium_het_vec_cumulative, new_equilibrium_het_vec_cumulative, ...
            old_equilibrium_absolute_het, new_equilibrium_absolute_het, ...
            fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, mathematica_flag);
        
        if(mathematica_flag) % save results in Mathematica format
            R = []; j_ctr=1; % empty cell
            for j=3 %% 1:5 % loop on all indices
                %                R{j_ctr} = ['{ ' num2str_delim(plot_x_vec{2}{j}, ', ', 3) ' }']; j_ctr = j_ctr+1; % bin values
                %                R{j_ctr} = ['{ ' num2str_delim(plot_y_vec{2,2}{j}, ', ', 3) ' }']; j_ctr = j_ctr+1; % y values
                R{j_ctr} = num2str_delim(plot_y_vec{2,2}{j}, ' ', 3); j_ctr = j_ctr+1; % y values
                
                %                R{j_ctr} = ''; j_ctr=j_ctr+1; % space
            end
            savecellfile(R', ...
                fullfile(fisher_wright_output_dir, 'cumulative_normalized_log', [tmp_dir '_mathematica_expansion_only.txt']));
        end % if mathematica flag
        
end % switch init_str


% Plot heterozygosity as function of number of generations which passed 
figure; hold on; compute_mode_legend_vec = cell(length(freq_struct), 1);
for i=1:length(freq_struct)
    plot( 1:num_generations-1, ...
        freq_struct{i}.old_total_het_at_each_generation_vec(1:num_generations-1), color_vec(i), 'linewidth', 2);
    compute_mode_legend_vec{i} = freq_struct{i}.compute_mode;
end
title(['OLD Heterozygosity at each #generations ' title_str]);
%compute_mode_legend_vec{end+1} = 'equilibrium'; 
%plot(1:num_generations+1, repmat(2*N*mu, num_generations+1, 1), 'k--'); % should be 4*N*mu
xlabel('Time (generations)'); ylabel('Total Old Heterozygosity Per-Site'); legend(compute_mode_legend_vec); 

figure; hold on; 
for i=1:length(freq_struct)
    plot( 1:num_generations-1, ...
        freq_struct{i}.total_het_at_each_generation_vec(1:num_generations-1), color_vec(i), 'linewidth', 2);
end
title(['NEW Heterozygosity at each #generations ' title_str]);

%plot(1:num_generations+1, repmat(2*N*mu, num_generations+1, 1), 'k--'); % should be 4*N*mu
xlabel('Time (generations)'); ylabel('Total New Heterozygosity Per-Site'); legend(compute_mode_legend_vec); 



return;




if(s == 0) % New plot. % Plot approximate total fraction of new vs. old heterozygosity 
    figure; 
    s_vec = -[0 logspace(-6, -1, 1001)];
    old_het_vec = vec2row(absorption_time_by_selection(abs(s_vec), 4*N_vec{1}(1)*mu, N_vec{1}(1), 0, 1, 'var')) .* ...
        (1+s_vec).^num_generations;
    new_het_vec = 2*mu* (1 - (1+s_vec).^num_generations) ./ (-s_vec);
    new_over_old_het_vec = new_het_vec ./ (new_het_vec + old_het_vec);
    
    semilogx(abs(s_vec), new_over_old_het_vec);
    xlabel('|s|'); ylabel('New het. / Total (New+Old) het.');
    title([title_str(1:max(find(title_str == ',', 3))) ' new het. proportion after expansion']);
    %    hold on; semilogx(abs(s_vec), repmat(num_generations/N, length(s_vec), 1), 'k--')
    my_saveas(gcf, fullfile(fisher_wright_output_dir, 'new_vs_old_het'), format_fig_vec);
end




if(length(total_het_vec) > 1)
    figure; plot(1:length(total_het_vec)-1, diff(total_het_vec)); title(['Change of Heterozygosity by #generations ' title_str]);
    xlabel('Time (generations)'); ylabel('\Delta H ');
    
    
    figure; hold on; plot(1:length(total_het_vec)-1, total_het_vec(2:end)./total_het_vec(1:end-1)-1);
    plot(1:num_generations, repmat(-1./(2.*N), num_generations, 1), 'k--');
    title(['Relative Change of Heterozygosity by #generations ' title_str]);
    xlabel('Time (generations)'); ylabel('\frac{\Delta H}{H}', 'interpreter', 'latex');
end



figure; hold on; plot(1:length(frac_alleles_kept_vec), prob_site_polymorphic .* frac_alleles_kept_vec);
% if(mu == 0)
%     plot(1:num_generations, prob_site_polymorphic .* frac_alleles_kept_vec, 'k--'); % repmat(2*N*mu*absorb_time_vec(2), num_generations, 1), 'k--'); % plot asymptote
% end
title(['Prob. Polymorphic by #generations ' title_str]);
xlabel('Time (generations)'); ylabel('Prob. Polymorphic');

if(new_sim)
    freq_struct{1}.x_vec = (0:2*N) ./ (2*N);
else
    freq_struct{1}.x_vec = (1:2*N-1) ./ (2*N);
end


% figure; hold on;
% if(length(absorb_time_vec) == 2*N-1)
%     absorb_time_vec = [0 absorb_time_vec' 0]';
%     fixed_time_vec = [0 fixed_time_vec' 0]';
%     loss_time_vec = [0 loss_time_vec' 0]';
% end
% plot(freq_struct.x_vec, absorb_time_vec); plot(freq_struct.x_vec, fixed_time_vec, 'g'); plot(freq_struct.x_vec, loss_time_vec, 'r');
% title('Expected absorption Time by Frequency');
% legend('time-to absorption', 'time-to-fixation', 'time-to-loss');
% xlabel('Derived Allele Freq.'); ylabel('Mean Time to Fixation/Loss (#generations)');



% figure; plot(absorb_time_vec) % Final plot
% N_vec = N .* expansion_factor .^ (0:num_generations);
% LOH_vec = 1-1 ./ (2.*N_vec);
% total_LOH_vec = prod(LOH_vec)
%
% approx_total_LOH_vec = 1 - (1/(2*N)) * expansion_factor / (expansion_factor-1) % approximation


% Plot heterozygosity in various ways (log and linear scale)
% Problem: there are two different 'plot_x_vec' and 'plot_y_vec'
%
% Input:
%
% Output:
% plot_x_vec - cell array of x-values to plot
% plot_y_vec - cell array of y-values (probability distributions) to plot
% legend_vec -
% dir_str - where to save results (?)
%
function [plot_x_vec plot_y_vec legend_vec dir_str] = ...
    fisher_right_heterozygosity_plot_internal( ...
    freq_struct, absorption_struct, simulation_struct, ...
    N_vec, s, mu, num_bins, title_str, display_x_vec, display_equilibrium_x_vec, ...
    old_equilibrium_het_vec_density, new_equilibrium_het_vec_density,  ...
    old_equilibrium_het_vec_cumulative, new_equilibrium_het_vec_cumulative, ...
    old_equilibrium_absolute_het, new_equilibrium_absolute_het, ...
    fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, mathematica_flag)



AssignGeneralConstants;
if(~exist('num_generations', 'var') || isempty(num_generations))
    num_generations = length(freq_struct{1}.x_vec)-1;
end
format_fig_vec = {}; % {'pdf'};

if(~mathematica_flag)
    num_bins = min(num_bins, 2*N_vec{1}(num_generations));
end


old_weight = zeros(length(freq_struct),1); new_weight = old_weight;
max_quantile_ind = old_weight; max_quantile = old_weight;
combined_x = cell(length(freq_struct), 1); combined_het = combined_x; 

for j=1:length(freq_struct)
    switch freq_struct{j}.compute_mode
        case 'simulation'
            old_weight(j) = absorption_struct{1}.prob_site_polymorphic_at_equilibrium; %  / ...
%                simulation_struct{1}.old_num_simulated_polymorphic_alleles_vec(1);
            new_weight(j) = 1; %  / simulation_struct{1}.new_num_simulated_polymorphic_alleles_vec(1);
        case 'numeric'
            old_weight(j)=1; new_weight(j)=1;
    end
    
    display_old_het_vec{j} = weighted_hist(freq_struct{j}.all_old_x_vec ./ (2*N_vec{1}(num_generations)), ...
        old_weight(j) .* freq_struct{j}.all_old_het_vec, display_x_vec);
    display_new_het_vec{j} = weighted_hist(freq_struct{j}.all_new_x_vec  ./ (2*N_vec{1}(num_generations)), ...
        new_weight(j) .* freq_struct{j}.all_new_het_vec, display_x_vec);
    [display_combined_x_vec{j} display_combined_het_vec{j}] = ...
        sum_hist(display_x_vec, display_old_het_vec{j}, display_x_vec, display_new_het_vec{j}, [], 0); % don't normalize sum
    
    
    max_quantile_ind(j) = find(cumsum(freq_struct{j}.all_old_het_vec) >= sum(freq_struct{j}.all_old_het_vec) * 0.9999, 1);
    max_quantile(j) = freq_struct{j}.all_old_x_vec(max_quantile_ind(j)) /  (2*N_vec{1}(num_generations));
    
    max_quantile_ind(j) = find(old_equilibrium_het_vec_cumulative >= old_equilibrium_het_vec_cumulative(end) * 0.999, 1);
    max_quantile(j) = max(max_quantile(j), display_equilibrium_x_vec(max_quantile_ind(j)));
    max_quantile_ind(j) = find(new_equilibrium_het_vec_cumulative >= new_equilibrium_het_vec_cumulative(end) * 0.999, 1);
    max_quantile(j) = max(max_quantile(j), display_equilibrium_x_vec(max_quantile_ind(j)));
    
    % Combine new and old heterozygosities together 
    combined_x{j} = [freq_struct{j}.all_old_x_vec freq_struct{j}.all_new_x_vec] ./ (2*N_vec{1}(num_generations));
    combined_het{j} = [old_weight(j) .* freq_struct{j}.all_old_het_vec new_weight(j) .* freq_struct{j}.all_new_het_vec];
    [combined_x{j} sort_perm] = sort(combined_x{j});
    combined_het{j} = cumsum(combined_het{j}(sort_perm));  
    
end % loop on computation method


% slash_inds = find(fisher_wright_output_file == '\');
% fisher_wright_output_file_same_dir = fisher_wright_output_file;
% fisher_wright_output_file_same_dir(slash_inds(end)) = '_';

plot_symbol_vec = cellstr(repmat('-', 3*length(freq_struct) + 2, 1));
absolute_het_vec = [];
for j=1:length(freq_struct)
    absolute_het_vec(3*j-2) = sum(old_weight(j) .* freq_struct{j}.all_old_het_vec)
    absolute_het_vec(3*j-1) = sum(new_weight(j) .* freq_struct{j}.all_new_het_vec);
    absolute_het_vec(3*j) = combined_het{j}(end); 

    switch freq_struct{j}.compute_mode
        case 'numeric'
            plot_symbol_vec{3*j-2} = '--'; plot_symbol_vec{3*j-1} = '--'; plot_symbol_vec{3*j} = '--';                        
        case 'simulation'                        
            plot_symbol_vec{3*j-2} = ':'; plot_symbol_vec{3*j-1} = ':'; plot_symbol_vec{3*j} = ':';                        
    end


end
plot_symbol_vec{end-1} = '-'; plot_symbol_vec{end} = '-'; % last two for equilibrium (analyitc) 


absolute_het_vec = [absolute_het_vec ...
    old_equilibrium_het_vec_cumulative(end) new_equilibrium_het_vec_cumulative(end)];


new_to_total_het = new_equilibrium_het_vec_cumulative(end) ./ ... % fraction of heterozygosity due to new alleles born after the expansion
    (old_equilibrium_het_vec_cumulative(end) + new_equilibrium_het_vec_cumulative(end));

nonzero_s = min(s, -0.0000000000000001); % avoid 0
for j=1:length(freq_struct)
    absolute_predicted_het_vec(3*j-2) = ... % heuristic prediction of heterozygosity
        old_equilibrium_absolute_het .* (1+s)^num_generations;
    absolute_predicted_het_vec(3*j-1) = ...
        2*mu* (1-(1+nonzero_s)^num_generations ) / (-nonzero_s);
    absolute_predicted_het_vec(3*j) = ... % independent of calculation method b
        old_equilibrium_absolute_het .* (1+nonzero_s)^num_generations +  2*mu* (1-(1+nonzero_s)^num_generations ) / (-nonzero_s);
end
absolute_predicted_het_vec = [absolute_predicted_het_vec ...
    old_equilibrium_absolute_het .* (1+nonzero_s)^num_generations +  2*mu* (1-(1+nonzero_s)^num_generations ) / (-nonzero_s) ... % combined
    old_equilibrium_absolute_het new_equilibrium_absolute_het];


plot_x_vec = cell(2,1); plot_y_vec = cell(2,1);
for log_x_flag = 1:1 % use log-scale on x-axis
    for normalized_flag = 0:1
        for cumulative_flag = 1:1 % 0 - plot density. 1 - plot cumulative
            figure; % plot heterozygosity density
            
            if(~cumulative_flag)
                cum_str = 'density'; legend_loc = 1;
                plot_x_vec{cumulative_flag+1} = {display_x_vec, display_x_vec, display_combined_x_vec, ...
                    display_equilibrium_x_vec, display_equilibrium_x_vec};
            else % plot cumulatives
                plot_x_vec{cumulative_flag+1}=[];
                cum_str = 'cumulative'; legend_loc = 2;
                for j=1:length(freq_struct)
                    plot_x_vec{cumulative_flag+1}{3*j-2} = freq_struct{j}.all_old_x_vec ./ (2*N_vec{1}(num_generations));
                    plot_x_vec{cumulative_flag+1}{3*j-1} = freq_struct{j}.all_new_x_vec ./ (2*N_vec{1}(num_generations));
                    plot_x_vec{cumulative_flag+1}{3*j} = combined_x{j};
                end
                plot_x_vec{cumulative_flag+1} = [plot_x_vec{cumulative_flag+1} ...
                    {display_equilibrium_x_vec, display_equilibrium_x_vec}];
            end
            
            
            if(~cumulative_flag) % first compute un-normalized curves
                plot_y_vec{cumulative_flag+1,normalized_flag+1} = ...
                    {display_old_het_vec, display_new_het_vec, display_combined_het_vec, ...
                    old_equilibrium_het_vec_density, new_equilibrium_het_vec_density};
            else
                % New! allow comparison of different computation modes !!!
                plot_y_vec{cumulative_flag+1,normalized_flag+1} = [];
                for j=1:length(freq_struct)
                    plot_y_vec{cumulative_flag+1,normalized_flag+1}{3*j-2} = cumsum(old_weight(j) .* freq_struct{j}.all_old_het_vec);
                    plot_y_vec{cumulative_flag+1,normalized_flag+1}{3*j-1} = cumsum(new_weight(j) .* freq_struct{j}.all_new_het_vec);
                    plot_y_vec{cumulative_flag+1,normalized_flag+1}{3*j} = combined_het{j};
                end
                plot_y_vec{cumulative_flag+1,normalized_flag+1} = [plot_y_vec{cumulative_flag+1,normalized_flag+1} ...
                    {old_equilibrium_het_vec_cumulative, new_equilibrium_het_vec_cumulative}];
            end
            if(~normalized_flag)
                norm_str = ' (not normalized)'; save_str = '';
            else % normalized data
                norm_str = '_normalized'; save_str = '_normalized';
                if(~cumulative_flag) % plot density
                    plot_y_vec{cumulative_flag+1,normalized_flag+1} = ...
                        {display_old_het_vec ./ integral_hist(display_x_vec, display_old_het_vec), ...
                        display_new_het_vec ./ integral_hist(display_x_vec, display_new_het_vec), ...
                        display_combined_het_vec ./ integral_hist(display_combined_x_vec, display_combined_het_vec), ...
                        old_equilibrium_het_vec_density ./ integral_hist(display_equilibrium_x_vec, old_equilibrium_het_vec_density), ...
                        new_equilibrium_het_vec_density ./ integral_hist(display_equilibrium_x_vec, new_equilibrium_het_vec_density)};
                else % plot cumulatives
                    for j=1:length(plot_y_vec{cumulative_flag+1,normalized_flag+1})
                        plot_y_vec{cumulative_flag+1,normalized_flag+1}{j} = ...
                            plot_y_vec{cumulative_flag+1,normalized_flag+1}{j} ./ plot_y_vec{cumulative_flag+1,normalized_flag+1}{j}(end);
                    end
                    %                     {cumsum(freq_struct{1}.all_old_het_vec) ./ sum(freq_struct{1}.all_old_het_vec), cumsum(freq_struct{1}.all_new_het_vec) ./ sum(freq_struct{1}.all_new_het_vec), ...
                    %                         combined_het ./ combined_het(end), ...
                    %                         old_equilibrium_het_vec_cumulative ./ old_equilibrium_het_vec_cumulative(end), ...
                    %                         new_equilibrium_het_vec_cumulative ./ new_equilibrium_het_vec_cumulative(end)};
                end
            end % if normalized flag
            
            for j=1:length(plot_x_vec{cumulative_flag+1}) % loop on different plots
                if(log_x_flag)
                    log_str = '_log'
                    semilogx(plot_x_vec{cumulative_flag+1}{j}, plot_y_vec{cumulative_flag+1,normalized_flag+1}{j}, ...
                        [color_vec(j) plot_symbol_vec{j}], 'linewidth', 2);
                else
                    log_str = '';
                    plot(plot_x_vec{cumulative_flag+1}{j}, plot_y_vec{cumulative_flag+1,normalized_flag+1}{j}, ...
                        [color_vec(j) plot_symbol_vec{j}], 'linewidth', 2);
                end
                hold on;
            end
            
            xlabel('Derived Allele Freq.'); ylabel(str2title(['Heterozygosity (' cum_str ') ' norm_str]));
            legend_vec = [];
            for j=1:length(freq_struct)
                legend_vec{3*j-2} = ['old-survived-' freq_struct{j}.compute_mode];
                legend_vec{3*j-1} = ['newly-born-' freq_struct{j}.compute_mode];
                legend_vec{3*j} = ['expansion-end-' freq_struct{j}.compute_mode];
            end
            legend_vec = [legend_vec {'old-equilibrium', 'new-equilibrium'}]';
            for j=1:length(legend_vec)
                legend_vec{j} = [legend_vec{j} ', H=' num2str(absolute_het_vec(j), 2) ...
                    ' (' num2str(absolute_predicted_het_vec(j), 2) ')'];
            end
            legend(legend_vec, legend_loc);
            title(title_str); % str2title(title_str)
            xlim([0 max(max_quantile)]);
            dir_str = [cum_str save_str log_str];
            my_saveas(gcf, fullfile(fisher_wright_output_dir, dir_str, tmp_dir), format_fig_vec);
            
        end % loop on cumulative flag
        
    end % loop on normalization flag
end % loop on log flag
