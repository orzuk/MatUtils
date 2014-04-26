% Here plot just .. (temporary function) 
function FisherwrightPlotResults2(het_struct, s_vec, fisher_wright_output_dir, tmp_dir)

AssignGeneralConstants;
plot_inds = 3:5;

for normalized_flag = 0:1
    figure;
    for j=1:length(plot_inds)
        %                                    figure;
        for k= [1 2:2:length(het_struct)] % loop on s
            %                                         semilogx(het_struct{k}.plot_x_vec{2}{plot_inds(j)}, ...
            %                                             het_struct{k}.plot_y_vec{2,normalized_flag+1}{plot_inds(j)}, ...
            %                                             [color_vec(k) symbol_vec{ceil(k/6)}], 'linewidth', 2);
            semilogx(het_struct{k}.plot_x_vec{2}{plot_inds(j)}, ...
                het_struct{k}.plot_y_vec{2,normalized_flag+1}{plot_inds(j)}, ...
                [color_vec(ceil((k+1)/2)) symbol_vec{j}], 'linewidth', 2);
            
            hold on;
        end
        xlabel('Derived Allele Frequency'); ylabel('Heterozygosity (cumulative)');
        %                                    title(str2word(',', het_struct{1}.legend_vec{plot_inds(j)}, 1) );
        %                                    legend([repmat('|s|=', length(s_vec), 1) num2str(-s_vec', 3)], 2);
        legend([repmat('|s|=', length(s_vec([1 2:2:end])), 1) num2str(-s_vec([1 2:2:end])', 3)], 2);
        
        if(normalized_flag)
            xlim([10^(-7) 1]); norm_str = '_normalized';
        else
            xlim([10^(-5) 1]); norm_str = '';
        end
        
        %                                 my_saveas(gcf, fullfile(fisher_wright_output_dir, ...
        %                                     ['heterozygosity_cumulative_all_s_' ...
        %                                     str2word(',', het_struct{1}.legend_vec{plot_inds(j)}, 1) norm_str]), 'pdf');
    end % loop on which distribution to plot
    title('Solid - expansion-end,  Dotted - old equilibrium (N=10^4), Dashed - new equilibium (N=1.5*10^6)');
    my_saveas(gcf, fullfile(fisher_wright_output_dir, ...
        ['heterozygosity_cumulative_all_s_three_distributions' ...
        norm_str]), 'pdf');
    
end % loop on normalized flag



for freq = [0.005 0.01]
    S_mat = zeros(length(plot_inds), length(het_struct));
    for i=1:length(plot_inds)
        for j=1:length(het_struct)
            [~, I] = min( abs( het_struct{j}.plot_x_vec{2}{plot_inds(i)} - freq ) );
            S_mat(i, j) = het_struct{j}.plot_y_vec{2,2}{plot_inds(i)}(I);
        end
    end
    S_mat2 = [abs(s_vec') S_mat'];
    R = [{'s', 'expansion-end', 'old-equilibrium(N=10^4)', 'new-equilibrium(N=1.5*10^6)'}' ...
        num2cell(S_mat2)']';
    
end % loop on freq.
savecellfile(R, ...
    fullfile(fisher_wright_output_dir, 'cumulative_normalized_log', [tmp_dir '_mathematica_expansion_only.txt']));



