% Show how many loci are needed to explain heritability for varrious effect
% sizes
% 
% Input: 
% fig_outfile - where to save the output figures 
% 
function plot_num_loci_needed(fig_outfile)

AssignGeneralConstants;

mu_vec  = [0.0001 0.001 0.01 0.1]; % try different prevalences

%grr_vec = 10.^(-6:0.1:0);
%grr_vec = log(interp(exp(grr_vec), 10))

grr_vec = -4:0.01:0; % take a finer resolution. On exponential scale 
% grr_pow_vec = ([exp(1):0.01:exp(2)]); %1:0.01:3; % effect size
f_vec = 0.01:0.01:0.99; % allele frequency
f_vec = -4:0.01:0; % move f also to exponential scale 

for mu = mu_vec
    for i=1:length(grr_vec)
        [~, ~, ~, h_add, V_add, ...
            ~, ~, ~, ~, ~, h_liab h_liab_vec] = ...
            genetic_relative_risk_to_heritability(10.^f_vec, repmat(1+10.^(grr_vec(i)), length(f_vec), 1), mu); % compute heritability
        num_loci(i,:) = 1./h_liab_vec;
    end
    
    figure; [C,h] = contour(f_vec, grr_vec, log10(num_loci));
    clabel(C,h);
    ylabels = get(gca, 'YTickLabel');
    ylabels = num2str(1+10.^(str2num(ylabels)), 6);
    set(gca, 'YTickLabel', ylabels);
    xlabels = get(gca, 'XTickLabel');
    xlabels = num2str(10.^(str2num(xlabels)), 6);
    set(gca, 'XTickLabel', xlabels);
    
    % imagesc(log(num_loci)); colorbar;
    
    title(['Num loci needed to explain all variance. Prevalence = ' num2str(mu)]);
    xlabel('Risk-Allele-Frequency');  ylabel('Genetic-Relative-Risk');
    my_saveas(gcf, [fig_outfile '_prevalence_' num2str(mu)], format_fig_vec); % save figs
end
