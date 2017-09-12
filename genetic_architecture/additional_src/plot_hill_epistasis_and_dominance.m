% function plot_hill_epistasis_and_dominance()

figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/two_locus_model';
a = 1;
res = 0.001;
p_vec = res:res:1-res;
q_vec = p_vec; n = length(p_vec);



for i=1:2
    if(i == 2)% Plot a 2-d
        p_vec = repmat(p_vec, n, 1); q_vec = p_vec';
    end
    Het_p_vec = 2.*p_vec.*(1-p_vec);
    Het_q_vec = 2.*q_vec.*(1-q_vec);
    
    
    V_A = a^2 .* (Het_p_vec + Het_q_vec - 4.*Het_p_vec.*Het_q_vec);
    V_AA = a^2 .* Het_p_vec .* Het_q_vec;
    V_G = a^2 .* (Het_p_vec + Het_q_vec - 3 .* Het_p_vec .* Het_q_vec);
    
    
    if(i==1) % plot 1-d
            figure; hold on;
        plot(p_vec, V_AA ./ V_G, 'linewidth', 2);
        plot(p_vec, V_A ./ V_G, 'g', 'linewidth', 2);
        plot(p_vec, V_G, 'r', 'linewidth', 2);
%        plot(p_vec, (V_AA + V_A) ./ V_G, 'k'); 
        xlabel('p,q'); ylabel('Variance');
        legend({'$\frac{V_{AA}}{V_G}$', '$\frac{V_{A}}{V_G}$', '$V_G$'}, 'location', 'east', 'interpreter', 'latex'); % legend
        legend('boxoff');
        title([repmat(' ', 1, 80) 'supp-fig9'], 'fontsize', 16, 'fontweight', 'bold');
        my_saveas(gcf, fullfile(figs_dir, 'two_locus_epistatic_variance'), {'epsc', 'pdf', 'fig', 'jpg'});
    else
            figure; mesh(p_vec, q_vec, V_AA ./ V_G); hold on;
        mesh(p_vec, q_vec, V_A ./ V_G);
        surf(p_vec, q_vec, V_G);
        
        legend({'$\frac{V_{AA}}{V_G}$', '$\frac{V_{A}}{V_G}$', '$V_G$'}, 'interpreter', 'latex'); % legend
        my_saveas(gcf, fullfile(figs_dir, 'two_locus_epistatic_variance_2d'), {'epsc', 'pdf', 'fig', 'jpg'});
        
    end
    
end
