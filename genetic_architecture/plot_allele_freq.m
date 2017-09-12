% Plot allele frequency distrinutions at equilibrium
% Input: 
% 
function plot_allele_freq(s_vec, N, two_side_flag, log_flag, cum_flag, scale_mode, weight_flag)

AssignGeneralConstants;
eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
my_symbol_vec = {'-', '--'};

num_s = length(s_vec);

% allele_freq_cumulative(x, s, N, two_side_flag, scale_mode, weight_flag)

legend_vec = [repmat('s= -10^{', length(s_vec), 1) num2str(log10(abs(s_vec')),3) ...
    repmat('}', length(s_vec), 1)];
legend_vec = cellstr(legend_vec);
legend_vec = strrep_cell(legend_vec, ' ', '');
legend_vec{1} = 's= 0'; % fix s=0
x_vec = (0:(2*N)) ./ (2*N);

figure;
for i=1:length(s_vec);
    f_vec{i} = allele_freq_cumulative(x_vec, s_vec(i), N, two_side_flag, scale_mode, weight_flag);    
    semilogx(x_vec, f_vec{i}, [eric_color_vec(ceil(i/2)) my_symbol_vec{mod_max(i+1,2)}], 'linewidth', 2); hold on;    
end
%legend(legend_vec, 'location', 'best'); legend('boxoff');
legend(legend_vec, 'position', [0.76 0.09 0.16 0.4]); legend('boxoff');

%    legend('frac. null', 'frac mutations captured');
xlabel('Derived allele frequency f'); ylabel('Proportion of alleles below f, \Psi_s(f)');
ylim([0 1.01*max(max_cell(f_vec))]);
if(log_flag)
    xlim([10^(-4) 2]);
end
add_faint_grid(0.5);

