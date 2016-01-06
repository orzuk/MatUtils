% Plot the relation between GRR and variance explained
function GRR_and_h_liab_plot()

AssignGeneralConstants();

f_vec = [0.001 0.01 0.1 0.2 0.4];

h_liab_vec = 0:0.001:0.2;
mu = 0.1; 
figure; hold on;
for i=1:length(f_vec)
    GRR_vec = heritability_to_genetic_relative_risk(h_liab_vec, 'liability', f_vec(i), mu);
    plot(h_liab_vec, GRR_vec, [color_vec(i)]); xlabel('h_{liab}^2'); ylabel('GRR');
    
end

legend(num2str_cell(num2cell(f_vec)))

