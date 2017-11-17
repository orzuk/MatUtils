
% function mutation_rates_analysis()

Assign24MammalsGlobalConstants; AssignRVASConstants;

M = load(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file_exons)); % hg18 (from single exons)
M2 = load(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file)); % from Samocha et al. (per gene) 

[common_genes, I, I2] = intersect(upper(M.gene_names), upper(M2.gene));

figure; 
for i=1:4
    switch i
        case 1
            x_vec = M.UniqueMutationRateTable(I, 1); y_vec = 10.^(M2.syn(I2)); class_str = 'synonymous';
        case 2
            x_vec = M.UniqueMutationRateTable(I, 2); y_vec = 10.^(M2.mis(I2)); class_str = 'missense';
        case 3
            x_vec = sum( M.UniqueMutationRateTable(I, 3:4), 2 );  y_vec = 10.^(M2.non(I2));  class_str = 'nonsense';
        case 4
            x_vec = sum( M.UniqueMutationRateTable(I, :), 2 );  y_vec = 10.^(M2.all(I2)); class_str = 'all';
    end
    subplot(2,2,i); plot(x_vec, y_vec, '*'); hold on;
    
    ratio(i) = round(median(x_vec ./ y_vec)); 
    beta  = polyfit(x_vec, y_vec, 1); plot(sort(x_vec), sort(x_vec) .* beta(1) + beta(2), 'r'); 
    title([class_str ', r=' num2str(corr(x_vec, y_vec)) ', \beta=' num2str(round(1/beta(1)))]); xlabel('HG18'); ylabel('Samocha');
end                        
     
