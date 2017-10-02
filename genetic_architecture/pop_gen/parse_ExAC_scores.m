% Parse constraint scores computed by ExAC
%function parse_ExAC_scores() 

Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
scores_data_files = {'/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'}; %  ... % NEW! exome aggregation data (much richer!!!)

for i=1:length(scores_data_files)
    if(~exist(file_name_to_mat(fullfile(spectrum_data_dir, scores_data_files{i})), 'file'))
        R3 = ReadDataFile(fullfile(spectrum_data_dir, scores_data_files{i}), []); 
    else
        R3 = load(file_name_to_mat(fullfile(spectrum_data_dir, scores_data_files{i}))); 
    end
    T = loadcellfile(fullfile(spectrum_data_dir, scores_data_files{i}));
end
    

% Plot constraint of missense vs. constraint on LOF 
figure; hold on; 
plot(R.mis_z, R.lof_z, '.'); xlabel('missense-Z'); ylabel('LOF-Z'); 
xlim([-10 15]); ylim([-10 15]); 
xlabel('Missense Z-score'); ylabel('LOF Z-score'); 
my_saveas(gcf, 'ExAC_Z_scores', {'jpg', 'pdf', 'epsc'}); 

figure; hist ( log(abs(R.mis_z ./ R.lof_z)), 100);
figure; hist ( (R.mis_z - R.lof_z), 100);

s_vec = normcdf((R.mis_z + R.lof_z) / 2); 
s_vec = min(0.1, max(0, 10.^(-6-log(s_vec)))); 

alpha_vec = abs(R.mis_z) ./ (abs(R.mis_z) + abs(R.lof_z)); 


figure; semilogx(s_vec, alpha_vec, '.'); 
xlabel('s'); ylabel('\alpha'); 
my_saveas(gcf, 'ExAC_alpha_s', {'jpg', 'pdf', 'epsc'}); 





