% Analyze IBD data from Barak Markus
function IBD_compute_statistics(IBD_file)

if(~exist(file_name_to_mat(IBD_file), 'file'))
    R = loadcellfile(IBD_file);
    save(file_name_to_mat(IBD_file), 'R'); % save as .mat
else
    load(file_name_to_mat(IBD_file));
end
IBD_vec = cell2mat(str2num_cell(R(:,7)));

IBD_vec = IBD_vec(IBD_vec < 0.4); % Remove siblings
figure; hold on; hist(IBD_vec, 500);
%title(['IBD Hist. Mean = ' num2str(mean(IBD_vec)) ', Std = ' num2str(std(IBD_vec))]);
xlabel('IBD'); ylabel('Freq.');
title(['IBD Hist. Mean = ' num2str(mean(IBD_vec), 3) ', Std = ' num2str(std(IBD_vec), 3)]);
