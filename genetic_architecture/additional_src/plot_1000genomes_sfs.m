% Plot SFS from 1000 genomes
R = loadcellfile('C:\users\user\Google Drive\lab\frequencies_matrix_file.txt');


for i=1:5
    R_mat{i} = R(i+4,2:end);
    R_lable{i} = R{i+4,1}
    num_snps(i) = find( isempty_cell(R_mat{i}), 1)-1;
    R_mat{i} = cell2num(R_mat{i}(1:num_snps(i)));
end

num_categories = length(R_mat); % 5

% plot cumulatives
log_flag = 1;
figure; % hold on;
for i=1:num_categories
    if(log_flag)
        semilogx(R_mat{i}, (1:length(R_mat{i})) ./ length(R_mat{i}), color_vec(i));     hold on;
    else
        plot(R_mat{i}, (1:length(R_mat{i})) ./ length(R_mat{i}), color_vec(i));     hold on;
    end
end
legend(R_lable, 4);


figure; % hold on;
for i=1:num_categories
    R_mat_cum{i} = cumsum(R_mat{i}) ./ sum(R_mat{i}); 
    if(log_flag)
        semilogx(R_mat{i}, R_mat_cum{i}, color_vec(i));     hold on;
    else
        plot(R_mat{i}, R_mat_cum{i}, color_vec(i));     hold on;
    end
end
legend(R_lable, 4);


