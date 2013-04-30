% function build_reed_solomon_code()

reed_solomon_dir = '../../compressed_sensing/IMA/ReedSolomon'; % where to save the results 
m = 4; % Set field size and irreducible polynomial
p = bin2dec('10011'); % Primitive polynomial is x^4 + x + 1  
% p = bin2dec('11001'); % Primitive polynomial is x^4 + x^3 + 1

% m = 2; p = bin2dec('111'); % Primitive polynomial is x^2 + x + 1

q = 2^m; % Field size 
deg = 1; % degree of polynomials

layered_code_mat = build_layered_reed_solomon_code(q, deg, p); % Generate code

% Draw figures and save to file
figure; colormap(bw); imagesc(layered_code_mat); % colorbar;
my_saveas(gcf, fullfile(reed_solomon_dir, ['reed_solomon_code_matrix_q_' num2str(q) '_deg_' num2str(deg)  '_poly_'  dec2bin(p)]), ...
    {'epsc', 'jpg', 'pdf'});
savecellfile(num2cell(layered_code_mat), ...
    fullfile(reed_solomon_dir, ['reed_solomon_code_matrix_q_' num2str(q) '_deg_' num2str(deg) '.txt']));



pools_overlap_mat = layered_code_mat * layered_code_mat';
figure; imagesc(pools_overlap_mat - diag(diag(pools_overlap_mat))); colorbar; 
title('Overlap of pools');
my_saveas(gcf, fullfile(reed_solomon_dir, ['reed_solomon_pools_overlap_q_' num2str(q) '_deg_' num2str(deg) '_poly_'  dec2bin(p)]), {'epsc', 'jpg', 'pdf'});
people_overlap_mat = layered_code_mat' * layered_code_mat;
figure; imagesc(people_overlap_mat - diag(diag(people_overlap_mat))); colorbar; 
title('Overlap of people');
my_saveas(gcf, fullfile(reed_solomon_dir, ['reed_solomon_individuals_overlap_q_' num2str(q) '_deg_' num2str(deg)]), {'epsc', 'jpg', 'pdf'});
people_overlap_mat = layered_code_mat' * layered_code_mat;
figure; hist(people_overlap_mat(1,:), 0:4); xlabel('number of overlapping pools');
ylabel('number of pools');

random_code_mat = rand(size(layered_code_mat)) > 0.5; % Save random pooling matrices  
figure; colormap(bw); imagesc(random_code_mat); % colorbar;
my_saveas(gcf, fullfile(reed_solomon_dir, ['random_code_matrix']), ...
    {'epsc', 'jpg', 'pdf'});
random_sparse_code_mat = rand(size(layered_code_mat)) < 0.0625; % Save random pooling matrices  
figure; colormap(bw); imagesc(random_sparse_code_mat); % colorbar;
my_saveas(gcf, fullfile(reed_solomon_dir, ['random_sparse_code_matrix']), ...
    {'epsc', 'jpg', 'pdf'});

