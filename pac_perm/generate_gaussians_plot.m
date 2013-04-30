% Plot nicely all the gaussians
% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest
ROSEN = 9; HEPATO =10; % New data's from Assif
Bhattacharjee=11;  YEOH = 12;  % New Lung data
BRAIN = 13; KIM = 14; % New aging related datasets

RAND_DATA = 15; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)


load('all_datas_vars.mat');

data_flag = OLD_VANT_VEER; % NEW_ROSSETA; % NEW_ROSSETA;  % Choose which data to plot

indexes_vec = [8 17 26 35];  % indexes_vec = [5 12 19 28 35];  % This is good for new Rosseta
indexes_vec = [3 6 9 12];   % This is good for OLD Rosseta

color_vec = 'brgmcky:o*+sdxv<>^';

figure; subplot(1,2,1); hold on;
%title('Rosseta', 'fontsize', 14); 

for color_ind = 1:length(indexes_vec)
    errorbar(data_output_strct{data_flag}.N_f_vec(indexes_vec(color_ind)), data_output_strct{data_flag}.f_vec_saddle_mean(indexes_vec(color_ind)), ...
    data_output_strct{data_flag}.f_vec_saddle_std(indexes_vec(color_ind)), 'LineWidth',2, 'Color', color_vec(color_ind) );
end

plot(data_output_strct{data_flag}.N_f_vec, data_output_strct{data_flag}.f_vec_saddle_mean, 'LineWidth',2, 'Color', 'k');

%errorbar(data_output_strct{data_flag}.N_f_vec(indexes_vec), data_output_strct{data_flag}.f_vec_model_sim_mean(indexes_vec), data_output_strct{data_flag}.f_vec_model_sim_std(indexes_vec), 'r' );
xlabel('n', 'fontsize', 14); % ylabel('f', 'fontsize', 14); %legend('Saddle point solution', 'simulations');
text(-3,0.108,'f_n^*', 'fontsize', 14); % This is the ylabel



y_min = 0; y_max = 0.105; x_min=0; x_max = data_output_strct{data_flag}.N_f_vec(end) * 1.05; AXIS([x_min x_max y_min y_max]); 

get(gca);
set(gca, 'FontSize',14);

set(gca, 'YTick', []); % Eliminate numbers on the Y axis



% Here plot different Gaussians
subplot(1,2,2); hold on;  xlabel('P_{n, \alpha}(f)', 'fontsize', 14); % ylabel('f', 'fontsize', 14);

f_res = 0.001;
f_vec = [f_res:f_res:1-f_res];

legend_str = {};

for color_ind = 1:length(indexes_vec) % length(data_output_strct{data_flag}.N_f_vec)
    %            plot(f_vec, (sqrt(Ngenes)/(sqrt(2*pi) * inf_limit_std(index,color_ind)  )), color_vec(color_ind));
    plot( (1/(sqrt(2*pi) * data_output_strct{data_flag}.f_vec_saddle_std(indexes_vec(color_ind))  )) .* ...
        exp ( -(f_vec -  data_output_strct{data_flag}.f_vec_saddle_mean(indexes_vec(color_ind))).^2 .* ...
        1 ./(2.0 *data_output_strct{data_flag}.f_vec_saddle_std(indexes_vec(color_ind)).^2)) , f_vec, 'LineWidth',2, 'Color', color_vec(color_ind));
    
        legend_str{color_ind} =  [num2str((data_output_strct{data_flag}.N_f_vec(indexes_vec(color_ind))))];
end




text(-3,0.108,'f', 'fontsize', 14); % This is the ylabel
legend(legend_str);

% % % legend([num2str(data_output_strct{data_flag}.N_f_vec(1)) ' samples'], [num2str(data_output_strct{data_flag}.N_f_vec(2)) ' samples'], ...
% % %     [num2str(data_output_strct{data_flag}.N_f_vec(3)) ' samples'], ...
% % %     [num2str(data_output_strct{data_flag}.N_f_vec(4)) ' samples'], [num2str(data_output_strct{data_flag}.N_f_vec(5)) ' samples']);
%title('Prob. density func. of f for different number of samples'); xlabel('f'); ylabel('Pr(f)');

y_min = 0; y_max = 0.105; x_min=0; x_max = 50; AXIS([x_min x_max y_min y_max]); 
get(gca);
set(gca, 'FontSize',14);


