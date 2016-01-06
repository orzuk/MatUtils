% Run different values of N for Crohn's
crohns_mu = 0.001; % 0.001;
crohns_lambda_mz =   54.2592; % 54.3; % 250; % 54.3;
crohns_lambda_s =  10.2633; % 10.3; % 35; % 10.3;
AssignGeneralConstants;


for N=1:10 % fit different values of N
    run_N = N
    crohns_mu_l = 1-(1-crohns_mu)^(1/N)
    [h_pathway(N) c_R(N)] = familial_risk_to_LP_parameters( ...
        [crohns_lambda_mz crohns_lambda_s], N,  [1 0.5], 'binary', crohns_mu); % fit model (heritability + shared environment)
    h_pathway_shared(N) = c_R(N) * (1-h_pathway(N)); % this is shared within a pathway
    [lambda_R{N} stat_struct{N}] = ...
        compute_k_of_N_liabilities_statistics(N, 1, crohns_mu_l, ...
        h_pathway(N), [h_pathway_shared(N) h_pathway_shared(N) 0 0 0 0], 6); % allow shared environment only for siblings
    [lambda_R_no_shared_env{N}] = ...
        compute_k_of_N_liabilities_statistics(N, 1, crohns_mu_l, ...
        h_pathway(N), 0, 6); % allow shared environment only for siblings
    %     [lambda_R_full_shared_env{N}] = ...
    %         compute_k_of_N_liabilities_statistics(N, 1, crohns_mu_l, ...
    %         h_pathway(N), h_pathway_shared(N), 6); % allow shared environment only for siblings
    
    ctr_x=1;
    for x_shared = [0 0.25 0.5 0.75 1] % loop on fraction of shared environment fr distant relatives
        [lambda_R_distant_shared_env{N,ctr_x}] = ...
            compute_k_of_N_liabilities_statistics(N, 1, crohns_mu_l, ...
            h_pathway(N), h_pathway_shared(N) .* [1 1 x_shared x_shared x_shared x_shared], 6); % allow shared environment only for siblings
        ctr_x=ctr_x+1;
    end
    
    r_MZ(N) = familial_risk_to_heritability(lambda_R{N}(1), 'liability', crohns_mu, 1);
    r_MZ_no_shared_env(N) = familial_risk_to_heritability(lambda_R_no_shared_env{N}(1), 'liability', crohns_mu, 1);
    h_shared(N) = r_MZ(N) - r_MZ_no_shared_env(N); % compute shared environment on liability scale
    
    h_all(N) = stat_struct{N}.h_liab_loci;
    h_pop(N) = stat_struct{N}.h_liab_pop;
    pi_phantom(N) = stat_struct{N}.pi_liab_phantom;
    if(h_shared(N) < 0.0000001)
        h_shared(N) = 0;
    end
    if(h_pathway_shared(N) < 0.0000001)
        h_pathway_shared(N) = 0;
    end
    if(c_R(N) < 0.0000001)
        c_R(N) = 0;
    end
end


all_lambda_R = cell2vec(lambda_R);
% all_lambda_R_half_shared_env = cell2vec(lambda_R_half_shared_env);
% all_lambda_R_full_shared_env = cell2vec(lambda_R_full_shared_env);
for ctr_x=1:5
    all_lambda_R_distant_shared_env{ctr_x} = cell2vec(lambda_R_distant_shared_env(:,ctr_x));
end

all_lambda_R(1,:) = crohns_lambda_mz; all_lambda_R(2,:) = crohns_lambda_s;

R_header = {'N', 'h_{pathway}^2', 'c_R',  '\lambda_MZ', '\lambda_s', ...
    '\lambda_grand (0, 1/4, 1/2, 3/4, full shared-env)', 'lambda_cousin (0, 1/4, 1/2, 3/4, full shared-env)', ...
    'h_all^2', 'h_pop^2', 'V_c', '\pi_{phantom}'};
pi_phantom(1)=0;

R = [(1:N)' h_pathway' c_R'  all_lambda_R(1:2,:)']; % ... % just sibs
for ctr_x=1:5
    R = [R all_lambda_R_distant_shared_env{ctr_x}(3:4,:)'];
end
%all_lambda_R_half_shared_env(3:4,:)' all_lambda_R_full_shared_env(3:4,:)' ...
R = [R h_all' h_pop' h_shared' pi_phantom']; % add other information
R = num2str_cell(num2cell(R,3), 3);


for i=1:size(R,1) % concatenate the five different ones
    R{i,6} = ['(' R{i,6} ', ' R{i,8} ', ' R{i,10} ', ' R{i,12} ', ' R{i,14} ')'];
    R{i,7} = ['(' R{i,7} ', ' R{i,9} ', ' R{i,11}  ', ' R{i,13}  ', ' R{i,15} ')'];
end
R = R(:, [1:7 16:end]);
R = [R_header' R']';
R = num2str_cell(R, 3);

IBD_dir = '../../common_disease_model/data/visscher';
crohns_outfile = fullfile(IBD_dir, 'crohns_new_many_N_with_quarters.txt');
figs_dir = fullfile(IBD_dir, 'figs');

savecellfile(R, crohns_outfile, [], 1);


plot_lambda = 0; % skip printing for now
if(plot_lambda) % Plot decrease in lambda_mz as function of familial relationship
    k_shared_vec = [0.125 0.25 0.5 1];
    for i=1:3
        switch i
            case 1
                cur_lambda_R = all_lambda_R; shared_str = 'NONE';
            case 2
                cur_lambda_R = all_lambda_R_half_shared_env;  shared_str = 'HALF';
            case 3
                cur_lambda_R = all_lambda_R_full_shared_env;  shared_str = 'ALL';
        end
        figure; hold on;
        for N=1:10
            plot(k_shared_vec, log2(cur_lambda_R(end-2:-1:1,N)), [color_vec(mod_max(N,5)) symbol_vec{ceil(N/5)}]);
        end
        xlabel('Kinship'); ylabel('log(\lambda_R)');
        title(['distant-relatives have ' shared_str  ' of sibling shared environment']);
        legend(num2str((1:10)'), 2);
        for N=1:10
            plot(k_shared_vec, log2(cur_lambda_R(end-2:-1:1,N)), ['.' color_vec(mod_max(N,5)) symbol_vec{ceil(N/5)}], ...
                'markersize', 12);
        end
        text(1, 0.1, 'MZ', 'fontweight', 'bold');
        text(0.5, 0.1, 'DZ', 'fontweight', 'bold');
        text(0.25, 0.1, 'grand', 'fontweight', 'bold');
        text(0.125, 0.1, 'cousin', 'fontweight', 'bold');
        my_saveas(gcf, fullfile(figs_dir, ['crohns_lambda_R_' shared_str '_shared_environment']), ...
            {'epsc', 'pdf', 'fig', 'jpg'});
    end
end % if plot lambda

