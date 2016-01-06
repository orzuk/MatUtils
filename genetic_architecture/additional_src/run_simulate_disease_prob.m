% Script for running simulate_disease_prob for different models
AssignGeneralConstants;
AssignStatsConstants;

% Flags saying what do we want to run
do_additive = 0; % look at reverse simulations with additive model
do_pairs = 0; % look at pairs of loci
do_examples = 0; % look at examples
test_h_lods = 1; % small testing of 2x2 tables 

N = 500; % number of loci
f = ones(N,1) .* 0.1; % MAF of each locus
d = ones(N,1); % what is this?
a = 3; % Maximal probability of disease (can be less than one)
generations = 5; % how far to go on family tree

iters = 10000; % loop over different models
N_vec = [1 100 100 100]; % num of loci again (?)

d_vec = [10000 1 1000 1];
a_vec = [1 1 1 1]; % set sigmoid parameters
b_vec = [1 1 1 1]; % set sigmoid parameters
c_vec = [1 1 1 1]; % set sigmoid parameters
labels_vec = {'Single-Locus', 'Common-Weak', 'Rare-Strong', 'Protective'}; % different models (all linear)
num_models = length(N_vec);
rho = zeros(num_models,1); %
RARE = 0; COMMON = 1;

% if(do_pairs) % Simulate with pairwise interactions. Try to make it work for general N (not only pairs) but require sampling ..
%     run_compute_genetic_architectures();    
% end % if do pairs















%%%%%%%%%%% Old Stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(do_additive) % look at additive model of disease
    for disease = [COMMON]
        r = zeros(num_models,generations); h = zeros(num_models,1);
        if(disease == RARE)
            f_vec = [0.001 0.2 0.001 0.999];
            rho(1) = 0.001; h(1) = 0.1; % rare mandelian disease
            disease_str = 'rare';
        else
            f_vec = [0.0742 0.2 0.002 0.998]; % The frequency of mendelian model was increased
            rho(1) = 0.1; h(1) = 0.7; % rare mandelian disease
            disease_str = 'common';
        end
        v = cell(num_models,1); P=v; S=v;
        for i=1:num_models
            f = repmat(f_vec(i), N_vec(i), 1); % common alleles
            d = repmat(d_vec(i), N_vec(i), 1); % doesn't really matter ..
            [v{i} d_norm] = normalize_allele_risk(f, d);
            [S{i} SW{i}] = simulate_disease_risk(f, d, iters, generations); % Get the disease risk (heavy part)
            %    [P{i} S{i} rho(i) h(i) r(i,:)] = simulate_disease_prob(f, d, a_vec(i), 1, 1, iters, generations);
            [xxx b_vec(i) c_vec(i)] = fit_risk_sigmoid(S{i}, rho(1), h(1)); % solve for a given rho
            [P{i} rho(i) h(i) r(i,:)] = ...
                simulate_disease_prob(f, d, 1, b_vec(i), c_vec(i), generations, S{i}, SW{i}); % run again with the fitted value
            %     figure; hist(P{i}, 50); title(['P ' labels_vec{i} ' dist.']);
            %     figure; hist(S{i}, 50); title(['Risk (S) ' labels_vec{i} ' dist.']);
        end
        r = [h r]; % concatenate with twins
        figure; hold on; % figure of risk decay with relatedness
        for i=1:num_models
            plot(0:generations, r(i,:), color_vec(i));
        end
        plot(0:generations, repmat(rho(1), generations+1,1), 'k--');
        set(gca, 'XTick', 0:generations); % set the xticks
        xlabel('Relative Degree'); ylabel('Prob. Disease'); legend([labels_vec 'population freq.']);
        for i=1:num_models
            plot(0:generations, r(i,:), [color_vec(i) '*']);
        end
        my_saveas(gcf, ['../../common_disease_model/figures/' disease_str '_disease_prob_relatives'], {'jpg', 'epsc', 'fig'});
        
        figure; hold on; % Figure showing the risk and sigmoid cut
        for i=1:num_models
            subplot(2,2,i); hold on; title(labels_vec{i});
            [h_ bins_vec] = hist_density(S{i}, 100,[],2); xlabel('S (Risk)'); ylabel('Freq.');
            plot(bins_vec, 1 ./ (a_vec(i)+exp(-bins_vec)), 'r'); % plot sigmoid
        end
        
        
        figure; hold on; % Figure shoiwing just how much is explained by each locus (when all loci the same its meaningless)
        for i=1:num_models
            plot(cumsum(sort(v{i} ./ sum(v{i}), 'descend')), color_vec(i));
        end
        legend(labels_vec); xlabel('loci'); ylabel('Frac. Var. Explained');
        my_saveas(gcf, ['../../common_disease_model/figures/' disease_str '_fraction_explained_by_loci'], {'jpg', 'epsc', 'fig'});
    end % loop on diseae type (RARE/COMMON)
end


if(do_examples) % What are these examples?
    [P S rho h r] = simulate_disease_prob(f, d, a, 1,1, generations);
    figure; hist(P, 50); title('P dist.');
    figure; hist(S, 50); title('Risk (S) dist.');
    
    N = 1; % Example with one rare allele: here we should have r = h/2:
    f = 0.001;
    d = 10000;
    iters = 1000 / f; % enough iteations to see f
    [P S rho h r] = simulate_disease_prob(f, d, a, 1, 1, generations);
    h_is = h
    r_is = r
    diff_is = h/2 - r
    should_be_the_same =  r .* (2.^(1:generations)')
    figure; hist(P, 50); title('P Mendelian dist.');
    figure; hist(S, 50); title('Risk (S) Mendelian dist.');
    
    N = 1200; % Try to loci and see if it's changed ...
    f2 = repmat(f, N, 1);
    d2 = repmat(d/N, N, 1);
    [P2 S2 rho2 h2 r2] = simulate_disease_prob(f2, d2, a, 1, 1, generations);
    h2_is = h2
    r2_is = r2
    diff2_is = h2/2 - r2
    
    N = 100; % Common alleles with small effect
    f_common = repmat(0.2, N, 1); % common alleles
    d_common = repmat(1, N, 1); % doesn't really matter ..
    [v_common d_common] = normalize_allele_risk(f_common, d_common);
    [P_common S_common rho_common h_common r_common] = ...
        simulate_disease_prob(f_common, d_common, a, 1, 1, generations);
    figure; hist(P_common, 50); title('P Common dist.');
    figure; hist(S_common, 50); title('Risk (S) Common dist.');
    
    N = 100; % Rare alleles with large effect model:
    f_rare = repmat(0.001, N, 1); % common alleles
    d_rare = repmat(1, N, 1); % doesn't really matter ..
    [v_rare d_rare] = normalize_allele_risk(f_rare, d_rare);
    [P_rare S_rare rho_rare h_rare r_rare] = simulate_disease_prob(f_rare, d_rare, a, generations);
    figure; hist(P_rare, 50); title('P Rare dist.');
    figure; hist(S_rare, 50); title('Risk (S) Rare dist.');
end

p=0.1;
for i=1:9
    p = p*0.9 + (1-p)*0.1;
end






if(test_h_lods)
    iters = 10000; 
    for maf = [0.01 0.05 0.1 0.2 0.5]
        for disease_freq = [0.01 0.05 0.1 0.2 0.5]
            P = zeros(iters,4); % Generate random P with given marginals 
            P(:,4) = rand(iters,1) * min(maf, disease_freq);
            P(:,3) = disease_freq - P(:,4); 
            P(:,2) = maf - P(:,4);
            P(:,1) = 1-sum(P,2); 
            
            test_inheritance_and_lods(iters,P); 
            title(['inheritance vs. lods ratio MAF=' num2str(maf) ...
                ' disease-freq.=' num2str(disease_freq)]);
        end
    end
end            
    




% x = 15;
% switch x
%     case {14, 15, 16}
%         y = 1234124
%         switch x
%
%             case 13
%                 yy = 1234
%             case 15
%                 yyy = 1341241234
%             case 16
%                 ttt = 134
%         end
%     otherwise
%         zzz = 1234
% end
