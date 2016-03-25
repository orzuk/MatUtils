% function compute_crohns_effect_sizes
% Parsing supp. tab. 3 from the paper [Franke et al., Nature Genetics, 2010]
% A script for parsing the excel file with Crohn's diseases by
% sub-populations, and estimating both discovery and replication effect
% sizes

% Parse file
R = loadcellfile('../../common_disease_model/data/Daly_Crohs_paper/Franke_et_al_Crohns_71_loci_meta_analysis_ng.717_supp_tab3.txt');
labels = R(1:2,:); R = R(3:73,:);  % seperate labels and data

load('disease_data.mat');
crohn_ind = strmatch('Crohn''s disease', data_params.trait_name);
crohn_snp_inds = strmatch('Crohn''s disease', data.Trait);


snp_ids = R(:,1)
chr = R(:,2)
pos = R(:,3);
rep_het = R(:,4)
populations = unique(labels(1,:)); populations = populations(~isempty_cell(populations));
num_populations = length(populations); num_snps = 71;
RAF_mat = zeros(num_snps, num_populations);
OR_mat = zeros(num_snps, num_populations);

for i=1:length(populations)
    RAF_ind = strmatch(populations{i}, labels(1,:));
    RAF_mat(:,i) = cell2mat(empty_cell_to_numeric_val(str2num_cell(R(:,RAF_ind)), -1));
    OR_mat(:,i) = cell2mat(empty_cell_to_numeric_val(str2num_cell(R(:,RAF_ind+1)), -1));
end

bad_inds = union(find(RAF_mat == 0), find(OR_mat == 0)); % bad indices in database;
RAF_mat(bad_inds) = -1; OR_mat(bad_inds) = -1; 
[RAF_mat OR_mat] = flip_allele(RAF_mat, OR_mat, 'Binary', 1); % Flip to ensure MAF < 0.5



discovery_populations = {'Adolescent - GWAS', 'Germany - GWAS', 'Cedar_1 - GWAS', ...
    'Belgium - GWAS', 'NIDDKJ - GWAS', 'WTCCC - GWAS'};% sample size should add to: 6,333 15,056
num_samples_discovery = [1689 6197; ...
    479 1145; ...
    925 2882; ...
    537 913; ...
    956 982; ...
    1747 2937];
num_cases_discovery = num_samples_discovery(:,1);
num_controls_discovery = num_samples_discovery(:,2);

replication_populations = {'Australia - repl', 'Belgium - repl', 'France - repl', 'Germany - repl', ...
    'Israel - repl', 'Italy - repl', 'Netherlands - repl', ...
    'New Zealand - repl', 'Spain - repl', 'Sweden - repl', ...
    'United Kingdom - repl', 'Cedar - repl', 'NIDDK - repl'};
num_samples_replication = [ 1357 1923; ...
    1282 1682; ...
    0 0; ... % 414 414*2;  ... % trios are not counted !!!!
    3808 2747;  ...
    444 376;  ...
    921 899;  ...
    1101 269;  ...
    514 457;  ...
    325 987; ...
    724 992; ...
    3243 2431;  ...
    1172 501; ...
    803 762];
num_cases_replication = num_samples_replication(:,1);
num_controls_replication = num_samples_replication(:,2);

[~, discovery_inds, discovery_inds2] = intersect(populations, discovery_populations) % strfind_cell(populations, 'repl')
num_cases_discovery = num_cases_discovery(discovery_inds2);
num_controls_discovery = num_controls_discovery(discovery_inds2);
%discovery_inds = setdiff(1:num_populations, replication_inds);
[~, replication_inds, replication_inds2] = intersect(populations, replication_populations) % strfind_cell(populations, 'repl')
num_cases_replication = num_cases_replication(replication_inds2);
num_controls_replication = num_controls_replication(replication_inds2);

mu = 0.002; % Set Crohn's Prevalence
[discovery_f_vec discovery_grr_vec] = ...
    meta_analysis_effect_size(vec2row(num_cases_discovery), vec2row(num_controls_discovery), ...
    RAF_mat(:,discovery_inds), OR_mat(:,discovery_inds), mu, 'Binary');

replication_flag_vec = ones(num_snps,1);
for i=1:num_snps % for each snp take only studies in which it participates
    cur_rep_inds = find(RAF_mat(i,replication_inds) > -1);
    if(~isempty(cur_rep_inds))
        [replication_f_vec(i) replication_grr_vec(i)] = ...
            meta_analysis_effect_size(vec2row(num_cases_replication(cur_rep_inds)), vec2row(num_controls_replication(cur_rep_inds)), ...
            RAF_mat(i,replication_inds(cur_rep_inds)), OR_mat(i,replication_inds(cur_rep_inds)), mu, 'Binary');
    else % just copy the iscovery effect size
        replication_flag_vec(i) = 0; % no replication for this SNP!!!
        replication_f_vec(i) = discovery_f_vec(i);
        replication_grr_vec(i) = 1; % temp discovery_grr_vec(i);
    end
end
figure; plot(discovery_grr_vec, replication_grr_vec, '.');
xlabel('discovery'); ylabel('replication');

[inter_snps I J] = intersect(strrep_cell(data.SNPs(crohn_snp_inds), ' ', ''), snp_ids); 

figure; hold on; plot(data.GRR(crohn_snp_inds(I)),  ...
    max(replication_grr_vec(J), 1./replication_grr_vec(J)), '.');
plot(1:2, 1:2, 'r'); 
xlim([1 1.35]);
xlabel('replication (published)'); ylabel('replication (computed)');

figure; hold on; plot(data.GRR(crohn_snp_inds(I)),  ...
    max(discovery_grr_vec(J), 1./discovery_grr_vec(J)), '.');
plot(1:2, 1:2, 'r'); 
xlim([1 1.35]);
xlabel('replication (published)'); ylabel('discovery (computed)');
