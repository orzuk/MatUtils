% Generate an excel table summarizing all parameters for a two-class model
%
% Input:
% two_class_stat_struct, ...
%    mutation_fraction_birth_vec, f_rare_vec, max_f_rare_vec, c_cumulative, ...
%    frac_null_by_freq_cumulative, ...
%    strong_prop, mutation_str, mutation_type_str, mutation_rate, ...
%    mu, theta, L, N, num_null_s, show_s_null, show_s_null_ind
%
% Output:
% R - cell-array with all values
%
function R = compute_two_class_table(two_class_stat_struct, ...
    mutation_fraction_birth_vec, f_rare_vec, max_f_rare_vec, c_cumulative, ...
    frac_null_by_freq_cumulative, ...
    strong_prop, mutation_str, mutation_type_str, mutation_rate, ...
    mu, theta, L, N, num_null_s, show_s_null, show_s_null_ind)


mutation_fraction_birth_vec = [mutation_fraction_birth_vec(1) ...
    mutation_fraction_birth_vec(1) .* strong_prop mutation_fraction_birth_vec(1) * (1-strong_prop) ...
    mutation_fraction_birth_vec(2:end)]; % comptue proportion of strong and neutral missense

R = [];
R_header = []; R_header{1,1} = 'Parameters:'; ctr=2;
R_header{ctr,1} = 'Gene-Length'; R_header{ctr,2} = num2str(L); ctr=ctr+1;
R_header{ctr,1} = 'Mut. Rate (per-nucleotide-per-generation)'; R_header{ctr,2} = ...
    num2str_log_pretty(mutation_rate); ctr=ctr+1;
R_header{ctr,2} = 'Prop.'; R_header{ctr,3} = 'Mut. Rate'; R_header{ctr,4} = 'Theta (=4N*mu)';  R_header{ctr,5} = 'Theta*L (=4N*mu*L)';ctr=ctr+1;
for j=1:length(mutation_str)
    R_header{ctr,1} = ['pi_' mutation_str{j}];
    R_header{ctr,2} = num2str(mutation_fraction_birth_vec(j),2);
    R_header{ctr,3} = num2str_log_pretty(mutation_rate*mutation_fraction_birth_vec(j));
    R_header{ctr,4} = num2str_log_pretty(4*N*mutation_rate*mutation_fraction_birth_vec(j));
    R_header{ctr,5} = num2str_log_pretty(L*4*N*mutation_rate*mutation_fraction_birth_vec(j));
    ctr=ctr+1;
end
R_header{ctr,1} = []; R_header{ctr,num_null_s+1} = []; % ctr=ctr+1; % add an empty line

R_header = [R_header' ['Selection:' '0' num2str_cell(num2cell(log10(show_s_null(2:end))))]']';
for j=3:num_null_s+1
    R_header{end,j} = ['10^(' R_header{end,j} ')'];
end

%             R_header{ctr,1} = 'Mut. Rate'; R_header{ctr,2} = num2str(mutation_rate); ctr=ctr+1;
%             R_header{ctr,1} = 'Mut. Rate'; R_header{ctr,2} = num2str(mutation_rate); ctr=ctr+1;

R = [R R_header];


freq_headers = {'Ave. Freq. of All Alleles', 'Median Freq. of All Alleles', ...
    'Ave. Freq. of Null Alleles', 'Median Freq. of Null Alleles'};
freq_table = zeros(num_null_s, 4);

freq_table(:,1) = two_class_stat_struct.mean_mixture_x(show_s_null_ind);
freq_table(:,2) = two_class_stat_struct.median_mixture_x(show_s_null_ind);
freq_table(:,3) = two_class_stat_struct.mean_x(show_s_null_ind);
freq_table(:,4) = two_class_stat_struct.median_x(show_s_null_ind);


%            R{end,num_null_s+1} = [];
R = [R' [freq_headers' num2cell(freq_table')]']'; % add freq information
R{end+1,1} = [];

for j=1:length(max_f_rare_vec)
    R_cumulative = {['Freq. of Alleles (f<' num2str(max_f_rare_vec(j)) ')'], 'Proportion of Null Alleles', ...
        'Null Alleles Under-Representation', ...
        'Effect-Size Correction Factor Quant. Traits', 'Effect-Size Correction Factor Disease Traits (approx.)'}';
    [~, f_index] = min((f_rare_vec - max_f_rare_vec(j)).^2);
    cum_table = zeros(num_null_s, 5);
    cum_table(:,1) = c_cumulative{1}(show_s_null_ind,f_index);
    cum_table(:,2) = frac_null_by_freq_cumulative{1}(show_s_null_ind,f_index);
    cum_table(:,3) = two_class_stat_struct.bayes_factor_vec(show_s_null_ind,f_index);
    cum_table(:,4) = frac_null_by_freq_cumulative{1}(show_s_null_ind,f_index); % correction factor is same as null proportion for null alleles
    cum_table(:,5) = frac_null_by_freq_cumulative{1}(show_s_null_ind,f_index); % correction factor is same as null proportion for null alleles
    
    
    R_cumulative = [R_cumulative num2cell(cum_table)']';
    R_cumulative{1,end+1} = [];
    R = [R' R_cumulative]';
end

R_expected_prop_carriers = [];
R_expected_prop_carriers{1,1} = 'Expected # alleles (missense+stop+frameshift) per-person';
R_expected_prop_carriers{2,1} = 'Expected # substituation alleles (missense+stop) per-person'; ctr=3;

for j=1:length(mutation_str)
    R_expected_prop_carriers{ctr,1} = ['Expected # ' mutation_str{j} ' Alleles per-person ' mutation_type_str{j}]; ctr=ctr+1;
end
R_expected_prop_carriers{end,1} = 'Expected # total-substituation alleles (missense+stop+synonymous) per-person';
prop_carriers_table = zeros(num_null_s, 9);

for j=1:6 % loop on strong/neutral missense, and then frameshift, stop
    prop_carriers_table(:,j+3) = theta * L * mutation_fraction_birth_vec(j+1);
    
    switch mutation_str{j+1}
        case {'missense-strong', 'stop', 'frameshift'} % all null mutations
            %   if(j ~= 2) % this is for all allele types except the neutral missense
            prop_carriers_table(:,j+3) = prop_carriers_table(:,j+3) .* ...
                two_class_stat_struct.bayes_factor_vec(show_s_null_ind,end); %
    end
end
prop_carriers_table(:,9) = sum(prop_carriers_table(:,[4 5 6 8]), 2); % all substitutions
prop_carriers_table(:,3) = sum(prop_carriers_table(:,[4 5]),2); % all carriers
prop_carriers_table(:,1) = sum(prop_carriers_table(:,[3 6 7]),2); % all missense carriers
prop_carriers_table(:,2) = sum(prop_carriers_table(:,[3 6]),2); % all carriers - no frameshift

R_expected_prop_carriers = [R_expected_prop_carriers num2cell(prop_carriers_table')]';
R_expected_prop_carriers = R_expected_prop_carriers(:,[3:end 1:2]);

R = [R' R_expected_prop_carriers]';
R{end+1,1} = [];

R_var = {'Time Null Allele Spends as Polmorphic (# generations)', 'Mean Heterozygosity of Null Alleles', ...
    'Proprtion of Var. Explained by Null for Quant. Traits (\beta=1)'};

grr_minus_one_vec = [2-1 5-1];
mu_vec = [0.01 0.1]; % different mu's per Eric's request
var_table = zeros(num_null_s, 3+2*length(grr_minus_one_vec)*length(mu_vec));
var_table(:,1) = phi_s_integral(1-1/(2*N), -4*N*show_s_null, 0) - phi_s_integral(1/(2*N), -4*N*show_s_null, 0); % heterozygosity
var_table(:,2) = phi_s_integral(0.999999999, -4*N*show_s_null, -1) - phi_s_integral(0.000000001, -4*N*show_s_null, -1); % heterozygosity
var_table(:,3) = theta * L * ...  % quantitative traits
    (phi_s_integral(0.999999999, -4*N*show_s_null, -1) - phi_s_integral(0.000000001, -4*N*show_s_null, -1));
for j=1:length(grr_minus_one_vec)
    for k=1:length(mu_vec)
        R_var{(k-1)*4+2*j+2} = ['Proprtion of Var. Explained by Null for Disease Traits (K=' num2str(mu_vec(k)) ...
            ',GRR=' num2str(grr_minus_one_vec(j)+1) ')'];
        R_var{(k-1)*4+2*j+3} = ['Proprtion of Var. Explained by Null for Disease Traits on liability scale (K=' num2str(mu_vec(k)) ...
            ',GRR=' num2str(grr_minus_one_vec(j)+1) ')'];
        var_table(:,(k-1)*4+2*j+2) = grr_minus_one_vec(j) .^ 2 .* theta * L * (mu_vec(k)/(1-mu_vec(k))) * ... % disease traits
            (phi_s_integral(0.999999999, -4*N*show_s_null, 'disease', grr_minus_one_vec(j)) - ...
            phi_s_integral(0.000000001, -4*N*show_s_null, 'disease', grr_minus_one_vec(j)));
        var_table(:,(k-1)*4+2*j+3) = heritability_scale_change(var_table(:,(k-1)*4+2*j+2), 'liability', mu_vec(k)); % heritability on liability scale
    end
end

R_var = [R_var' num2cell(var_table)']';
R = [R' R_var]';
R = num2str_cell(R, 3); % convert to strings

for k=1:size(R,2)
    max_field_len = max(length_cell(R(:,k)));
    for j=1:size(R,1)
        if(k == 1)
            R{j,k} = [R{j,k} repmat('.', 1, max_field_len - length(R{j,k}))];
        else
            R{j,k} = [repmat(' ', 1, max_field_len - length(R{j,k})) R{j,k}];
        end
    end
end
