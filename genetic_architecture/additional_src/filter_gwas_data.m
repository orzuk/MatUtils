% Filter association with data according to desired format.
% Determine if trait is binary (qualitative, with OR) or QTL (with beta).
% Also split diseases by severity (from Eliane)
%
% Input:
% data - structure with various fields with data from gwas results
% output_data_file - where to save data with added fields completed
% catalog_type - from NHGRI of from MIT turks (default)
%
% Output:
% filtered_data - data after filtering and throwing away bad and quantitative traits
% filtered_inds - indices of retained studies
%
function [filtered_data filtered_inds] = ...
    filter_gwas_data(data, output_data_file, catalog_type)

DISCOVERY = 1; REPLICATION = 2; COMBINED = 3;

MAX_ALLOWED_OR = 25; % maximum OR/GRR which makes sense for a binary locus
MAX_ALLOWED_BETA = 1; % maximum beta (in units of st.d.) which make sense for a QTL locus

quantitative_traits = {'AIDS progression', 'Adiponectin levels', 'Adiposity', ...
    'Angiotensin-converting enzyme activity', 'Biomedical quantitative traits', ...
    'Blood pressure', 'Body mass index', 'Bone mineral density', ...
    'Electrocardiographic conduction measures', 'Electrocardiographic traits', ...
    'Fibrinogen', 'Folate pathway vitamin levels', ...
    'HDL cholesterol',  'Height', ...
    'Hip bone size',  'LDL cholesterol', ...
    'Menarche (age at onset)', 'Menarche and menopause (age at onset)', ...
    'PR interval',  'Plasma E-selectin levels', 'Plasma eosinophil count', ...
    'QT interval', 'QT interval prolongation',   'Systolic blood pressure', ...
    'Triglycerides', 'Warfarin maintenance dose',  'Weight'};  % list of quantitative traits in table
neutral_traits = {'Angiotensin-converting enzyme activity', 'Electrocardiographic traits', 'PR interval', ...
    'Osteoporosis', 'Osteoarthritis', 'Serum metabolites', 'Pulmonary function', ...
    'Adiponectin levels', 'Height', 'Hair morphology', 'Serum iron concentration', ...
    'Weight', 'Response to citalopram treatment', 'Response to Hepatitis C treatment', ...
    'Menarche and menopause (age at onset)', 'Menarche (age at onset)', 'Systolic blood pressure', ...
    'Biomedical quantitative traits', 'LDL cholesterol', 'Fibrinogen', 'Arterial stiffness', ...
    'QT interval', 'Warfarin maintenance dose', 'Folate pathway vitamin levels', 'Bone mineral density', ...
    'Electrocardiographic conduction measures', 'Plasma eosinophil count', 'HDL cholesterol', ...
    'Response to treatment for acute lymphoblastic leukemia', 'Blood pressure', ...
    'Body mass index', 'Male-pattern baldness', 'Hip bone size', 'Serum IgE levels', ...
    'Response to statin therapy', 'Soluble ICAM-1', 'Blond vs. brown hair color', 'Blue vs. green eyes', ...
    'Burning and freckling','Freckles','Red vs. non-red hair color','Triglycerides', ...
    'Blue vs. brown eyes', 'Skin pigmentation', 'F-cell distribution'}; % Stratify traits by how bad is it to get them
lethal_traits = {'Multiple sclerosis', 'Inflammatory bowel disease (early onset)', ...
    'Myopia (pathological)', 'Amyotrophic lateral sclerosis' ... % 'Parkinson''s disease',
    'Acute lymphoblastic leukemia (childhood)', 'Glioma (high-grade)', 'Narcolepsy' ...
    'Myocardial infarction (early onset)', 'Hirschsprung''s disease', 'Kawasaki disease', ...
    'Asthma (childhood onset)', 'Celiac disease'}; % serious traits happening early in life
other_traits = {};
%    'Drug-induced liver injury (flucloxacillin)', ... % this is binary but mult. model is problematic


filtered_data = data; num_entries = length(filtered_data.SNPs);
filtered_data.Trait = empty_cell_to_empty_str(filtered_data.Trait);
filtered_data.RAF = cell2mat(empty_cell_to_numeric_val(str2nums_cell(data.Risk_Allele_Frequency,1,1), -1));

filtered_data.OR = zeros(num_entries,3); % New: Set, (discovery, replication, combined)
filtered_data.Power_OR = zeros(num_entries,1); % another effect size used for power correction (currently onley for breast cancer)
filtered_data.OR_interval = zeros(num_entries,2);
for i=1:num_entries
    tmp_discovery_nums = str2nums(data.OR_or_beta_discovery{i});
    tmp_replication_nums = str2nums(data.Replication_Effect_Size{i});
    tmp_combined_nums = str2nums(data.Combined_Effect_Size{i});
    
    
    if(~isempty(tmp_discovery_nums))
        filtered_data.OR(i,1) = tmp_discovery_nums(1);
    end
    
    if(~isempty(tmp_replication_nums))
        filtered_data.OR(i,2) = tmp_replication_nums(1); % set replication effect size
        for j=2:min(3,length(tmp_replication_nums)) % fill also confidence intervals
            filtered_data.OR_interval(i,j-1) = tmp_replication_nums(j);
        end
    end
    if(~isempty(tmp_combined_nums)) % combined effect size
        filtered_data.OR(i,3) = tmp_combined_nums(1);
    end
end
%filtered_data.OR = cell2mat(empty_cell_to_numeric_val(str2num_cell(data.OR_or_beta), -1));
%%new_quant_inds = setdiff(1:length(filtered_data.OR), new_inds); % all others are assumed to be quantitative

[inter_traits quant_inds] = intersect_all(filtered_data.Trait, quantitative_traits);
switch catalog_type
    case 'NHGRI'
        new_inds = find( (filtered_data.OR(:,2) > 1) & (filtered_data.OR(:,2) < 10)); % filter by odds ratios
        new_quant_inds = find( (filtered_data.OR(:,2) > -1) & (filtered_data.OR(:,2) < 1)); % filter by quantitative effect. Assume 0<beta<1 (Problem: some give effect size less than one!!!!)
        new_quant_inds = intersect(new_quant_inds, find(abs(filtered_data.OR(:,2)) > epsilon)); % remove snps with no effect size
    case 'MIT-TURK'
        quant_inds = union(quant_inds, strmatch('QTL', filtered_data.trait_type)); % use existing annotation
        quant_inds = union(quant_inds, strmatch('Quantitative', filtered_data.trait_type));
        new_inds = find( (filtered_data.OR(:,2) > 0) & (filtered_data.OR(:,2) < MAX_ALLOWED_OR)); % filter by odds ratios. Large values probably not OR!!!
        new_quant_inds = quant_inds; % Trust all quant inds % find( (filtered_data.OR(:,2) > -1) & (filtered_data.OR(:,2) < 1)); % filter by quantitative effect. Assume 0<beta<1 (Problem: some give effect size less than one!!!!)
end

quant_inds = union(quant_inds, new_quant_inds);
binary_inds = setdiff(1:length(filtered_data.OR(:,2)), quant_inds);
filtered_data.trait_type_num = zeros(num_entries,1); % set default as binary
filtered_data.trait_type_num(quant_inds) = 1;
filtered_data.trait_type = cell(num_entries, 1);
[neutral_traits_unique neutral_inds] = intersect_all(filtered_data.Trait, neutral_traits); % Compute how bad is each trait
[lethal_traits_unique lethal_inds] = intersect_all(filtered_data.Trait, lethal_traits);

for i=1:num_entries %   num_filtered % loop over all binary traits
    if(filtered_data.trait_type_num(i))
        filtered_data.trait_type{i} = 'Quantitative';
    else
        filtered_data.trait_type{i} = 'Binary';
    end
    filtered_data.lethality{i} = 1; % defaut lethality intermidiate
end
for i=1:length(neutral_inds)
    filtered_data.lethality{neutral_inds(i)} = 0; % netrual traits
end
for i=1:length(lethal_inds)
    filtered_data.lethality{lethal_inds(i)} = 2; % really lethal traits
end
filtered_data.lethality = vec2column(cell2mat(filtered_data.lethality));


%good_inds1 = find(isnumeric_cell(data.Risk_Allele_Frequency) & (~isempty_cell(data.Risk_Allele_Frequency)));
%good_inds2 = find(isnumeric_cell(data.OR_or_beta)  & (~isempty_cell(data.OR_or_beta)));
bad_inds = setdiff(1:length(filtered_data.OR(:,2)), union(new_inds, new_quant_inds)); % if new_inds and new_quant_inds didn't pick then, these are bad/unknown and should be thrown away
good_inds1 = find(filtered_data.RAF > -1); % require a valid RAF in [0,1]
good_inds2 = find(filtered_data.OR(:,2) > -MAX_ALLOWED_BETA); % require a valid effect size can't have |beta|>1??)
good_inds = intersect(good_inds1, good_inds2);
good_inds = setdiff(good_inds, bad_inds); % retention criteria: require both valid RAF, OR, and appearing in new_inds or new_quant_inds

empty_snp_inds = find(isempty_cell(filtered_data.SNPs));
good_inds = setdiff(good_inds, empty_snp_inds);
filtered_inds = good_inds; % alternative
empty_traits_inds = find(isempty_cell(data.Trait));
filtered_inds = setdiff(filtered_inds, empty_traits_inds); % remove empty traits 

%% filtered_inds = good_inds(new_inds); % Here remove quantitative ones:
num_filtered = length(filtered_inds);



% New: also save original data as tab-delimited .txt file, with lethality and quant fields
data.trait_type = filtered_data.trait_type;
data.lethality = filtered_data.lethality;
cell_to_mat=0;
R = WriteDataFile(data, output_data_file, cell_to_mat, 0); % why is this so slow? 

filtered_data = struct_by_inds(filtered_data, filtered_inds); % rewrite filtered data

%filtered_inds = filtered_inds(good_inds);


% filtered_data.RAF = cell2mat(str2num_cell(filtered_data.Risk_Allele_Frequency));
% filtered_data.OR = cell2mat(str2num_cell(filtered_data.OR_or_beta));
filtered_data.Trait = empty_cell_to_empty_str(filtered_data.Trait);
unique_traits = unique(filtered_data.Trait); % Find manually which ones are quantitative or qualitative traits



if(~isfield(filtered_data, 'VarianceofTrait')) % prepare variance vector 
    filtered_data.VarianceofTrait = ones(length(filtered_data.trait_type_num), 1); 
    bmi_inds = strmatch('Body', filtered_data.Trait); 
    filtered_data.VarianceofTrait(bmi_inds) = 4.5; % empirical BMI st.d. 
end


% Convert continuous traits effect size to binary traits effect sizes. Also
% flip effect size and RAF
filtered_data.Beta = zeros(length(filtered_data.trait_type_num), 1);
for i=1:length(filtered_data.trait_type_num)
    
    
    
    if(filtered_data.trait_type_num(i)) % convert continuous QTL to binary
        switch filtered_data.Effect_Size_Units{i}    % First convert values to true units 
            case {'beta', 'Beta', 'Beta (95% CI)', 'Beta (Standard error)', 'Beta (standard error)'}
                % do nothing 
            case {'%SD', '% SD', 'Percentage of standard deviation'}
                filtered_data.OR(i,:) = filtered_data.OR(i,:) ./ 100; % transfer from %                 
            case {'cm', 'cm/s increase', 'kg/m^2', 'ml/min/1.73 m2 increase', ...
                    'per allele change in age (weeks) at menarche', 'ug/ml', 'unit increase'}
                filtered_data.OR(i,:) = filtered_data.OR(i,:) ./ filtered_data.VarianceofTrait(i); % normalize by trait st.d.
        end        
        [filtered_data.RAF(i), filtered_data.OR(i,2) flipped_flag] = ...  % negative beta - flip risk allele !
            flip_allele(filtered_data.RAF(i), filtered_data.OR(i,2), 'QTL', 2);
        if(flipped_flag) % If flipped, need to flip also two others
            for j=[1 3]
                if(filtered_data.OR(i,j) ~= 0)
                    [~, filtered_data.OR(i,j)] = ...  % replication (?) GRR<1 - flip risk allele !
                        flip_allele(1-filtered_data.RAF(i), filtered_data.OR(i,j), 'QTL');
                end
            end
        end
        for j=1:3
            filtered_data.Beta(i,j) = filtered_data.OR(i,j); % copy a continuous value
            if(filtered_data.Beta(i,j) < MAX_ALLOWED_BETA) % make sure we got st.d. units
                filtered_data.OR(i,j) = ... % convert the already flipped (high) OR
                    beta_to_genetic_relative_risk(filtered_data.Beta(i,j), ...
                    filtered_data.RAF(i), 0.01); % we don't know the prevalence!!! put here 0.01!!!
            else
                filtered_data.OR(i,j) = 1; % here we don't know the units!! set effect size at zero!!!
            end
            if(filtered_data.OR(i,j) < 1) % negative effect sizes (from beta < 0).
                filtered_data.OR(i,j) = 1; % again we don't know units - set effect size at zero
                why_is_effect_size_negative_beta_too_high = 888
                %                 iii_is = i
                %                     error(888)
            end
            if(~isreal(filtered_data.OR(i,j)))
                effect_size_please_get_real = 999
                iii_is = i
                error(9999)
            end
        end % loop on 3 possible effect sizes
        
    else % here correct BINARY traits. Correct both RAF and MAF!!!
        [filtered_data.RAF(i), filtered_data.OR(i,2) flipped_flag] = ...  % replication (?) GRR<1 - flip risk allele !
            flip_allele(filtered_data.RAF(i), filtered_data.OR(i,2), 'binary', 2);
        if(flipped_flag) % If flipped, need to flip also two others
            for j=[1 3]
                if(filtered_data.OR(i,j) ~= 0)
                    [~, filtered_data.OR(i,j)] = ...  % replication (?) GRR<1 - flip risk allele !
                        flip_allele(1-filtered_data.RAF(i), filtered_data.OR(i,j), 'binary');
                end
            end
        end
    end
end % loop on all retained snps

%good_inds = setdiff(1:length(filtered_data.OR), bad_inds); % Throw away bad inds
%filtered_data = struct_by_inds(filtered_data, good_inds);

filtered_data.MAF = min(filtered_data.RAF, 1-filtered_data.RAF); % set minor-allele-frequency
flipped_inds = find(filtered_data.RAF > 0.5);
flipped_binary_inds = intersect(flipped_inds, binary_inds);
flipped_quant_inds = intersect(flipped_inds, quant_inds);
filtered_data.OR_with_MAF = filtered_data.OR;
filtered_data.OR_with_MAF(flipped_inds,:) = 1 ./ filtered_data.OR_with_MAF(flipped_inds,:);
%filtered_data.OR_with_MAF(flipped_quant_inds,:) = -filtered_data.OR_with_MAF(flipped_quant_inds,:);

filtered_data.Beta_with_MAF = filtered_data.Beta; % we've got beta only for quantinds
filtered_data.Beta_with_MAF(flipped_quant_inds,:) = -filtered_data.Beta_with_MAF(flipped_quant_inds,:);
%filtered_data.Beta_with_MAF(flipped_binary_inds,:) = -filtered_data.Beta_with_MAF(flipped_binary_inds,:);



