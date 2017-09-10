%function load_and_norm_data(sample_names, array_names, working_path, diploid_ind, chip)
function load_and_norm_data(sample_names, array_names, working_path, diploid_ind, chip)

diploid_ind = sort(diploid_ind);
num_pm_probes = 20;
num_samples = length(sample_names);
cdf_f_name = get_cdf_f_name(chip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first load all cell files and save their mean intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
median_pm_intens_vec = zeros(1, num_samples);
for i = 1:num_samples
    array_name = char(array_names{i});
    if(~is_file_pres(array_name))
        array_name = [working_path array_name];
    end
    [cell_path, cell_name] = full_f_name_into_dir_f_name(array_name);
    [cdf_path, cdf_name] = full_f_name_into_dir_f_name(cdf_f_name);
    ProbeStructure = celintensityread(cell_name, cdf_name, 'CELPath', cell_path, 'CDFPath', cdf_path);
    median_pm_intens_vec(i) = median(ProbeStructure.PMIntensities);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose baseline array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_med = median(median_pm_intens_vec);
[baseline_dummy baseline_ind] = min( (median_pm_intens_vec - baseline_med).^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load baseline array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array_name = char(array_names{baseline_ind});
if(~is_file_pres(array_name))
    array_name = [working_path array_name];
end
[cell_path, cell_name] = full_f_name_into_dir_f_name(array_name);
[cdf_path, cdf_name] = full_f_name_into_dir_f_name(cdf_f_name);
baseline_ProbeStructure = celintensityread(cell_name, cdf_name, 'CELPath', cell_path, 'CDFPath', cdf_path, ...
    'PMOnly', 'false');

% find probe sets that have 20 probes
probe_set_start_ind = find(baseline_ProbeStructure.ProbeIndices==0);
probe_set_len_vec = diff(probe_set_start_ind);

first_probe_set_ind = min(find(probe_set_len_vec==num_pm_probes));
first_probe_ind = probe_set_start_ind(first_probe_set_ind);
num_probes = length(baseline_ProbeStructure.PMIntensities(first_probe_ind:end));
snp_id = baseline_ProbeStructure.ProbeSetIDs(first_probe_set_ind:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now normalize data with the baseline array and calc allele ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:num_samples
    array_name = char(array_names{i});
    if(~is_file_pres(array_name))
        array_name = [working_path array_name];
    end
    [cell_path, cell_name] = full_f_name_into_dir_f_name(array_name);
    [cdf_path, cdf_name] = full_f_name_into_dir_f_name(cdf_f_name);
    ProbeStructure = celintensityread(cell_name, cdf_name, 'CELPath', cell_path, 'CDFPath', cdf_path, ...
        'PMOnly', 'false');
    % calc allele ratio

    probe_intens_mat = [(reshape(ProbeStructure.PMIntensities(first_probe_ind:end), num_pm_probes, num_probes / num_pm_probes))' ...
        (reshape(ProbeStructure.MMIntensities(first_probe_ind:end), num_pm_probes, num_probes / num_pm_probes))'];
    pa_vec = [1:10]; ma_vec = [21:30]; pb_vec = [11:20]; mb_vec = [31:40];
    allele_ratio_vec  = calc_probe_allele_ratio(probe_intens_mat, 1, ...
        pa_vec, ma_vec, pb_vec, mb_vec);

    %%%%%%%%%%%%
    %    supress normalization now
    %%%%%%%%%%%%
    supress_norm_flag = 1;
    if(~supress_norm_flag)
        [NormData_probes MedStructure] = affyinvarsetnorm...
            ([baseline_ProbeStructure.PMIntensities(first_probe_ind:end),  ...
            ProbeStructure.PMIntensities(first_probe_ind:end)], 'baseline', 1);
    else
        NormData_probes = ProbeStructure.PMIntensities(first_probe_ind:end);
    end
    NormData = mean( (reshape(NormData_probes, num_pm_probes, length(NormData_probes) / num_pm_probes)))';
    eval(['save ' working_path sample_names{i} '_' chip '_norm.mat NormData allele_ratio_vec;']);

    % get summation of data of diploid samples
    if(i==diploid_ind(1))
        MeanNormData = NormData; 
    else
        for j = 1:length(diploid_ind)
            if(i==diploid_ind(j))
                MeanNormData = MeanNormData + NormData; 
            end
        end
    end
end

MeanNormData = MeanNormData ./ num_samples; % calc mean on diploid samples

for i=1:num_samples   % 3nd loop: transfer to application input format
    i

    eval(['load(''' working_path sample_names{i} '_' chip '_norm.mat'');']);
    copy_num_vec = (2 .* NormData ./ MeanNormData)'; % Just copy to the correct name and transpose!

    if(size(copy_num_vec,1) == 1)
        copy_num_vec = copy_num_vec';
    end
    eval(['allele_ratio_vec_' chip '=allele_ratio_vec;']);
    eval(['copy_num_vec_' chip '=copy_num_vec;']);
    eval(['snp_id_' chip '=snp_id;']);

    eval(['save ' fullfile(working_path, sample_names{i}) '_' chip '.mat ' ...
        'snp_id_' chip ' copy_num_vec_' chip ' allele_ratio_vec_' chip ';']);
end

