function [OutputFilesNames NormalMeanIntensities r_mat] = ...
    PreLearnMoGChromFromData(user_dir, chip_type)

AssignAllGlobalConstants;
max_outlier_frac = 0.05;
% chip_type = lower(SNPChipAnnotStruct.chip); % 'hind'; % should be input ?

% load the matrix withh all the data for all the samples
load(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'array_outlier_vec', 'NormalMeanIntensities', ...
    'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB', 'AllData_M');
% check if the data was already moved
moved_flag = exist('NormalMeanIntensities')

if(~moved_flag)
    % First calculate the mean intensities of each sample
    [NormalMeanIntensities SampleNames mix_gauss_params] = GetNormalMeanIntensities(user_dir, chip_type, NormalizedSNPsCopyMatA, NormalizedSNPsCopyMatB);

    NormalMultIntensities = 2 ./ NormalMeanIntensities; smooth_window = 100;  % width of window for smoothing

    num_samples = length(NormalMeanIntensities);
    % Move the data to equate the mean to 2
    smoothed_data = zeros(size(NormalizedSNPsCopyMatA), 'single');
    for cur_sample = 1:num_samples
        doing_sample = cur_sample

        NormalizedSNPsCopyMatA(:,cur_sample) = NormalizedSNPsCopyMatA(:,cur_sample) .* NormalMultIntensities(cur_sample);
        NormalizedSNPsCopyMatB(:,cur_sample) = NormalizedSNPsCopyMatB(:,cur_sample) .* NormalMultIntensities(cur_sample);

        smoothed_data(:,cur_sample) = single(smooth(NormalizedSNPsCopyMatA(:,cur_sample) + NormalizedSNPsCopyMatB(:,cur_sample), smooth_window));

    end

    save(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'NormalizedSNPsCopyMatA', 'NormalizedSNPsCopyMatB', ...
        'NormalMeanIntensities','-append');
    % Learn a new MoG model for the data - only on samples with low %
    % outliers
    good_ind = find(array_outlier_vec<max_outlier_frac);
    smoothed_data = smoothed_data(:,good_ind);
    num_of_iterations =25; D=4; INIT_P = [0.05 0.85 0.05 0.05]; INIT_M = [1.35 2 2.5 3]; INIT_S = [0.1 0.1 0.1 0.1];
    [AllData_P, AllData_M, AllData_S, AllData_LogLike]=MixtureOfGaussiansGivenInitMex(smoothed_data(:),D,num_of_iterations, ...
        INIT_P, INIT_M, INIT_S); % make sure all are the same type

%    LL = MixtureOfGaussiansGetLikelihood(smoothed_data(:), AllData_P, AllData_M, AllData_S)

    % Plot data and learned gaussians
    labels_vec = {'x', 'Prob. Density'}; legends_vec = {'data', 'MoG fit'}; color_vec = {}; axes_vec = {};
    num_bins = 500;
    fig_handle = MixtureOfGaussiansDraw1dGaussians(smoothed_data(:), AllData_P, AllData_M, AllData_S, labels_vec, legends_vec, color_vec, axes_vec, num_bins);
    saveas(fig_handle, fullfile(user_dir, ['Raw_data_moved_histogram_' chip_type '.fig']));
    % save the SHIFTED matrix withh all the data for all the samples
    save(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']), 'AllData_M','AllData_S','AllData_P', 'mix_gauss_params','-append');

end
% What's next? perform a 2-d MoG EM ?
r_mat = zeros(4);
r_mat(2,:) = AllData_M ./ AllData_M(2);

OutputFilesNames = ''; % should write the file name in here
