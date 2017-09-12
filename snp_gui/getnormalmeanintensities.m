% This function goes over all the samples and find for each sample what is
% the mean intensity when the copy number is 'normal' - 2.
function [NormalMeanIntensities SampleNames mix_gauss_params] = GetNormalMeanIntensities(user_dir, chip_type)

AssignAllGlobalConstants;

% load the matrix withh all the data for all the samples
ALL=load(fullfile(user_dir, ['AllSamplesMat_' chip_type '.mat']));


num_samples = length(ALL.SampleNames); SampleNames = ALL.SampleNames;
smooth_window = 100;  % width of window for smoothing
num_of_iterations =25;
num_starting_points = 10;
max_num_of_Gaussians = 4;
min_normal_gauss = 1.65; % distance of normal from 2
max_normal_gauss = 2.15;
min_normal_p = 0.3;
max_pair_normal_distance = 0.1; % 0.05
mix_gauss_params = [];
NormalMeanIntensities = zeros(1, num_samples);
labels_vec = {'x', 'Prob. Density'}; legends_vec = {'data', 'MoG fit'}; color_vec = {}; axes_vec = {};
for cur_sample = 1:num_samples
    doing_sample = cur_sample
    smoothed_data = smooth(ALL.NormalizedSNPsCopyMatA(:,cur_sample) + ALL.NormalizedSNPsCopyMatB(:,cur_sample), smooth_window);

    [P,M,S, Dim, LogLike] = MixtureOfGaussiansFindModelDimension(smoothed_data, max_num_of_Gaussians, ...
        num_of_iterations, num_starting_points);
    mix_gauss_params.S{cur_sample} = S; mix_gauss_params.M{cur_sample} = M;
    mix_gauss_params.P{cur_sample} = P; mix_gauss_params.Dim{cur_sample} = Dim;
    mix_gauss_params.LogLike{cur_sample} = LogLike;


    % Look for the maximum
    if(length(M) == 1)
        NormalMeanIntensities(cur_sample) = M;
    else
        num_x_vals = 1000000; res = (max(M) - min(M)) * 1.1 / num_x_vals;

        x_vec = [min(M)- (max(M) - min(M))*0.05:res:max(M) + (max(M) - min(M))*0.05];
        y_vec = zeros(1,length(x_vec));
        for i=1:length(M)
            y_vec = y_vec + (P(i) ./ (sqrt(2*pi)*S(i))) .* exp(-(x_vec - M(i)).^2 ./ (2.*S(i)^2));
        end
        y_diff = diff(y_vec);
        max_inds_vec = find( (y_diff(1:end-1) > 0) & (y_diff(2:end) < 0) );
        if(~isempty(max_inds_vec))
            max_x_vec = x_vec(max_inds_vec);
            max_y_vec = y_vec(max_inds_vec);
            max_close_inds = find((max_x_vec >min_normal_gauss) & (max_x_vec < max_normal_gauss) );
            if(~isempty(max_close_inds))
                max_inds_vec = max_inds_vec(max_close_inds);
                [max_val max_ind] = max( max_y_vec(max_close_inds) );
                NormalMeanIntensities(cur_sample) = x_vec(max_inds_vec(max_ind));
            else
                [tmp_min_val tmp_min_ind] = min( min( abs(max_x_vec - min_normal_gauss), abs(max_x_vec - max_normal_gauss) ) );
                NormalMeanIntensities(cur_sample) = x_vec(max_inds_vec(tmp_min_ind(1)));
                TTT = 13123123;
            end
        else
            TTTT = 123123;
        end
    end

    % Check which Gaussians are close to each other
    % % %     diff_inds = find(diff(M) < max_pair_normal_distance);
    % % %     if(~isempty(diff_inds))
    % % %         for i=length(diff_inds)
    % % %             M(diff_inds(i)) = ( P(diff_inds(i)) * M(diff_inds(i)) + P(diff_inds(i)+1) * M(diff_inds(i)+1) ) / ( P(diff_inds(i)) + P(diff_inds(i)+1) );
    % % %             S(diff_inds(i)) = ( P(diff_inds(i)) * S(diff_inds(i)) + P(diff_inds(i)+1) * S(diff_inds(i)+1) ) / ( P(diff_inds(i)) + P(diff_inds(i)+1) );
    % % %             P(diff_inds(i)) = P(diff_inds(i)) + P(diff_inds(i)+1);
    % % %             M = M(setdiff(1:length(M), diff_inds(i)+1));
    % % %             S = S(setdiff(1:length(S), diff_inds(i)+1));
    % % %             P = P(setdiff(1:length(P), diff_inds(i)+1));
    % % %             if(i < length(diff_inds))
    % % %                 diff_inds(i+1:end) = diff_inds(i+1:end)-1;
    % % %             end
    % % %         end
    % % %     end
    % % % %%    MixtureOfGaussiansDraw1dGaussians(smoothed_data, P, M, S, labels_vec, legends_vec, color_vec, axes_vec);
    % % %     big_p_ind = find(P>min_normal_p & M > min_normal_gauss);
    % % %     if(~isempty(big_p_ind))
    % % %         [P_max P_ind] = min(M(big_p_ind)); P_ind = big_p_ind(P_ind);
    % % %     else
    % % %         %close_inds = find(abs(M-2) < max_normal_distance);
    % % %         close_inds = find(M > min_normal_gauss & M < max_normal_gauss);
    % % %         [P_max P_ind] = max(P(close_inds)); P_ind = close_inds(P_ind);
    % % %     end
    % % % %%    max_inds = find( abs(M-M(P_ind)) < max_pair_normal_distance);
    % % %     max_inds = P_ind;
    % % %
    % % %     NormalMeanIntensities(cur_sample) = sum(M(max_inds) .* P(max_inds)) ./ sum(P(max_inds));

end