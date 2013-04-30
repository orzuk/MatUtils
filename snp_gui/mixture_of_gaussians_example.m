function mixture_of_gaussians_example()

addpath(fullfile('..','src'));
% load('E:\Libi\tools\SNP_tool\data\Leukemia_all_hg18\mixture_gaussians_example_xba.mat', ...
%     'smoothed_data', 'max_num_of_Gaussians', 'num_of_iterations', 'num_starting_points', 'sample_name');

load('mixture_gaussians_example_xba.mat', ...
    'smoothed_data', 'max_num_of_Gaussians', 'num_of_iterations', 'num_starting_points', 'sample_name');

figure; hist(smoothed_data, 500); title(sample_name);
[P,M,S, Dim, LogLike] = MixtureOfGaussiansFindModelDimension(smoothed_data, max_num_of_Gaussians, ...
	num_of_iterations, num_starting_points);
MixtureOfGaussiansDraw1dGaussians(smoothed_data, P, M, S, ...
    {'1', '2', '3'}, {'1', '2', '3'}, 'rgb', [], 500);

M % doesn't make sense
t=1;
