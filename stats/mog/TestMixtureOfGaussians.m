% A script for testing the MixtureOfGaussians package functions.
% First generate data from a MOG distribution, and then fit the sample
% distribution to a MOG model. Finally plot the data and the output model.
%
% Written by Liat Ein-Dor and Or Zuk 10/2006
%

% Generate data
num_g = 2; % # gaussians
mu = [0, 1]; % means
sigma = [1, 0.2]; % std.s 
prior = [0.5 0.5]; % probs. for each gaussian
n_points = 500;
g_data = randn(1,n_points);
ind_pick = rand(1,n_points) > prior(1);
g_data(find(ind_pick)) = g_data(find(ind_pick)) * sigma(2);
g_data(find(ind_pick)) = g_data(find(ind_pick)) + mu(2);

iters = 100; % Number of iterations for the EM algorithm
start_points =20; % Number of different random starting points for the EM algorithm
num_bins = 100; 

% Plot  data
[height,bin_loc]=hist(g_data,num_bins); 
figure; hold on;  hist(g_data, num_bins); xlabel('x'); ylabel('freq');

[P,M,S, LogLike] = MixtureOfGaussiansEM(g_data',num_g,iters,start_points) % % Learn model using EM. Last argument is dummy

% Plot the resulting mixture model
clear y; y = zeros(num_g, num_bins); 
for i=1:num_g
    y(i,:)=P(i)*1/(sqrt(2*pi)*S(i))*exp(-(bin_loc-M(i)).^2/(2*S(i)^2));
end
g_min = min(g_data);  g_max = max(g_data); g_gap = g_max-g_min;
x_vec = [g_min :g_gap*0.01:g_max-g_gap*0.01];
MOG_fit_vec = sum(y,1);
for i=1:num_g
   plot(x_vec, y(i,:) .* n_points .*  (bin_loc(2)-bin_loc(1)), 'g', 'LineWidth', 2);
end
plot(x_vec, MOG_fit_vec .* n_points .* (bin_loc(2)-bin_loc(1)), 'r', 'LineWidth', 2);
legend_vec = cell(num_g+2,1); legend_vec{1} = 'Data'; legend_vec{end} = 'Mog Fit';
for i=1:num_g
    legend_vec{i+1} = ['Gauss. ' num2str(i)]; 
end
legend(legend_vec); 


% Now try learning a multi-dimensional model 
dim = 2; num_g = 3;
mu_2d = [0.5 3; 3 0.5; 2 2];
sigma_2d{1} = [0.07 0; 0 1.7]; sigma_2d{2} = [ 1.7 0.017; 0.017 0.07]; sigma_2d{3} = [ 0.7 -0.01; -0.01 0.7]; prior = [0.35 0.35 0.3];

g_data = MixtureOfGaussiansMultiDimSimulateData(prior, mu_2d, sigma_2d, n_points); 
[S,M,P, LogLike]=MixtureOfGaussiansMultiDimGivenInit(g_data,num_g,iters, (1-prior) ./ sum( 1-prior), mu_2d+2*rand(num_g,dim),  sigma_2d);
C = MixtureOfGaussiansMultiDimClassify(g_data, P, M, S); 

% Plot the data and the  MOG models 
colorvec = 'bgrmcxok:.'; labels_vec = []; legends_vec = []; 
figure; hold on; 
for m=1:num_g
    plot(g_data(1,find(C == m)), g_data(2,find(C == m)), [colorvec(m), '.']);
end
MixtureOfGaussiansDraw2dGaussians(mu_2d, sigma_2d, labels_vec, legends_vec, 'kkk',[],1); % draw original model 
MixtureOfGaussiansDraw2dGaussians(M, S, labels_vec, legends_vec, 'ccc',[],1); % draw learned model
legend('Class 1', 'Class 2', 'Class 3', ...
    'True Model 1', 'True Model 2', 'True Model 3', ...
    'Learned Model 1', 'Learned Model 2', 'Learned Model 3');

