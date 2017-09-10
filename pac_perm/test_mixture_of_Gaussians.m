% Test the mixture_of_Gaussians function

% Generate data
num_g = 2;
miu = [0, 1];
sigma = [1, 0.2];
prior = [0.5 0.5];
n_points = 2000;

g_data = randn(1,n_points);
ind_pick = rand(1,n_points) > prior(1);
g_data(find(ind_pick)) = g_data(find(ind_pick)) * sigma(2);
g_data(find(ind_pick)) = g_data(find(ind_pick)) + miu(2);

iters = 100;


% Plot the data
[height,bin_loc]=hist(g_data,100);%bin_loc is Fisher_Zs_possible_vec
figure; hold on;  hist(g_data, 100); xlabel('x'); ylabel('freq');

num_g = 2; % Give more gaussians

[S,M,P, LogLike]=mixture_of_Gaussians(g_data',num_g,iters,111) % last argument is dummy

% Plot the resulting mixture
clear y;
for i=1:num_g
    y(i,:)=P(i)*1/(sqrt(2*pi)*S(i))*exp(-(bin_loc-M(i)).^2/(2*S(i)^2));
end

g_min = min(g_data);  g_max = max(g_data); g_gap = g_max-g_min;


x_vec = [g_min :g_gap*0.01:g_max-g_gap*0.01];

MOG_fit_vec = sum(y,1);
plot(x_vec, MOG_fit_vec .* n_points .* (bin_loc(2)-bin_loc(1)), 'r');




