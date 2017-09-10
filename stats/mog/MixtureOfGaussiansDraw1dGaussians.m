% Plot histogram and graph for 1d MoG model. 
% This is only for one-dimensional Gaussians
%
% The input: 
% data - vector of data generated
% PVecs - prior probs.
% MuVecs - vectors of means.
% SigmaMats - standard deviation matrices.
% labels_vec - labels for x and y axis
% legends_vec - legends of gaussians
% color_vec - colors of gaussians (optional)
% axes_vec - plot axes (optional)
% num_bins - number of bins used in Gaussian histogram
% new_fig_flag - if to start a new fig (default) or not
% plot_data_flag - if to plot data (default, then everything scaled to #points) or only MoG distributions
%
function fig_handle = MixtureOfGaussiansDraw1dGaussians(data, PVecs, MuVecs, SigmaMats, ...
    labels_vec, legends_vec, color_vec, axes_vec, num_bins, new_fig_flag, plot_data_flag, varargin)

if(~exist('new_fig_flag', 'var'))
    new_fig_flag=1;
end
if(isempty(new_fig_flag))
    new_fig_flag=1;
end
if(~exist('plot_data_flag', 'var'))
    plot_data_flag=1;
end
if(new_fig_flag)
    fig_handle = figure; 
else
    fig_handle = gcf;
end
hold on;
    
n_points = length(data);
[height,bin_loc]=hist(data,num_bins);%bin_loc is Fisher_Zs_possible_vec
if(plot_data_flag)
    hist(data, num_bins);
    
end

xlabel('x'); ylabel('freq');

% Plot the resulting mixture model
num_gaussians = length(SigmaMats); y=zeros(num_gaussians, length(bin_loc)); 
for i=1:num_gaussians
    y(i,:)=PVecs(i)*1/(sqrt(2*pi)*SigmaMats(i))*exp(-(bin_loc-MuVecs(i)).^2/(2*SigmaMats(i)^2));
end
g_min = min(data);  g_max = max(data); g_gap = g_max-g_min;
x_vec = [g_min :g_gap*(1/num_bins):g_max-g_gap*(1/num_bins)];
MOG_fit_vec = sum(y,1);

%h{i} = h{i} ./ (sum(h{i}) * bin_size); % normalize to some to one

if(isempty(color_vec))
    color_vec = 'gr';
end
if(length(color_vec) < num_gaussians+1)
    color_vec = [repmat(color_vec(1), 1, num_gaussians) color_vec(2)];
end
bin_size = bin_loc(2)-bin_loc(1);
for i=1:num_gaussians
    if(~plot_data_flag)
        n_points = 1/bin_size;%      PVecs(i)/sum(y(i,:));
    end
    plot(x_vec, y(i,:) .* n_points .* bin_size, color_vec(i), 'LineWidth', 2);
end
% if(~plot_data_flag)
%     n_points = 1 / sum(MOG_fit_vec); 
% end
plot(x_vec, MOG_fit_vec .* n_points .* bin_size, color_vec(num_gaussians+1), 'LineWidth', 2);

legend(legends_vec); xlabel(labels_vec{1}); ylabel(labels_vec{2});
