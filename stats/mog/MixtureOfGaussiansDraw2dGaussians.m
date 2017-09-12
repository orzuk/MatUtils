% Plot elipses representing 2d MOG components. 
% This is only for two-dimensional Gaussians
%
% The input:
% MuVecs - vector of means
% SigmaMats - st.d. matrices (cell arrays)
% labels_vec - a vector of labels
% legends_vec - a vector of legends
% color_vec - gaussians colors
% centers_flag - plot also gaussian centers (optional. Default is 0)
%
function MixtureOfGaussiansDraw2dGaussians(MuVecs, SigmaMats, ...
    labels_vec, legends_vec, color_vec, axes_vec, centers_flag, varargin)

num_gaussians = length(SigmaMats);
if(~exist('color_vec', 'var'))
    color_vec = 'bgrkmc';
end
if(isempty(color_vec))
    color_vec = 'bgrkmc';
end

if(~exist('centers_flag', 'var'))
    centers_flag = 0;
end
if(isempty(centers_flag))
    centers_flag = 0;
end


%figure;
hold on;

theta = [0:0.05:2*pi]'; x = [cos(theta) sin(theta)]';
for m=1:num_gaussians
    cur_color = color_vec(mod(m-1, 6)+1);
    A = sqrtm(SigmaMats{m});
    y = A * x + repmat(MuVecs(m,:),  length(x), 1)';
    if(centers_flag) % plot also centers
        y = [y MuVecs(m,:)'];
    end
    plot(y(1,:), y(2,:), [cur_color '.']);
end

% if(centers_flag) % plot also centers
%     for m=1:num_gaussians
%         cur_color = color_vec(mod(m-1, 6)+1);
%         plot(MuVecs(m,1), MuVecs(m,2), ['*' cur_color]);
%     end
% end

if(~isempty(labels_vec))
    xlabel(labels_vec{1}); ylabel(labels_vec{2});
end

if(~isempty(legends_vec))
    legend(legends_vec);
end
title('Gaussians Centers and 1 St.d.s');

if(exist('axes_vec', 'var'))
    if(~isempty(axes_vec))
        axis(axes_vec);
    end
end
