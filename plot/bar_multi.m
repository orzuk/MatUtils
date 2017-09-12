% Plot bars one on top of each other
% 
% Input: 
% y - a matrix of values to plot 
% new_fig_flag - if to create new figure
% legend_vec - contains legends
% 
function bar_multi(x, y, new_fig_flag, legend_vec)

AssignGeneralConstants();
if(exist('new_fig_flag', 'var') && new_fig_flag)
    figure; hold on;
end
n = size(y,2); m = size(y,1);
if(isempty(x))
    x = 1:m;
end

y=cumsum(y,2);
bar(x, y(:,n), color_vec(1));
for i=n-1:-1:1
    bar(x, y(:,i), color_vec(n+1-i));
end

if(exist('legend_vec', 'var'))
    legend(legend_vec);
end
