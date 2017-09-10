% Draw a histogram, such that the total area is one.
% The function also enables to draw many overlapping histograms together
% It also gets rid of nans in the data 
%
% Input:
% x_vec - vector or matrix of values to plot (can also be cell-array)
% bins (optional) - number of bins or vector of bins
% color_vec (optional) - which color to plot each histogram
% plot_flag (optional) - whether to plot (1, default), bar (2) or not to plot at all (0)
% scaling (optional) - enable the histogram sum to be different from one
% do_smooth (optional) - perform gaussian smoothing to make density look nicer
% style_vec (optional) - plot style (e.g. dashed) 
%
% Output: 
% h - frequency in each bin
% bin_locs - bins centers 
% 
function [h, bins_loc] = ...
    hist_density(x_vec, bins, color_vec, plot_flag, scaling, do_smooth, style_vec, varargin)

if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 1;
end

if(~iscell(x_vec))
    x_vec = vec2column(x_vec);
    m = size(x_vec, 2);
    x = cell(m,1);
    for i=1:m
        x{i} = x_vec(:,i);
    end
    x_vec = x;
end
if(~exist('color_vec', 'var') || isempty(color_vec))
    color_vec = 'brgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmc';
end
if(~exist('do_smooth', 'var'))
    do_smooth = 0; 
end
m = length(x_vec); 
if(~exist('style_vec', 'var') || isempty(style_vec))
    style_vec = '-';
end
if(plot_flag)
    hold on;
end
if(length(bins) == 1)
    bins = repmat(bins, m, 1);
end
h = cell(1,m); bins_loc = cell(1,m); % get a cell array of all distributions
% if(exist('bins', 'var'))
%     if(~iscell(bins))
%         bins = {bins};
%     end
% end
for i=1:m
    x_vec{i} = x_vec{i}(~isnan(x_vec{i})); % remove nans
    x_vec{i} = x_vec{i}(abs(x_vec{i}) < Inf); % remove infs
    if(exist('bins', 'var'))
        if(is_row(bins))  % set of bins
            [h{i}, bins_loc{i}] = hist(x_vec{i}, bins);
        else % number of bins
            [h{i}, bins_loc{i}] = hist(x_vec{i}, bins(i));
        end
    else
        [h{i}, bins_loc{i}] = hist(x_vec{i});
    end
    % bin_size = bins_loc{i}(2)-bins_loc{i}(1); % size of bins
    if(do_smooth) % smooth histograms
        h{i} =  gauss_smooth(h{i}, do_smooth, 1); % assume bins_loc are equally spaced !!!
    end
    h{i} = normalize_hist(bins_loc{i}, h{i}); % normalize to sum to one
%    h{i} = h{i} ./ (sum(h{i}) * bin_size); % normalize to sum to one
    if(exist('scaling', 'var') && (~isempty(scaling)) )
        h{i} = scaling .* h{i};
    end

    switch plot_flag
        case 1 % just plot
            plot(bins_loc{i}, h{i}, [color_vec(i) style_vec], 'LineWidth', 2); % draw a thick line
        case 2 % bars
            bar(bins_loc{i}, h{i}, color_vec(i)); % how to set the color
        otherwise
            % do nothing ..
    end
    bins_loc{i} = vec2row(bins_loc{i}); % make bins vectors always row vectors 
    
end
if(m == 1) % for one hist do not return a cell array
    h = h{1};
    bins_loc = bins_loc{1};
end


    