% Set axis to log scale.
% Matlab problem: loglog sometimes forget it's on log scale and goes back
% to standard scale
%
% Input:
% axes_str - string represnting axis (x,y)
% log_base - basis of the log (default is e)
%
function my_loglog(axes_str, log_base)


if(~exist('log_base', 'var') || isempty(log_base)) % set default: natural logarithm
    log_base = exp(1); 
end
new_ticks = 1; % how to display ticks on log scale 
for i=1:length(axes_str)
    switch axes_str(i)
        case {'X', 'x'}% change x axis ticks
            if(~new_ticks)
                %             xlabels = get(gca, 'XTickLabel');
                %             xlabels = num2str(exp(str2num(xlabels)), 3);
                %             set(gca, 'XTickLabel', xlabels);
                
            else            % New: make 'nice' numbers on log scale
                x_ticks = get(gca, 'XTick');
                num_ticks = length(x_ticks);
                tick_range = [get_first_digits(log_base.^(x_ticks(1)), 'down') ...
                    get_first_digits(log_base.^(x_ticks(end)), 'up')];
                tick_diff = (tick_range(2) - tick_range(1)) / (num_ticks-1);
                new_x_ticks = tick_range(1):tick_diff:tick_range(2);
                set(gca, 'XTick', log(new_x_ticks) ./ log(log_base));
                xlabels = num2str_cell(num2cell(new_x_ticks));
                set(gca, 'XTickLabel', xlabels);
            end
        case {'Y', 'y'} % change y axis ticks
            if(~new_ticks)
                ylabels = get(gca, 'YTickLabel');
                ylabels = num2str(exp(str2num(ylabels)), 3);
                set(gca, 'YTickLabel', ylabels);
            else % New: make 'nice' numbers on log scale
                y_ticks = get(gca, 'YTick');
                num_ticks = 11; % always use 10 different ticks length(y_ticks);
                tick_range = [get_first_digits(log_base.^(y_ticks(1)), 'down') ...
                    get_first_digits(log_base.^(y_ticks(end)), 'up')];
                tick_diff = (tick_range(2) - tick_range(1)) / (num_ticks-1);
                new_y_ticks = tick_range(1):tick_diff:tick_range(2);
                set(gca, 'YTick', log(new_y_ticks)./log(log_base));
                ylabels = num2str_cell(num2cell(new_y_ticks));
                set(gca, 'YTickLabel', ylabels);
            end
    end
end

