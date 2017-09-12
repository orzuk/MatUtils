% Internal function - fit monotonic curve to data (slm interface)
% Input: 
% x - x-values
% y - y-values. Muse be monotonically decreasing? 
function [x_fit, y_fit] = fit_monotonic_curve(x, y, params)
% use slm tool - fit a decreasing function
if(~isfield(params, 'x_fit'))
    x_fit = x;
else
    if(isscalar(params.x_fit)) % number of points - set logarithmic grid
        x_fit = logspace(log10(min(x)), log10(max(x)), params.x_fit);
    else % points given
        x_fit = params.x_fit;
    end
end
slm = slmengine(double(log(x)), y,'plot','off','knots', params.knots, 'decreasing','on'); % 'leftslope',0,'rightslope',0);
y_fit = slmeval(double(log(x_fit)), slm);

