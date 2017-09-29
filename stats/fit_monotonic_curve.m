% Internal function - fit monotonic curve to data (interface to slm functions)
% Input:
% x - x-values
% y - y-values. Muse be monotonically decreasing?
% params - structure with fitting parameters
%
% Output:
% x_fit - values at which to calculate fitted y
% y_fit - fitted y values 
%
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

if(~isfield(params, 'fit_log'))
    params.fit_log = [0, 0]; % 0; default: no log on x and y 
end
prescription = slmset('plot','off', 'knots', params.knots);  % set parameters

if(isfield(params, 'direction'))
    switch params.direction
        case {'up', 'increasing'}
            prescription.Increasing = 'on';
        case {'down', 'decreasing'}
            prescription.Decreasing = 'on';
    end
end

if(isfield(params, 'min'))
    prescription.MinValue = params.min;
end
if(isfield(params, 'max'))
    prescription.MaxValue = params.max;
end
if(params.cum)
    y_vec = cumsum(y);
else
    y_vec = y;
end
if(params.fit_log(1)) % log-transform on x
    x_vec = double(log(max(0.5, x))); x_fit = double(log(max(0.5, x_fit)));
else
    x_vec = double(x);
end
if(params.fit_log(2)) % log-transform on y
    y_vec = log(y_vec); % double(log(max(0.5, x))); x_fit = double(log(max(0.5, x_fit)));
else
%    y_vec = double(x);
end


copy_fields = intersect(fields(prescription), fields(params)); % add all additional fields 
for i=1:length(copy_fields)
    eval(['prescription.' copy_fields{i} ' = params.' copy_fields{i} ';']);
end

slm = slmengine(x_vec, y_vec, prescription); %'plot','off','knots', params.knots, 'increasing','on', 'minvalue', 0); % 'concavedown', 'on',  % 'leftslope',0,'rightslope',0);
y_fit = slmeval(double(x_fit), slm);
if(params.cum)
    y_fit = [y_fit(1) diff(y_fit)]; % go back to density
end
if(params.fit_log(2))
    y_fit = exp(y_fit);
end

%else
%    slm = slmengine(double((x)), y, 'plot','off','knots', params.knots, 'decreasing','on', 'minvalue', 0); % 'leftslope',0,'rightslope',0);
%    y_fit = slmeval(double((x_fit)), slm);
%end
