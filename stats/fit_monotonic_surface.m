% Internal function - fit monotonic surface to data (slm or gridfit interface - iterative)
% Constraints (here are very specific - need to generalize them later):
% z(x, y) DECREASING with x for any fixed y
% \int_{t=0}^x z(t, y)dt INCREASING with y for any fixed x
%
% Input:
% x - vector of x coordinates
% y - vector of y coordinates
% z - vector of function of x and y: z(x,y)
%
function [x_fit, y_fit, z_fit] = fit_monotonic_surface(x, y, z, params) % constraints

num_x = length(x);
num_y = length(y);

if(~isfield(params, 'plot')) 
    params.plot = 0; 
end

if(~isfield(params, 'x_fit'))
    x_fit = x;
else
    if(isscalar(params.x_fit)) % number of points - set logarithmic grid
        x_fit = logspace(log10(min(x)), log10(max(x)), params.x_fit);
    else % points given
        x_fit = params.x_fit;
    end
end
if(~isfield(params, 'y_fit'))
    y_fit = y;
else
    if(isscalar(params.y_fit)) % number of points - set logarithmic grid
        y_fit = logspace(log10(min(y)), log10(max(y)), params.y_fit);
    else % points given
        y_fit = params.y_fit;
    end
end

z_fit = zeros(length(y_fit), length(x_fit)); z_fit(2,1)=-1; % create z grid for fitted data
z_fit0 = zeros(size(z));
for i=1:num_y % First fit each y seperately monotonically
    [~, z_fit0(i,:)] = fit_monotonic_curve(x, z(i,:), params);
end
z_fit0 = normalize(z_fit0, 2); % set to sum to one
z_fit_cum0 = cumsum(z_fit0, 2); % here take cumsum of ROWs
%z_fit_cum = z_fit;
params.x_fit = y_fit; % switch roles
[x_mesh, y_mesh] = meshgrid(x_fit, y_fit);
z_fit_cum = interp2(double(x), y, z_fit_cum0, double(x_mesh), y_mesh);
%for i=1:num_x % First fit each y seperately monotonically
%    if(mod(i, 100)==0)
%        run_i = i
%    end
%    [~, z_fit_cum(:,i)] = fit_monotonic_curve(y, -z_fit_cum0(:,i), params);
%end
%z_fit_cum=-z_fit_cum;
params.x_fit = x_fit; % get back to x
z_fit = max(0, [z_fit_cum(:,1) diff(z_fit_cum, [], 2)]); % update surface
z_fit = normalize(z_fit, 2); % set to sum to one

%issortedtol(x, tol, direction)
[num_y, num_x] = size(z_fit); % update size
tol=10e-6;
while(~(issortedtol(z_fit, tol, 'decreasing') && issortedtol(z_fit_cum', tol)))  % add tolerance
    params.x_fit = x_fit; % switch roles
    for i=1:num_y % First fit each y seperately monotonically
        [~, z_fit(i,:)] = fit_monotonic_curve(x_fit, z_fit(i,:), params);
    end
    z_fit_cum = cumsum(z_fit, 2); % get cumulative of ROWs
    z_fit_cum  = cummax(z_fit_cum, 1);
    %     params.x_fit = y_fit; % switch roles
    %     for i=1:num_x % next fit for each x seperately monotonically
    %         if(mod(i, 100)==0)
    %             run_i = i
    %         end
    %         [~, z_fit_cum(:,i)] = fit_monotonic_curve(y_fit, -z_fit_cum(:,i), params);
    %     end
    %     z_fit_cum=-z_fit_cum;
    z_fit = max(0, [z_fit_cum(:,1) diff(z_fit_cum, [], 2)]); % update surface
    z_fit = normalize(z_fit, 2); % set to sum to one
end % while not sorted
params.x_fit = x_fit; % get back to x


% Next: draw surface:
% Plot data
if(params.plot)
    [x_g, y_g] = meshgrid(x, y);
    figure; plot3(x_g(:), y_g(:), z(:), '+');
    %hold on;
    %plot3(x, y, z, '.'); % plot original points
    
    % Plot surface
    figure;
    surf(log(1+x_fit),y_fit,z_fit_cum);
    shading interp
    colormap(jet(256))
    camlight right
    lighting phong
    xlabel('f'); ylabel('s'); zlabel('\Psi_f');
end % if plot_flag
% for i=1:num_y
%     if(~issorted(-z_fit(i,:)))
%         rrr = i
%     end
% end
% for i=1:num_x
%     if(~issorted(-z_fit_cum(:,i)))
%         ccc = i
%     end
% end






