% Internal function - fit monotonic surface to data (slm or gridfit interface - iterative)
% Constraints (here are very specific - need to generalize them later):
% z(x, y) DECREASING with x for any fixed y
% \int_{t=0}^x z(t, y)dt INCREASING with y for any fixed x
%
% Input:
% x - vector of x coordinates
% y - vector of y coordinates
% z - vector of function of x and y: z(x,y)
% params - structure with smoothing parameters, indluding x_fit, y_fit
%
% Output:
% x_fit - vector of fitted x coordinates
% y_fit - vector of fitted y coordinates
% z_fit - matrix of function of x and y: z(x,y)
%
function [x_fit, y_fit, z_fit] = fit_monotonic_surface(x, y, z, params) % constraints

num_x = length(x); num_y = length(y); num_z = length(z);


if(~isfield(params, 'plot'))
    params.plot = 0;
end

if(~isfield(params, 'x_fit'))
    x_fit = x; params.x_fit = x_fit;
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

if((num_z == num_x) && min(size(z))==1) % here we're given vectors of x,y,z of the same size
    if(isscalar(y))
        y = repmat(y, size(x));
    end
    x_unique = unique(x); num_x = length(x_unique);
    y_unique = unique(y); num_y = length(y_unique);
    z_fit0 = zeros(num_y, num_x);
    params.x_fit = x_unique; x_fit = x_unique;
    one_vec = 1;
else
    z_fit0 = zeros(size(z));
    one_vec = 0;
end


%z_fit = zeros(length(y_fit), length(x_fit)); % z_fit(2,1)=-1; % create z grid for fitted data
for i=1:num_y % First fit each y seperately monotonically
    params.cum = 0; params.fit_log = [1 0]; params.min = realmin; params.RightMinValue = realmin;  % should we give up direction constraint? params.direction = 'down'; % params.direction = 'up';
    if(one_vec) % here we're given vectors of x,y,z of the same size
        I = find(y == y_unique(i)); % get indices
        fit_x_vec = x(I); fit_z_vec = z(I) .* double(max(1,x(I))); fit_z_vec_unweighted = z(I);
    else
        fit_x_vec = x; fit_z_vec = z(i,:) .* double(max(1,x)); fit_z_vec_unweighted = z(i,:);
    end
    %       [~, z_fit0(i,:)] = fit_monotonic_curve(x(I), z(I), params);
    [~, z_fit0(i,:)] = fit_monotonic_curve(fit_x_vec, fit_z_vec, params);
    z_fit0(i,:) = max(realmin,  z_fit0(i,:) ./ double(max(1,params.x_fit)));      % normalize
    max_ind = find(fit_z_vec_unweighted>0, 1, 'last'); max_val = fit_x_vec(max_ind); % x(I(max_ind)); % find last value
    fit_again = 1;
    if(fit_again)
        max_ind = min(find(params.x_fit > max_val, 1), size(z_fit0, 2)-1);
        max_ind2 = find(diff(z_fit0(i,2:end))>0, 1);
        if(~isempty(max_ind2))
            if(isempty(max_ind))
                max_ind = max_ind2; 
            else
                max_ind = min(max_ind, max_ind2);
            end
        end        
        if(~isempty(max_ind))
            params.cum = 0; params.direction = 'down'; params.min = 0; params.fit_log = [1 1]; params.RightMinValue = 0; x_fit = params.x_fit; params.x_fit = params.x_fit(max_ind:end);
            %            [~, z_fit0(i,max_ind:end)] = fit_monotonic_curve(x_fit(2:end), z_fit0(i,2:end), params);  % fit again       % force z to be monotonic
            deg = 5;
            ppp = polyfit(log(x_fit(2:max_ind)), log(z_fit0(i,2:max_ind)), deg); % fit quadratic
%             figure;
%             loglog(fit_x_vec ./ max(x), (fit_z_vec_unweighted)); hold on;
%             loglog(x_fit ./ max(x_fit), z_fit0(i,:), 'r', 'linewidth', 2);

            z_fit0(i,max_ind:end) = 0; zz1 = zeros(size(x_fit)); 
            for jj=0:deg
                z_fit0(i,max_ind:end) = z_fit0(i,max_ind:end) + ppp(deg+1-jj) .* log(params.x_fit).^jj;
                zz1 = zz1 + ppp(deg+1-jj) .* log(x_fit).^jj;
            end
            z_fit0(i,max_ind:end) = exp(z_fit0(i,max_ind:end)); zz1 = exp(zz1); 
%                exp(ppp(3) + ppp(2) .* log(params.x_fit) + ppp(1) .* log(params.x_fit).^2 );
            
            % new fit 
%             gaussEqn = 'a*exp(b*x)+d*x^c' % +f*x^g'
%             p_fit = fit(log(x_fit(2:max_ind))', log(z_fit0(i,2:max_ind))', gaussEqn); 
%             z_fit00 = exp(feval(p_fit, log(params.x_fit))); 
%             loglog(params.x_fit ./ max(params.x_fit), z_fit00, 'r--');  hold on; 
             params.x_fit = x_fit; % return back
% %              loglog(params.x_fit ./ max(params.x_fit), zz1, 'g--', 'linewidth', 2);
% %              title(['s = ' num2str(y_unique(i))]); legend('data', 'fit0', 'fit1'); 
% %              loglog(fit_x_vec(max_ind) ./ max(x), max(fit_z_vec_unweighted(max_ind), min(zz1(zz1>0))), '*k'); 
% %              loglog(params.x_fit ./ max(params.x_fit), z_fit0(i,:), 'm-.', 'linewidth', 2);
        end
    end % fit again
    if(params.plot)
        figure; semilogx(fit_x_vec ./ max(x), cumsum(fit_z_vec));
        hold on; semilogx(params.x_fit ./ max(params.x_fit), cumsum(z_fit0(i,:) .* double(max(1,params.x_fit))), 'r--');
        legend({['original, s=' num2str(y_unique(i))], ['fitted, s=' num2str(y_unique(i))]}, 'location', 'southeast');
    end
    % % %     else
    % % %         %    [~, z_fit0(i,:)] = fit_monotonic_curve(x, z(i,:), params);
    % % %         [~, z_fit0(i,:)] = fit_monotonic_curve(x, z(i,:) .* double(max(1,x)), params);
    % % %         z_fit0(i,:) = z_fit0(i,:) ./ double(max(1,x));
    % % %         if(params.plot)
    % % %             figure; semilogx(x, cumsum(z(i,:) .* double(max(1,x))));
    % % %             hold on; semilogx(x, cumsum(z_fit0(i,:) .* double(max(1,x))), 'r--');
    % % %             legend(['original, s=' num2str(y(i))], ['fitted, s=' num2str(y(i))]);
    % % %         end
    % % %     end
end % loop on y
z_fit0 = normalize(z_fit0, 2); % set to sum to one
z_fit_cum0 = cumsum(z_fit0, 2); % here take cumsum of ROWs
%z_fit_cum = z_fit;
params.x_fit = y_fit; % switch roles
[x_mesh, y_mesh] = meshgrid(x_fit, y_fit);
if(one_vec)
    if(length(y_unique)>1)
        z_fit_cum = interp2(double(x_unique), y_unique, z_fit_cum0, double(x_mesh), y_mesh);
        z_fit = interp2(double(x_unique), y_unique, z_fit0, double(x_mesh), y_mesh);
    else
        z_fit_cum = z_fit_cum0; z_fit = z_fit0; 
    end
else
    z_fit_cum = interp2(double(x), y, z_fit_cum0, double(x_mesh), y_mesh);
    z_fit = interp2(double(x), y, z_fit0, double(x_mesh), y_mesh);
end
params.x_fit = x_fit; % get back to x
%z_fit = max(0, [z_fit_cum(:,1) diff(z_fit_cum, [], 2)]); % update surface
z_fit = normalize(z_fit, 2); % set to sum to one
z_fit = max(z_fit, realmin); 
% Add monotonicity in s:
% z_fit_cum = cumsum(z_fit, 2); % get cumulative of ROWs
% z_fit_cum  = cummax(z_fit_cum, 1); % enforce monotonicity with s
% z_fit = max(0, [z_fit_cum(:,1) diff(z_fit_cum, [], 2)]); % update surface
% z_fit = normalize(z_fit, 2); % set to sum to one

return;  % up to here simple fitting - no row-column joint information. Fitting looks good (?)

[num_y, num_x] = size(z_fit); % update size
tol=10e-6; ctr=1;
while(~(issortedtol(z_fit, tol, 'decreasing') && issortedtol(z_fit_cum', tol)))  % add tolerance
    params.x_fit = x_fit; % switch roles
    for i=1:num_y % First fit each y seperately monotonically
        [~, z_fit(i,:)] = fit_monotonic_curve(x_fit, z_fit(i,:) .* double(max(1,x_fit)), params); z_fit(i,:) = z_fit(i,:) ./ double(max(1,x_fit)); % NEW! force f*psi(f) to be monotonic
        %        [~, z_fit(i,:)] = fit_monotonic_curve(x_fit, z_fit(i,:), params);
    end
    z_fit = normalize(max(0, z_fit), 2); % set to sum to one
    
    %   return;
    z_fit_cum = cumsum(z_fit, 2); % get cumulative of ROWs
    z_fit_cum  = cummax(z_fit_cum, 1); % enforce monotonicity with s
    
    z_fit = max(0, [z_fit_cum(:,1) diff(z_fit_cum, [], 2)]); % update surface
    z_fit = normalize(z_fit, 2); % set to sum to one
    
    %   return % finish before iterative procedure !!
    run_iter = ctr
    ctr=ctr+1;
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






