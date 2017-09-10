% Plot average selection coefficient as a function of effect size
%
% Input:
% selection_func_str - what kind of selection coefficient we use
% params_vec - parameter of dependence of selection by effect size
% fig_outfile - where to save output file
%
function plot_stabilizing_selection(selection_func_str, params_vec, fig_outfile)

AssignGeneralConstants;
beta_vec = 0:0.001:1; % effect size (we use absolute value so take only beta>0)

N = length(params_vec);
s_beta_vec = zeros(N, length(beta_vec));
switch lower(selection_func_str)
    case {'polynomial', 'power', 'monomial'} % here parameter k is power |z|^k
        param_str = 'k';
    case 'gaussian' % here parameter sigma is st.d. of Gaussian z ~ N(0,sigma)
        param_str = '\sigma';
end
for i=1:N
    switch selection_func_str
        case {'polynomial', 'power', 'monomial'} % here parameter k is power |z|^k
            s_beta_vec(i,:) = (1 - (2^(i/2) * gamma((i+1)/2) / sqrt(pi)) * ...
                hypergeom(-i/2, 1/2, -beta_vec.^2/2)); % ./ ...
            %        (1 - (2^(k/2) * gamma((k+1)/2) / sqrt(pi))); % selection for each beta
        case 'gaussian' % here parameter sigma is st.d. of Gaussian z ~ N(0,sigma)
            s_beta_vec(i,:) = exp(-beta_vec.^2 ./ (2.*(params_vec(i)^2+1)));
    end
end

%s_beta_vec = s_beta_vec .* 0.5; % Take maximal possible RAF!!!!
RAF = 1; % 1: selection for a carrier. 0.5 - assume RAF f=0.5 (maximum possible)
for log_scale = 0:1 % plot in linear and log scales
    full_figure(0);
    if(log_scale) % Plot 0.5*(s_beta_vec-1): assume RAF f=0.5 (maximum possible)
        semilogy(100*beta_vec, RAF.*(s_beta_vec-1), 'linewidth', 2); % plot the reduction in selection (deviation from one)
        %    plotyy(100*beta_vec, 0.5.*(s_beta_vec-1), 0, 0, @semilogy); %  'linewidth', 2);
        
    else
        plot(100*beta_vec, 100*RAF.*(s_beta_vec-1), 'linewidth', 2); % plot the reduction in selection (deviation from one)
    end
    
    %     for k=1:N
    %         if(log_scale)
    %             plot(100*beta_vec, log(-100*(s_beta_vec(k,:)-1)), color_vec(k), 'linewidth', 2); % plot the reduction in selection (deviation from one)
    %         else
    %             plot(100*beta_vec, 100*(s_beta_vec(k,:)-1), color_vec(k), 'linewidth', 2); % plot the reduction in selection (deviation from one)
    %         end
    %     end
    title('Mean selection coefficient as a function of effect size', 'fontweight', 'bold');
    xlabel('Effect Size (|\beta|) %', 'fontweight', 'bold', 'fontsize', 11);
    legend_vec = num2str_cell(num2cell(params_vec));
    for k=1:N
        legend_vec{k} = [param_str '=' legend_vec{k}];
    end
    if(log_scale)
        ylabel('Selection (s)', 'fontweight', 'bold', 'fontsize', 11);
        log_str = ''; % 'log';
        y_lim = [-10^(-2) -10^(-6)];
    else
        ylabel('Selection (s) %', 'fontweight', 'bold', 'fontsize', 11);
        log_str = '_linear';
        y_lim = [-1 0];
    end
    ylim(y_lim); xlim([0 20]); legend(legend_vec, 3);
    y_ticks = get(gca, 'YTick'); % Get right axis ticks % [-10^(-2) -10^(-6)]; % get(gca, 'YTick');    
    y_tick_labels = get(gca, 'YTickLabel');
    if(log_scale)
        y_tick_labels = [repmat('-10^{', length(y_tick_labels), 1) ...
            y_tick_labels repmat('}', length(y_tick_labels), 1)];
    end
    axesPosition = get(gca,'Position');          %# Get the current axes position
    hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
        'Color','none',...           %#   ... with no background color
        'YLim',y_lim,...     %#   ... and a different scale
        'YAxisLocation','right',...  %#   ... located on the right
        'XTick',[],...               %#   ... with no x tick marks
        'YTick', [], ...
        'YTickLabel', [], ...
        'Box','off');                %#   ... and no surrounding box
    %    y_ticks = get(gca, 'YTick'); % Get right axis ticks % [-10^(-2) -10^(-6)]; % get(gca, 'YTick');
    if(log_scale)
        y_ticks = linspace(y_lim(1), y_lim(2), length(y_ticks));
    end
    my_yticklabels(hNewAxes, y_ticks,  mat2cell(y_tick_labels, ones(length(y_ticks),1)), ...
        'YAxisLocation', 'right');
    
    if(exist('fig_outfile', 'var') && (~isempty(fig_outfile)))
        my_saveas(gcf, [fig_outfile log_str], format_fig_vec);
    end
end


