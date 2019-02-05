% Internal function for setting defaults: density un-weighted
function plot_params  = set_default_plot_params(plot_params)
if(~isfield(plot_params, 'new_fig')) % normalize distribution
    plot_params.new_fig = 1;
end
if(~isfield(plot_params, 'ylabel_str'))
    plot_params.ylabel_str = '\Psi';
end
if(~isfield(plot_params, 'normalize')) % normalize distribution
    plot_params.normalize = 0;
end
if(~isfield(plot_params, 'cum')) % cumulative
    plot_params.cum = 0;
end
if(~isfield(plot_params, 'log')) % default: plot semilogx
    plot_params.log = [1 0];
end
if(isscalar(plot_params.log))
    plot_params.log = [plot_params.log plot_params.log];
end
if(~isfield(plot_params, 'weighted')) % plot log-log
    plot_params.weighted = 0;
end
if(~isfield(plot_params, 'xlim')) % plot lim
    plot_params.xlim = [10^(-4) 1];
end
if(~isfield(plot_params, 'hist')) % hold distribution as histgoram
    plot_params.hist = 0;
end
if(~isfield(plot_params, 'font_size')) % set default font size 
    plot_params.font_size = 14;
end
if(~isfield(plot_params, 'color_vec'))  % New: Set also color vec !
    AssignRVASConstants;
    plot_params.color_vec = num2cell(color_vec); 
end
