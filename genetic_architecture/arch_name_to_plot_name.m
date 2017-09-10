% Make architecture name nice for plotting
function plot_name_arch = arch_name_to_plot_name(arch)

plot_name_arch = arch;
fi = strfind(plot_name_arch, '(');
if(~isempty(fi))
    plot_name_arch = plot_name_arch(1:fi-1);
end
if(~isempty(strfind(arch, '-of-sig')) && isempty(strfind(arch, 'moids')))
    plot_name_arch = [plot_name_arch 'moids']; % add 'sigmoid'
end
