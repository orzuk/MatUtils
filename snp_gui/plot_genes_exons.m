%function [fig_ind, xmin, xmax, ymin, ymax, legend_vec, legend_cell]= plot_genes_exons(chr_num, loc_start, loc_end, genes_chr, gene_symbols, genes_loc_start, exon_start_cell, exon_end_cell)
function [fig_ind, xmin, xmax, ymin, ymax, legend_vec, legend_cell]= plot_genes_exons(chr_num, loc_start, loc_end, genes_chr, gene_symbols, genes_loc_start, exon_start_cell, exon_end_cell)

chr_ind = find(genes_chr == chr_num);
loc_ind = find(genes_loc_start(chr_ind) >= loc_start & genes_loc_start(chr_ind) <= loc_end);
ret_ind = chr_ind(loc_ind);
fig_ind = figure;
xmin = loc_start;
xmax = loc_end;
ymin = 0;
ymax = 1;
axis([xmin xmax ymin ymax]);
hold on;
% plot each gene
num_plot_genes = length(ret_ind);
legend_vec = [];
legend_cell = {};
legend_ind = 1;
exon_legend_flag = 0;

for a = 1:num_plot_genes
    gene_symbol = char(gene_symbols{ret_ind(a)});
    exon_start_vec = exon_start_cell{ret_ind(a)};
    exon_end_vec = exon_end_cell{ret_ind(a)};
    num_exons = length(exon_end_vec);
    for b = 1:num_exons
        ind_start = plot([exon_start_vec(b) exon_start_vec(b)] ,[ymin ymax] ,'--y');
        ind_end = plot([exon_end_vec(b) exon_end_vec(b)] ,[ymin ymax] ,'--g');
        if(~exon_legend_flag)
            legend_vec(legend_ind) = ind_start;
            legend_cell{legend_ind} = 'Exon Start';
            legend_ind = legend_ind+1;
            legend_vec(legend_ind) = ind_end;
            legend_cell{legend_ind} = 'Exon End';
            legend_ind = legend_ind+1;
            exon_legend_flag = 1;
        end
    end
end

y_loc = -0.01;
plot_x_labels_vertical(genes_loc_start(ret_ind), y_loc, gene_symbols(ret_ind));
