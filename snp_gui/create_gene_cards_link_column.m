%function link_str_column = create_gene_cards_link_column(gene_list)
function link_str_column = create_gene_cards_link_column(gene_list)


num_links = size(gene_list, 1);

link_str_column = cell(num_links, 1);
link_str_column(:) = {'=HYPERLINK("[http://www.genecards.org/cgi-bin/carddisp.pl?gene='};
link_str_column = concat_cell_str_column(link_str_column, gene_list);
temp_column = cell(num_links, 1);
temp_column(:) = {']", "Gene Cards")'};

link_str_column = concat_cell_str_column(link_str_column, temp_column);
