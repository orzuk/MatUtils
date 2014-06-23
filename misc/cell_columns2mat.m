% Put together a few cell array columns in one cell array matrix
function c = cell_columns2mat(c_cols)

m = length(c_cols); c_lens = length_cell(c_cols); n = max(c_lens); 
c = cell(n,m); 
for i=1:m
    c(1:c_lens(i),i) = c_cols{i};
end

