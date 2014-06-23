% Convert cell entries to doubles
function d = double_cell(c)
n = length(c); d = cell(n,1);
for i=1:n
    d{i} = double(c{i});
end
    
    