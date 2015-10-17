% Transfer cell content to upper-case
function l = upper_cell(c)
n = length(c); l = cell(n,1);
for i=1:n
    l{i} = upper(c{i});
end
