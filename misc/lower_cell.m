% Transfer cell content to lower-case
function l = lower_cell(c)
n = length(c); l = cell(n,1);
for i=1:n
    l{i} = lower(c{i});
end


