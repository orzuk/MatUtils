% Check if cell entries are empty 
function v = isempty_cell(c)
n = length(c); v = zeros(n,1);
for i=1:n
   v(i) = isempty(c{i}); 
end
