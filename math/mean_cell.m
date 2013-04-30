% Mean of the elements in each member of a cell array
function M = mean_cell(c)

n = length(c);
M = zeros(n,1); 
for i=1:n
    M(i) = mean(c{i});
end
