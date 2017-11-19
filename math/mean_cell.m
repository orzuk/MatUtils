% Mean of the elements in each member of a cell array
function M = mean_cell(c)

[n, k] = size(c); %n = length(c);
M = zeros(n,k); 
for i=1:n
    for j=1:k
        M(i,j) = mean(c{i,j});
    end
end
