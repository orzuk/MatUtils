% Write the cell array as a string with the correct matlab syntax
function s = cell2str(c)
n  = length(c);

for i=1:n
    c{i} = ['''' c{i} ''''];
end

s = ['{' cell2vec(c, ', ') '}'];
% s = [s(1:end)  '}']; % remove last ','

