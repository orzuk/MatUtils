% Convert strings to numbers in a cell array 
function nc = str2num_cell(s)
nc = s;
for i=1:length(nc)
    if(isa(nc{i}, 'char'))
        nc{i} = str2num(s{i});
    end
end
