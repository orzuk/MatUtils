% Like str2num for an array and deal with empty strings (from Tal Shay)
function N = my_str2num(S)

[n m] = size(S);
N = nan*ones(size(S));
for i = 1:n
    for j = 1:m
        if (~isempty(S{i,j}))
            t = str2num(S{i,j});
            if (~isempty(t))
                N(i,j) = t; 
            end
        end
    end
end