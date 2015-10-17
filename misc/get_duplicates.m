% Like matlab's unique, but gives all values' indices and multiplicities
%
% Input: 
% v - a vector
% 
% Output: 
% vals - unique values 
% inds - indices of unique values
% num_dups - number of times each value appears
% 
function [vals inds num_dups] = get_duplicates(v)

[v sort_perm] = sort(v);
[vals I J] = unique(v);

num_unique = length(vals);
num_dups = hist(J,1:num_unique);

ctr=1; inds = cell(num_unique,1);
for i=1:num_unique
    inds{i} = sort_perm(ctr:I(i));
    ctr = I(i)+1;
end

