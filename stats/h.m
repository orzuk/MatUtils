% Binary entropy function 
% ENTROPY Entropy log base 2
% H = entropy(v)
function H = h(v)

H = entropy([v 1-v]');
