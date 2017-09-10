% Split a vector into a cell array of vectors 
% 
% Input: 
% v - one vector
% lens - the cumulative lengths (start points) of the vectors in c (so we know where to start from each element)
%
% Output: 
% c - a cell array with many vectors 
function c = vec2cell(v, lens)

n = length(lens); % get number of segments  
c = cell(n,1); 
lens = [0 vec2row(lens)]; % add a zero at the beggining 
for i=1:n
    c{i} = v(lens(i)+1:lens(i+1)); % use cumulative sum so we know where to start each vector
end  

