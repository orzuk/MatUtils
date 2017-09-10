% Creat a vector of values from a histgoram 
% 
% Input: 
% vals - unique values
% counts - how many time each value appears (assumed integers)
% 
% Output: 
% v - vector such that each value vals(i) appears count(i) times
% 
function v = hist_to_vals(vals, counts)

v = zeros(1, sum(counts)); ctr=1;
for i=1:length(vals)
    v(ctr:(ctr+counts(i)-1)) = vals(i); ctr=ctr+counts(i); 
end

