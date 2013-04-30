% Compute cdf(t) for a histogram
function c = cdf_hist(x, p, t)

i = find(x < t, 1, 'last'); 
if(~isempty(i))
    c = sum(p(1:i)) / sum(p); 
else
    c = 0;
end

