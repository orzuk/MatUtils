% Compute median for a histogram
function m = median_hist(x, p)

m = quantile_hist(x, p, 0.5);


