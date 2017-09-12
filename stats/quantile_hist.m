% Compute quantile for a histogram
function v = quantile_hist(x, p, q)

[dummy med_ind] = min(abs(cumsum(p)./sum(p) - q));
v = x(med_ind);
