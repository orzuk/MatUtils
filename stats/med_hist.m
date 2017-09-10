% Compute median for a histogram
function v = med_hist(x, p)

[dummy med_ind] = min(abs(cumsum(p)./sum(p) - 0.5));
v = x(med_ind);

