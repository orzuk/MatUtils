% Compute noncentrality parameter n_{a,b} for a type-1-error and b
% type-2-error in chi-square test
function NCP = two_types_errors_to_non_centrality_parameter(a, b)

NCP = fminbnd(@(x) ( ncx2cdf( chi2inv(1 - a,1), 1, x) - b ).^2, ...
    0.1,  chi2inv(1 - a,1) * 10 );
