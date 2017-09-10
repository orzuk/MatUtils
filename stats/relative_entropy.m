% Compute the relative entropy between distributions P and Q
function KL = relative_entropy(P, Q)

KL = sum(P .* log(P)) - sum(P .* log(Q)); 