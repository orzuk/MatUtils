% Inverse standard Gaussian Phi
function x = PhiInv(p)

x = sqrt(2)*erfinv(2*p-1);
