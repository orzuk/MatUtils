% Temp function for debug
function r = phi_s(x,S) 

r = (1-exp(-S*(1-x))) / (x.*(1-x).*(1-exp(-S)));
