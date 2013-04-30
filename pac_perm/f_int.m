% Perform an integral by deviding into small areas to increase accuracy
function F = f_int(c, x, sigma)

F = PInt(x,c,sigma) .* exp(-c.^2./2) ./ sqrt(2*pi);



