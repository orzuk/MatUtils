% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = Pxx(x, c, sigma)

F = (-1 ./ (sqrt(2*pi) .* sigma.^3)) .* ( (x-c) .* exp(-(x-c).^2 ./ (2.*sigma.^2)  ) + (x+c) .* exp( -(x+c).^2 ./ (2.*sigma .^2)  )  );

