% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = Px(x, c, sigma)

F = (1 ./ (sqrt(2*pi) .* sigma)) .* ( exp(-(x-c).^2 ./ (2.*sigma.^2)  ) + exp( -(x+c).^2 ./ (2.*sigma .^2)  )  );

