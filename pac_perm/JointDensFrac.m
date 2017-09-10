% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = JointDensFrac(y, sig, C_alpha, one_side_flag)

if(one_side_flag)
     F = (1/sqrt(2.*pi)) .* exp(-y.*y./2) .* (1 - normcdf( (C_alpha.* sqrt(1+sig.^2)-y)./sig ));
else
     F = (2/sqrt(2.*pi)) .* exp(-y.*y./2) .* (1 - normcdf( (C_alpha.* sqrt(1+sig.^2)-y)./sig ) + normcdf( (-C_alpha.* sqrt(1+sig.^2)-y)./sig ));
end