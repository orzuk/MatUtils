% The function tries to find the 'best' curve, i.e. the curve which is as
% far as possible from the points given in the input. The function
% currently works for quadratuc curves passing through the origin: Y=a*X^2+b*X
function grad = AvoidPointsGaussianFieldsGrad(ab,x, y, sigma) 
a = ab(1); b = ab(2); 

grad(1) = -sum(  (x.^2./sigma^2).*(y-a*x.^2-b.*x) .*  exp( -(y-a*x.^2-b.*x).^2./(2*sigma^2) ) );  % df/da
grad(2) = -sum(  (x./sigma^2).*(y-a*x.^2-b.*x) .*  exp( -(y-a*x.^2-b.*x).^2./(2*sigma^2) ) );  % df/db


%% The 3-d case: use tanh: Y = (a*X^2+b*X)*tanh(c*X)
%% c = ab(3); 
%% grad(1) = -sum(  (x.^2.*tanh(c.*x)./sigma^2).*(y-tanh(c.*x).*(a*x.^2+b.*x)) .*  ...
%%    exp( -(y- tanh(c.*x).*(a*x.^2+b.*x)).^2./(2*sigma^2) ) );  % df/da
%% grad(2) = -sum(  (x.*tanh(c.*x)./sigma^2).*(y-tanh(c.*x).*(a*x.^2+b.*x)) .*  ...
%%    exp( -(y- tanh(c.*x).*(a*x.^2+b.*x)).^2./(2*sigma^2) ) );  % df/db
%% grad(3) = -sum(  x.*(a.*x.^2+b.*x) .* (1-tanh(c.*x).^2) .*(y-tanh(c.*x).*(a*x.^2+b.*x)) .*  ...
%%    exp( -(y- tanh(c.*x).*(a*x.^2+b.*x)).^2./(2*sigma^2) ) );  % df/dc


