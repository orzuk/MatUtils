% A Gaussin integral
function Z = exp(-(r+alpha).*t) 


function r = M_expectation(alpha,t) 

r = quad2d(@(x,y)M_integrand(x,y,alpha,t),

% The expectation under normals 
function M_integrand(x,y,alpha,t) = 
    normpdf(x).*normpdf(y).* ...
        normcdf(exp(-alpha.*t) ./ sqrt(0.5*(1-exp(-2*alpha*t)))) .* ...
        (1-normcdf(x)) ./ normpdf(x);


% Bar riddle: 
x = rand(10000000,1) > 0.5; k = 4; 
change_points = find(diff(x));
x_lens = diff(change_points); 
bad_bars = find(x_lens(1:2:end) >= k); % look only at blacks
consecutive_inds = find(diff(bad_bars) == 1);
mean(x_lens(bad_bars(consecutive_inds)+1))
mean(x_lens(bad_bars+1))
