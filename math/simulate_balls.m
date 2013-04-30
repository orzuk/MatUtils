% Add one ball to the pile at each time with prob. proportional to existing pile
function simulate_balls(n)

x = rand(n,1); 
y = zeros(n,1); y(1) = 1; 
white = 1;
for i=3:n
    y(i) = x(i) < (white / (i-1));
    white = white + y(i); 
end

frac_white = cumsum(y) ./ (1:n)'; 

figure; plot(frac_white, '.'); 
