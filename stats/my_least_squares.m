% Find the best a and b such that the sum of squares of y - (ax+b) is minimized 
function [a, b] = my_least_squares(x, y)

a = ( sum(x.*y) - sum(x)*mean(y) ) / var(x);
b = mean(y) - a*mean(x);












