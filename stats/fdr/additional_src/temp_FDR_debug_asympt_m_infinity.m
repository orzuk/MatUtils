% Temp debuggind script for FDR 
p = 0.399;
m = 100000; 
q = 0.7;
r = rand(m,1); 
k = floor((1-p)*m);
r(1:k) = 1 - sqrt(1-r(1:k));  % The non-null density: F(x) = 2x-x^2
%figure; hist(r, 500);
x_vec = [1:m]./m;
figure; hold on; plot(x_vec, sort(r), '.'); plot(x_vec, q.*x_vec, 'r');
plot(x_vec, p.*x_vec + (1-p).*(1 - sqrt(1-x_vec)), 'g');

figure; hold on; plot(x_vec, sort(r(1:k)), '.'); 
x_vec = [1:k]./k;
plot(x_vec, 1 - sqrt(1-x_vec), 'g');




r = sort(r);
R = FDR_main(r,q,'bh95'); 
R_emp = size(R,1) / m
x_emp = r(size(R,1))

x = (2*q-q*p-1) / (q^2*(1-p)) % The cutoff from the equation: x/F(qx) = (1-p)/(1-pq)

FDR_theo = p*x / (p*x + (1-p) * (2*x-x^2))
% R_theo = (p*x + (1-p) * (2*x-x^2))
R_theo = (p*x + (1-p) * (1-sqrt(1-x)))

x_vec = [1:m] ./ m;
figure; hold on; plot(x_vec, r, '.'); plot(x_vec, q.*x_vec, 'r');
plot(x_vec, p.*x_vec + (1-p).*(1 - sqrt(1-x_vec)), 'g');


plot(R_emp, x_emp, '+m');
plot(x, R_theo, 'om');

