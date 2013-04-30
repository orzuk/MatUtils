% Perofrm de-convolution of functions using the Discrete-Fourier-Transform
sigma = 0.0010;  g_sigma = sqrt(1+sigma^2); num_points = 900;
% n=[ones(1,1009)];
% x = rand(1,100000);
% z = sigma .* randn(1,100000);
% y=x+z;  % Noisy observation
% [g, bins] = hist(y, 100); g = g ./ 2000;
t=linspace(-4,4,num_points); t2 = linspace(-1,1,1000); bin_size = t(2)-t(1);

correct_q=exp(-t.^2/2) ./ sqrt(2*pi); 
n=exp(-t.^2/(2.*sigma.^2)) ./ sqrt(2*pi*sigma); 
%n =ones(1,num_points) ./ 8;
g=exp(-t.^2/(2.*g_sigma.^2)) ./ sqrt(2*pi*g_sigma); 

N = fftshift(fft(n,1028));
G = fftshift(fft(g,1028));

% remove zeros 
N1 = N(2:end);
G1 = G(2:end);
Q1=[0 G1./(N1+eps)];
Q=G./(N+eps);

q=ifft(fftshift(Q1));
%figure; plot(q)

% New : Try to solve linear equations : 
b = g';
A = zeros(num_points);
for i=1:num_points
    A(i,:) = exp( -(t(i)-t).^2./(2.*sigma.^2) );
end
A = A ./ (sqrt(2*pi)*sigma);

qq = linsolve(A,b) ./ bin_size;

figure; hold on; plot(t, g, 'r');  plot(t, qq);  legend('noisy g', 'reconstructed q');
figure; hold on;  plot(t, qq); % legend('noisy g', 'reconstructed q');