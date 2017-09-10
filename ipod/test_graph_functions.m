% Test all functions
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, '/ph2users/fezuk\pacperm');
path(path, '/ph2users/fezuk\pacperm\numeric');
path(path, 'E:\Research\ipod\graph');



% types of graphs
global ER LATTICE BA;

graph_ER = 0; graph_LATTICE=1; graph_BA=2; graph_SO=3; graph_WS=4


% Generate a random graph
p=0.5;d=4; N=2^d;
E_ER=generate_random_graph(graph_ER, N, p);
E_grid=generate_random_graph(graph_LATTICE, N, d);


iters=1000; 
alpha=0.8; 
X_init = zeros(1,N);X_init([1])=1; 
k_init = 1000*N; % total amount of money at the beginning
lambda_0 = 0; % very low. This is the flip rate when there is no positive influence
lambda_1 = 1; % 0.8; % The influence on flip rate of the strategies
money_factor =  0.00000000002*N;  % What the hell is this?
tv_factor = 0.00000000002*N; % What the hell is this? make it negligible

do_several_tests=0;


% Try the backtrack program, X_init is not needed here
 t_mean_vec =calc_all_mean_times(E_ER, X_init, alpha, lambda_1);



[flip_t_avg t_mean t_std] = calc_time_to_revolution(E_ER, X_init, alpha, iters, ...
					 k_init, lambda_0, lambda_1, ...
					 money_factor, tv_factor)

figure; hold on; plot(hamming_weight(1:2^N-1), t_mean_vec(2:end),  '*'); title('Time to rev for different set sizes');
 
flip_t_avg_vec = cumsum(flip_t_avg(end:-1:1));         
flip_t_avg_vec = flip_t_avg_vec(end:-1:1);                  
plot( flip_t_avg_vec, 'r'); plot( flip_t_avg_vec, '*r');                    
legend('Exact', 'Simulated'); 


                % do_several_tests takes longer... do
				% them only if specified by
				% do_several_tests variable
if ((~exist('do_several_tests'))  || (do_several_tests == 0))
  break
end

k_iter_vec =   0:0.1*N:0.5*N        %   0.99 % 0.39:0.3:0.99;
t_mean_vec = zeros(1,length(k_iter_vec)); t_std_vec = t_mean_vec;

ii=1; 
flip_t_avg= zeros(length(k_iter_vec), ceil(alpha*N));
for k_iter=k_iter_vec
  iters=25; 
  k_init=k_iter; % Increase the amount of people need to be influenced
  X_init = zeros(1,N); 
  %k_init = 10*N; % total amount of money
  lambda_0 = 0.01; % very low
  lambda_1 = 1;  
  money_factor = 0.2*N;
  tv_factor = 0.2*N;

  [flip_t_avg(ii,:)  t_mean_vec(ii) t_std_vec(ii)] = calc_time_to_revolution(E_grid, X_init, alpha, iters, ...
					   k_init, lambda_0, lambda_1, ...
					   money_factor, tv_factor);
  ii=ii+1;
end

color_vec = 'brgkcm';

figure; hold on; XMAX=0;
for ii=1:length(k_iter_vec)
    plot(cumsum(flip_t_avg(ii,:)), [1:length(flip_t_avg)]./N, color_vec(ii)); 
    XMAX = max(XMAX, max(cumsum(flip_t_avg(ii,:))));
end
XMIN=0; YMIN=0; YMAX=alpha*1.1;      
AXIS([XMIN XMAX YMIN YMAX]); 
title(['time to influence vertices \alpha=' num2str(alpha)]); xlabel('time'); ylabel('fraction influenced');


figure; hold on; errorbar(k_iter_vec, t_mean_vec, t_std_vec); 
title('Time needed to reach influence'); xlabel('\alpha'); ylabel('t');

