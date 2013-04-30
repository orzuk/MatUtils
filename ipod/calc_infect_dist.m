% A script for studying a model of two competing infections 
% and determining fixation distribution.
% We calculate the distribution of red infected vertices
% assuming that we start with two adjacent on

xglobal ER LATTICE BA; competeing

graph_ER = 0; % Erdos-Renyi
graph_LATTICE=1; % Z^d lattice
graph_BA=2; % Barabasi-Alberts preferential attachment model
graph_SO=3; graph_WS=4;
graph_DD=5;  % Given degree distribution

d=30;
N=d^2; % Number of vertices

ttt=cputime;
graphs_num = 1; % 10 graphs for each p
iters=2500; % Number of iterations for each graph
alpha=0.49999; % Fraction of occupied vertices

k_init = 1000*N; % total amount of money at the beginning
lambda_0 = 0; % very low. This is the flip rate when there is no positive influence
lambda_1 = 1; lambda_2 = 1*lambda_1; % 0.8; % The influence on flip rate of the strategies
money_factor =  0.00000000002*N;  % What the hell is this?
tv_factor = 0.00000000002*N; % What the hell is this? make it negligible

p_res = 0.5;

p_vec = p_res:p_res:1-p_res;
t_mean_p_vec = zeros(1,length(p_vec));

graph_type = graph_LATTICE; % Select which graph to work with

i=1;
for p=p_vec % Loop over p-values
    do_p = p
    for g=1:graphs_num  % loop over ER graphs with a given p


        % Generate a random graph
        if(graph_type == graph_ER)
            E_ER=generate_random_graph(graph_ER, N, 4/N); % Hope it is connected ...
            graph_str = 'Erdos-Renyi';
        else
            E_ER=generate_random_graph(graph_LATTICE, N, 2); % Do a 2-d lattice for comparison. This is like p=4*d/(N-1)
            graph_str = 'Z_2 lattice';
        end

        X_init = zeros(1,N);

        X_init_1 = X_init;X_init_1(round((N-d)/2))=1;
        X_init_2 = zeros(1,N);%X_init_2(round((N-d)/2))=1;
        %    X_init_2(1:d)=1; X_init_2(d:d:d*d)=1;
        X_init_2(round((N-d)/2)+1)=1;

        % Do simulations
        %         [flip_t_avg t_mean t_std] = calc_time_to_revolution(E_ER, X_init, alpha, iters, ...
        % 					 k_init, lambda_0, lambda_1, ...
        % 					 money_factor, tv_factor);

        [flip_t_avg t_mean t_std X_two_col X_order W_flip] = simulate_two_revolutions(E_ER, X_init_1, X_init_2, alpha, iters, lambda_1, lambda_2);

        %[flip_t_avg t_mean t_std]=calc_time_to_revolution_fast(E_ER, X_init, alpha, iters, lambda_1);

        t_mean_p_vec(i) = t_mean_p_vec(i)+t_mean;
    end
    i=i+1;
end

t_mean_p_vec = t_mean_p_vec ./ graphs_num; % correct to give the average

time_elapsed = cputime - ttt

% Plot results
figure; imagesc(reshape(X_two_col, sqrt(N), sqrt(N))); colorbar; title('points occupied by two colors');

figure; imagesc(reshape(X_order, sqrt(N), sqrt(N))); colorbar;  title('Ordering of points traversed');

red_frac = sum(X_two_col == 1) / length(X_two_col)
blue_frac = sum(X_two_col == 2) / length(X_two_col)

%figure; hold on; plot(p_vec(1:end), t_mean_p_vec(1:end), '*'); title('rev. time for different ps'); xlabel('p'); ylabel('mean time');
%Y = (rand(d) > 0.5)+1; Y(1)=0; figure; imagesc(Y); colorbar;





% Plot distributions
num_flips = size(W_flip,2)
W_flip = cumsum(W_flip, 2) ./ repmat([1:num_flips],  size(W_flip,1), 1);
W_flip_mean = mean(W_flip); W_flip_std = std(W_flip);

figure; subplot(3,4,1); plot([1:num_flips], W_flip_mean); title('mean red frac');
subplot(3,4,5); plot([1:num_flips], W_flip_std); title('std red frac');
subplot(3,4,9); errorbar([1:num_flips], W_flip_mean, W_flip_std);  title(graph_str);

for i=[2:4,6:8]
    subplot(3,4,i); plot([1:num_flips], W_flip(i,:)); title('typical red frac');
end

N_vec = round(num_flips .* [1:3] ./ 3);
for i=10:12
    subplot(3,4,i); hist(W_flip(:,N_vec(i-9)), 150); title(['red frac dist. N = ' num2str(N_vec(i-9))]);
end










