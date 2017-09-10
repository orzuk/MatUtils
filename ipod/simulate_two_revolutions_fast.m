% Simulate two competing percolations - faster algorithm. 
% Each one is starting from some
% random vertices. We assume that a vertex can be either un-influenced,
% 'red' or 'blue', and this is affected by its neighborhood.
%
function [flip_t_avg t_mean t_std X X_order] = ...
    simulate_two_revolutions_fast(E_init, X_init_1, X_init_2, alpha, iters, lambda_1, lambda_2)


E = E_init;


% Check that the init vector is legal
if(max(X_init_1 + X_init_2) > 1)
    problems_with_input = 999999
    return;
end



t_vec = zeros(1, iters);  % hold time-to-revolution for each iteration


N = length(X_init_1); % number of vertices
N_on_1 = sum(X_init_1); % initial number of infected vertices
N_on_2 = sum(X_init_2); % initial number of infected vertices
N_on = N_on_1 + N_on_2; % Total on vertice s


lambda_vec_1 = ones(1,N); % an initial lambda vector with insignificant
lambda_vec_2 = ones(1,N); % an initial lambda vector with insignificant


% values
lambda_tot_1 = sum(lambda_vec_1); lambda_tot_2 = sum(lambda_vec_2);

tot_flip_t_vec = zeros(1,ceil(alpha*N)-N_on);
tot_flip_t_squared_vec = zeros(1,ceil(alpha*N)-N_on);

for iter=1:iters

    X_order = zeros(1,N);
    c_vec = zeros(1,N); % holds the total money invested on each
    iter_is=iter;
    E = E_init; % We change E by adding more and more edges

    ttt = cputime;
    iter_is = iter;
    t = 0; % time inside simulation
    delta_t = 0; % assume no-time since last adoption event

    X = X_init_1 + 2*X_init_2; % initialize adoption vector : 0 - not-infected, 1- 'red',  2 - 'blue'

    flip_t_vec = zeros(1,ceil(alpha*N)-0*N_on); % holds successive adoption times
    lambda_total_vec = zeros(1,ceil(alpha*N)-N_on); % Hold successive total 'flip rate'. Should increase at the beginning
    % and then at the end decrease, at least if the tv factor is negligible


    % Compute the initial lambda_vec
    lambda_vec_1 = f_influence(E, X_init_1, c_vec, 0, 0, lambda_1,  0, 0); % update lambda vec
    lambda_vec_2 = f_influence(E, X_init_2, c_vec, 0, 0, lambda_1,  0, 0); % update lambda vec

    lambda_vec_1(find(lambda_vec_1));
    lambda_vec_2(find(lambda_vec_2));

    X_off = find(X==0); X_on_1 = find(X == 1); X_on_2 = find(X == 2);
    lambda_total_1 = sum(lambda_vec_1(X_off));         lambda_total_2 = sum(lambda_vec_2(X_off));

    % Find the edges which are in the boundary
    cross_edges_1 = E(X_off, X_on_1);
    cross_edges_2 = E(X_off, X_on_2);
    lam_tot_1 = full(sum(sum(cross_edges_1)));
    lam_tot_2 = full(sum(sum(cross_edges_2)));


     % Deal with sparsity by ourselves
    [I_1 J_1 S_1] = find(cross_edges_1);     [I_2 J_2 S_2] = find(cross_edges_2);
    
    % Loop until alpha*N are infected by either one of the colors
    for i=N_on+1:ceil(alpha*N) % loop over each adoption event until revolution

        % choose the vertex we flip randomly according to the lambdas, and
        % update the time.

        % First choose which color to flip
        %X_off = find(X==0);


        % Now choose the color of the flip
        r = rand(1) * (lam_tot_1 + lam_tot_2);
        if(r < lambda_total_1)

            flipped_edge = ceil(r);
            
            flipped_ind = I_1(flipped_edge)
            
            bad_ind = find(J_1 == flipped_ind);
            

            X(flipped_ind)=1; % mark vertex as 'red'
            cross_edges_1(flipped_ind, X_on_1) = 0; cross_edges_2(flipped_ind, X_on_2) = 0;
            cross_edges_1(flipped_ind, X_off) = 1;

            X_on_1 = [X_on_1 flipped_ind];
            X_off = setdiff(X_off, flipped_ind);


        else

            flipped_edge = ceil(r - lam_tot_1);

            [I J S] = find(cross_edges_2);
            flipped_ind = J(flipped_edge)

            X(flipped_ind)=2; % mark vertex as 'red'
            cross_edges_2(flipped_ind, X_on_2) = 0; cross_edges_1(flipped_ind, X_on_1) = 0;
            cross_edges_2(flipped_ind, X_off) = 1;

            X_on_2 = [X_on_2 flipped_ind];
            X_off = setdiff(X_off, flipped_ind);


        end
        lam_tot_1 = full(sum(sum(cross_edges_1)));
        lam_tot_2 = full(sum(sum(cross_edges_2)));

        delta_t = exprnd(1/ (lam_tot_1 + lam_tot_2));

        X_order(flipped_ind) = i; % Record the ordering

        E(flipped_ind,flipped_ind)=X(flipped_ind); % Mark 'full' by a self-loop in the adjancy matrix

        t = t + delta_t;
        flip_t_vec(i)=delta_t; % record time of adoption event


    end % loop on adopters


    flip_t_vec = flip_t_vec(N_on+1:end);

    tot_flip_t_vec = tot_flip_t_vec + flip_t_vec;
    tot_flip_t_squared_vec = tot_flip_t_squared_vec + flip_t_vec .* flip_t_vec;
    t_vec(iter) = t;

    if(mod(iter, 100) == 0)
        sprintf('Finished %ld iterations', iter)
    end

    iter_time = cputime-ttt;
end % loop on iters
flip_t_avg = tot_flip_t_vec ./ iters;
flip_t_std = sqrt((tot_flip_t_squared_vec - (flip_t_avg .^ 2)) ./ iters);
%  plot([N_on+1:ceil(alpha*N)], flip_t_avg, flip_t_std);

% Plot the times histogram
% figure; hold on; hist(t_vec, 50); xlabel('t'); ylabel('Density'); title(['time distribution. # iters = ' num2str(iters)]);
flip_t_avg_cumsum = cumsum(flip_t_avg);
XMIN=0; XMAX = flip_t_avg_cumsum(end); YMIN=0; YMAX=alpha*1.1;
% figure; plot(cumsum(flip_t_avg), [1:length(flip_t_avg)]./N); AXIS([XMIN XMAX YMIN YMAX]);
% title(['time to influence vertices \alpha=' num2str(alpha)]); xlabel('time'); ylabel('fraction influenced');

% plot the resulting graph:
%figure;

% Make the self-loops also for the initial vertices
for i=1:length(X_init_1)
    E(i,i)=E(i,i) || X_init_1(i);
end


%%% Don't need to draw all the time ...
%%%[x, y, labels] = draw_dot(E); title('Infected (gray) at the end');

% plot(sort(t_vec))

% return mean and std
t_mean = mean(t_vec);
t_std = std(t_vec);



