% Calculate the time until a fraction alpha of the popultaion is infected
% Do it in the fastest way, here there's no money and tv. Only percolation
% from neighbors
function [flip_t_avg t_mean t_std]=calc_time_to_revolution_fast(E_init, X_init, alpha, iters, lambda)

E = E_init;


t_vec = zeros(1, iters);  % hold time-to-revolution for each iteration



N = length(X_init); % number of vertices
N_on = sum(X_init); % initial number of infected vertices

lambda_vec = ones(1,N); % an initial lambda vector with insignificant
% values
lambda_tot = sum(lambda_vec);

tot_flip_t_vec = zeros(1,ceil(alpha*N)-N_on);
tot_flip_t_squared_vec = zeros(1,ceil(alpha*N)-N_on);

for iter=1:iters

    iter_is=iter;
    E = E_init; % We change E by adding more and more edges

    ttt = cputime;
    iter_is = iter;
    t = 0; % time inside simulation
    delta_t = 0; % assume no-time since last adoption event

    X = X_init; % initialize adoption vector
    c_vec = zeros(1,N); % holds the total money invested on each
    % individual so far

    flip_t_vec = zeros(1,ceil(alpha*N)-0*N_on); % holds successive adoption times
    lambda_total_vec = zeros(1,ceil(alpha*N)-N_on); % Hold successive total 'flip rate'. Should increase at the beginning
    % and then at the end decrease, at least if the tv factor is negligible



    % Compute the initial lambda_vec
    lambda_vec = f_influence(E, X, c_vec, 0, 0, lambda, ...
        0, 0); % update lambda vec

    % Loop until alpha*N are infected
    for i=N_on+1:ceil(alpha*N) % loop over each adoption event until revolution


        % choose the vertex we flip randomly according to the lambdas, and
        % update the time.

        [flipped_ind delta_t lambda_total]= randomize_vertex_chosen(X,lambda_vec);

        X(flipped_ind)=1; % mark vertex as adopter
        E(flipped_ind,flipped_ind)=1; % Mark 'full' by a self-loop in the adjancy matrix


        %      sprintf('%d: %f invested in flipping %d\n', i, c_vec(flipped_ind), flipped_ind)
        t = t + delta_t;
        flip_t_vec(i)=delta_t; % record time of adoption event
        lambda_total_vec(i) = lambda_total;


        % Now update lambda_vec locally by looking at the neighbors of the
        % flipped vertex
        flipped_neighbors = find(E(flipped_ind,:));
        flipped_neighbors_off = intersect(find(X==0),  flipped_neighbors );
        lambda_vec(flipped_neighbors_off) = lambda_vec(flipped_neighbors_off)+lambda;


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

for i=1:length(X_init)
    E(i,i)=E(i,i) || X_init(i);
end


%%% Don't need to draw all the time ...
%%%[x, y, labels] = draw_dot(E); title('Infected (gray) at the end');

% plot(sort(t_vec))

% return mean and std
t_mean = mean(t_vec);
t_std = std(t_vec);



