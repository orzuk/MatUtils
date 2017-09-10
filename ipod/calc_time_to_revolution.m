% Calculate the time until a fraction alpha of the popultaion is infected
% E is the graph (as a sparse matrix)
% X is the initial set of active ('on') vertices. Usually assume that all
% vertices are at the beginning 'inactive'
% alpha is the desired fraction of active vertices
% iters is the number of iteration
% k_init is the total money at the beginning
% lambda_0 is the 'basal' flip rate and lambda_1 is the rate reulsting from
% neighbors, money, global field etc.
function [flip_t_avg t_mean t_std]=calc_time_to_revolution(E_init, X_init, alpha, iters, ...
    k_init, lambda_0, lambda_1, ...
    money_factor, ...
    tv_factor)

E = E_init;
%figure; graph_draw(E); title('Infected (gray) at the beginning'); 



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
    k_global = 0; % initially no publicity investmenet
    k_left = k_init; % initialize budget
    X = X_init; % initialize adoption vector
    c_vec = zeros(1,N); % holds the total money invested on each
    % individual so far

    flip_t_vec = zeros(1,ceil(alpha*N)-0*N_on); % holds successive adoption times
    lambda_total_vec = zeros(1,ceil(alpha*N)-N_on); % Hold successive total 'flip rate'. Should increase at the beginning
                                                    % and then at the end decrease, at least if the tv factor is negligible
    
    
    % Loop until alpha*N are infected
    for i=N_on+1:ceil(alpha*N) % loop over each adoption event until revolution
       
        
        
        % update the strategy. k_vec are the local moneys, k_global the
        % global 'tv' money, k_left is the total
        % money left, which is updated
        if (delta_t * k_global > k_left) % ran out of money
            k_left = 0;
            delta_t = k_left / k_global; % k_global must be positive because
            % of if condition
        else
            k_left = k_left - k_global * delta_t;
        end
        
       
        [k_vec k_global] = update_strategy(E,X,alpha,k_left, lambda_vec, lambda_tot);
       
        
        
        k_left = k_left - sum(k_vec);
        c_vec = c_vec + k_vec;
        
        lambda_vec = f_influence(E, X, c_vec, k_global, lambda_0, lambda_1, ...
            money_factor, tv_factor); % update lambda vec
        % choose the vertex we flip randomly according to the lambdas, and
        % update the time. we use the fact that
        % a sum of exponential random variables
        % is exponential with parameter equal to
        % sum of individual parameter
       
%         lamlam = lambda_vec
        
        [flipped_ind delta_t lambda_total]= randomize_vertex_chosen(X,lambda_vec);
        
        %lambda_total
        
        X(flipped_ind)=1; % mark vertex as adopter
        E(flipped_ind,flipped_ind)=1; % Mark 'full' by a self-loop in the adjancy matrix


        %      sprintf('%d: %f invested in flipping %d\n', i, c_vec(flipped_ind), flipped_ind)
        t = t + delta_t;
        flip_t_vec(i)=delta_t; % record time of adoption event
        lambda_total_vec(i) = lambda_total;
        
       
    end % loop on adopters
    
    
    flip_t_vec = flip_t_vec(N_on+1:end);
%     size(tot_flip_t_vec)
%     size(flip_t_vec)
    
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
[x, y, labels] = draw_dot(E); title('Infected (gray) at the end'); 

% plot(sort(t_vec))

% return mean and std
t_mean = mean(t_vec);
t_std = std(t_vec);



