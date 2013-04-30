% Find the (unique) minimax solution for the game of infection when each one start
% chooses a vertex, and tries to maximize its infection at the end. The
% output is the vector of probabilities saying with what probability to
% start with which vertex
function p_vec = find_revs_game_minimax_solution(E, game_flag)

N = length(E) % number of verticess


% calculate the game matrix
P_mat = calc_two_revs_game_matrix(E, game_flag) 


trans_flag = is_transitive(P_mat)

% In the case of transitive matrix, we can order the vertices
if(trans_flag == 1)
    ordered_vertices_strategies =sum( max(sign(P_mat), 0));
    [val ordered_vertices_strategies] = sort(-ordered_vertices_strategies);
    ordered_vertices_strategies
end

% first check for pure strategy
[min_max i_min] = min(max(P_mat))
[max_min j_min] = max(min(P_mat))

if(min_max == max_min)
    % Print the minima solution
    game_value_is = min_max
    pure_strategies_are = j_min
    
    p_vec = zeros(1,N); p_vec(j_min)=1;
else % Here we must find the mixed solution using linear programming
    
    % Prepare f, A and b
    A = zeros(2*N,N+1);
    A(1:N,1:N) = P_mat;
    A(1:N,N+1)=-1;
    A(N+1:2*N,1:N) = -eye(N);
    
    f = zeros(N+1,1); f(end)=1;
    
    b = zeros(2*N,1);
    
    Aeq = ones(1,N+1); Aeq(end) = 0; beq = 1;
    
   [p_vec, game_val]= linprog(f,A,b, Aeq, beq)
end


