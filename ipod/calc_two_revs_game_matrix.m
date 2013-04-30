% Calculate the matrix for the zero sum game between two players who start
% infections, where each players profit is the number of sites he managed
% to infect
% The game_flag indicates which kind of game is it. 0 is infection
% by first-passage percolation, and 1 is infection by random walk.
function P_mat = calc_two_revs_game_matrix(E, game_flag)

N = length(E) % number of verticess

P_mat = zeros(N); % The payment matrix



if(game_flag == 0) % Here do complicated First-Passage-Percolation

    % Now start the recursion formulas for computing the profit for any initial
    % state of the two players. Note that the state space here is huge: ~3^N. In
    % general, for k players, the state space is roughly of size (k+1)^N
    % Each state represents a division of the vertex set to 3 sets. A is the
    % set of vertices with first infection (denoted 1), B the set of vertices
    % with second infection (denoted 2), and V \ {A U B} is the set of yet
    % unoccipied vertices (denoted 0).
    is_occupied_vec = zeros(1,3^N); % 0 means free, 1 means occupied
    color_occupied_vec = zeros(1,3^N); % 0 means A, 1 means B

    phi_edges_vec = zeros(1,2^N); % For each state store the edges with one vertex in the set and one in its complement

    ham_vec = hamming_weight([0:2^N-1]);  % Hold the "hamming vec", i.e. number of infected for both subsets

    starts_ind_vec = zeros(1,2^N);

    % First build the two vectors
    ind = 1;
    for i=0:2^N-1
        is_occupied_vec(ind:ind+2^ham_vec(i+1)-1) = i;  % Just store the indexes of the occupied vertices

        % Now comes the more difficult part of saying which are in A
        color_occupied_vec(ind:ind+2^ham_vec(i+1)-1) = stretch_to_indexes([0:2^ham_vec(i+1)-1], i);

        starts_ind_vec(i+1)=ind;

        ind = ind + 2^ham_vec(i+1);
    end

    ham_vec3 = hamming_weight(is_occupied_vec); % The hamming of the true states


    % The payment (for A) vec, i.e. expected #vertices of A
    pay_vec = zeros(1,3^N);

    % Start with the 'initial condition', which is the payment at the end
    end_indexes = find(ham_vec3 == N);
    pay_vec(end_indexes) = hamming_weight(color_occupied_vec(end_indexes));

    % Now start the recursion. There's no alpha here. We assume that we wait
    % until all vertices are occupied
    for num_on=N-1:-1:2  % loop over the number of on vertices in BACKWARDS order

        cur_states_indexes = find(ham_vec == num_on); % Find the relevant states for this stage


        cur_states_start_indexes = starts_ind_vec(cur_states_indexes); % Indicate the starting of the indexes in the big 3-array

        % loop over all possible neighbors, assume for now that N <=32
        for i=1:N



            % Find when to add another vertex
            add_i_states_indexes = bitor(cur_states_indexes-1, 2^(i-1))+1; % The indexes of the sets with one more vertex
            add_i_states_start_indexes = starts_ind_vec(add_i_states_indexes); % Indicate the starting of the indexes in the big 3-array

            % Now we need to compute the degree-into the set {A U B}
            i_neighbors = sum(2.^(find(E(:,i))-1)); % This just depends on the graph and stays the same
            d_i_indexes = hamming_weight(bitand(i_neighbors, cur_states_indexes-1)) .* (1- mod( floor((cur_states_indexes-1) ./ 2^(i-1)), 2));

            % How do we seperate between edges coming from A and edges coming from B ?


            % Loop on the entire 'color-vec' index set - very slowly ... :
            for AB_ind = 1:length(cur_states_indexes)

                cur_colors = color_occupied_vec(add_i_states_start_indexes(AB_ind):add_i_states_start_indexes(AB_ind)+2^(num_on+1)-1);

                % We need to sort based on the i-th place
                I_bit = bitand(color_occupied_vec(add_i_states_start_indexes(AB_ind):add_i_states_start_indexes(AB_ind)+2^(num_on+1)-1), 2^(i-1));
                vec_to_sort = color_occupied_vec(add_i_states_start_indexes(AB_ind):add_i_states_start_indexes(AB_ind)+2^(num_on+1)-1);
                I_bit = bitand(vec_to_sort, 2^(i-1));
                vec_to_sort = bitand(vec_to_sort, bitcmp(2^(i-1),32)) + bitshift(I_bit, 30-i) ;
                [ val ind] = sort( vec_to_sort);


                A_on = stretch_to_indexes([0:2^num_on-1], cur_states_indexes(AB_ind)-1);
                B_on = strech_to_indexes(2^num_on-1-[0:2^num_on-1], cur_states_indexes(AB_ind)-1);


                % count number of neighbors of vertex i in A and in B
                d_A_i_indexes = hamming_weight(bitand(i_neighbors,  A_on)) .* ...
                    (1- mod( floor((cur_states_indexes(AB_ind)-1) ./ 2^(i-1)), 2));
                d_B_i_indexes = hamming_weight(bitand(i_neighbors, B_on)) .* ...
                    (1- mod( floor((cur_states_indexes(AB_ind)-1) ./ 2^(i-1)), 2));

                pay_vec(cur_states_start_indexes(AB_ind):cur_states_start_indexes(AB_ind)+2^num_on-1) = ...
                    pay_vec(cur_states_start_indexes(AB_ind):cur_states_start_indexes(AB_ind)+2^num_on-1) + ...
                    d_B_i_indexes .* pay_vec(add_i_states_start_indexes(AB_ind) + ind(1:2^num_on)-1) + ...
                    d_A_i_indexes .* pay_vec(add_i_states_start_indexes(AB_ind) + ind(2^num_on+1:2^(num_on+1))-1); %  val(2^num_on+1:2^(num_on+1));

                % Pick one of the colorings, where all is colored by A
                phi_edges_vec(cur_states_indexes(AB_ind)) = phi_edges_vec(cur_states_indexes(AB_ind)) + d_B_i_indexes(1);

            end % loop on AB_ind

        end % loop on i

        relevant_phi_edges_vec = phi_edges_vec(cur_states_indexes);

        % Normalize by the set boundary set sizes
        for AB_ind = 1:length(cur_states_indexes)
            pay_vec(cur_states_start_indexes(AB_ind):cur_states_start_indexes(AB_ind)+2^num_on-1) = ...
                pay_vec(cur_states_start_indexes(AB_ind):cur_states_start_indexes(AB_ind)+2^num_on-1) ./ ...
                phi_edges_vec(cur_states_indexes(AB_ind));
        end

        num_on_is = num_on
    end  % loop on steps

    % Now 'take out' P_mat from the pay_vec:
    I_VEC = floor(log2(cur_states_indexes-1));
    J_VEC = floor(log2(cur_states_indexes - 1-2.^I_VEC));

    I_VEC = I_VEC+1;
    J_VEC = J_VEC+1;


    P_mat_up = pay_vec(cur_states_start_indexes+1);
    P_mat_down = pay_vec(cur_states_start_indexes+2);


    for i=1:length(P_mat_up)
        P_mat(I_VEC(i), J_VEC(i)) = P_mat_up(i)-N/2;  % Subtract to get a zero-sum game
        P_mat(J_VEC(i), I_VEC(i)) = P_mat_down(i)-N/2;
    end



else % Here do simple random walk, which can be easily solved using a linear system of equations


    % First try the 'out of A' approach

    % loop on strategies
    for i=1:N
        for j=i+1:N   
            A = eye(N); b_vec = zeros(N,1); b_vec(i) = 1; 
            
            for k=setdiff([1:N], [i,j]) 
                A(k,find(E(k,:))) = A(k,find(E(k,:)))-1./sum(E(k,:));
            end

%             A
%             b_vec
%             solsol = inv(A)*b_vec
%             size(solsol)
            
            % Solve linear system and average over all vertices
            P_mat(i,j) = sum(inv(A)*b_vec);
            
        end
    end
            
    
    % Make symmetryc and zero-sum
    P_mat = P_mat+P_mat' +N/2 .* eye(N) - N/2;
    
    
    
            
            
   
            
    
    % Second, try the 'into A' aprroach



    sol_mat = zeros(N^2,1);

    for k=1:N % This is the target vertex
        A = eye(N^2);


        for i=setdiff(1:N, k)
            for j=setdiff(1:N, k)
                A((i-1)*N+j, (find(E(i,:))-1)*N+j ) = A((i-1)*N+j, (find(E(i,:))-1)*N+j ) - 1./(sum(E(i,:))+sum(E(j,:)));
                A((i-1)*N+j, (i-1)*N + find(E(j,:))) = A((i-1)*N+j, (i-1)*N + find(E(j,:))) - 1./(sum(E(i,:))+sum(E(j,:)));
            end
        end


%         A
%         sum(A,2)



        b_vec = zeros(N^2,1); b_vec((k-1)*N+1:k*N)=1; b_vec((k-1)*N+k) = 0.5;
%         b_vec'

        sol_mat = sol_mat + inv(A)*b_vec;
    end

    P_mat2 = reshape(sol_mat, N, N) - N/2

    P_mat = P_mat2;

    % See if two approached are equivalent
   


end


