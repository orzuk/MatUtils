% We generate all the 'almost clicks' (with all possible orientations)
% such that G is NOT contained in them (i.e. not a sub-model of them).
% Note: this is an exponentially heavy function - do not run with big
% graphs.
function [Cij_dags C_all_edges] = get_all_Cij_for_dag(G)

N = length(G);
Cij_dags  = {};
Cij = zeros(N); % The 'almost clique' that we keep changing
dags_ind=1;
is_sub = 1; % Assume it is
% Go over all vertices pairs
ord_vec = zeros(1,N); orient_bits = zeros(1,N-2);
for i=1:N
    for j=i+1:N
        non_ij_vec = [1:i-1,i+1:j-1,j+1:N];
        ord_vec(1)=i; ord_vec(2) = j;
        % Now go over all 2^(N-2) orientations
        for orientation = 0:2^(N-2)-1
            % Loop over the bits of the orientation and decide if
            % to put a V-struct or not:
            Cij = zeros(N); % The 'almost clique' that we keep changing

            % First determine the ordering of the variables: first
            % i and j, then the non-V vertices and then the V vertices:
            for k=1:N-2
                orient_bits(k) = bitget(orientation, k);
            end
            num_non_vs = sum(orient_bits);
            ord_vec(3:2+num_non_vs) = find(orient_bits); % The non-Vs
            ord_vec(3+num_non_vs:N) = find(orient_bits == 0); % The Vs
            %                    ord_vec_now_is = ord_vec
            % Correct to eliminate i and j
            ord_vec(3:N) = non_ij_vec(ord_vec(3:end));
            %                    ord_vec(find(ord_vec(3:end) >= j)+2)=  ord_vec(find(ord_vec(3:end) >= j)+2)+1;
            %                    ord_vec(find(ord_vec(3:end) >= i)+2) =  ord_vec(find(ord_vec(3:end) >= i)+2)+1;

            %                    ord_vec_is = ord_vec

            % Finished ordering, now draw the edges
            Cij(i,[1:i-1,i+1:j-1,j+1:N]) = 1;  % From i to everywhere except j
            Cij(j,ord_vec(3+num_non_vs:N)) = 1; % From j to all the Vs
            Cij(ord_vec(3:2+num_non_vs), [j,ord_vec(3+num_non_vs:N)]) = 1; % From non Vs to Vs and to j
            Cij(ord_vec(3:2+num_non_vs),ord_vec(3:2+num_non_vs)) = ones(num_non_vs)-tril(ones(num_non_vs)); %From non Vs to themselves (clique)
            Cij(ord_vec(3+num_non_vs:N),ord_vec(3+num_non_vs:N)) = ones(N-num_non_vs-2)-tril(ones(N-num_non_vs-2)); %From non Vs to themselves (clique)

            % Check if its a valid graph
            if(~acyclic(Cij))
                i_is = i
                j_is = j
                orient_is = orientation
                figure; graph_draw(Cij);
                return;
            end
            if(~zuk_is_submodel(G, Cij))
                Cij_dags{dags_ind} = Cij;
                dags_ind = dags_ind+1;
            end

        end
    end
end

% Get all the edges which might be missing and represent them as a graph
C_all_edges = Cij_dags{1}+Cij_dags{1}';

for i=2:length(Cij_dags)
    C_all_edges  = bitand(C_all_edges,  Cij_dags{i}+Cij_dags{i}');
end
for i=1:N
    C_all_edges(i,i)=1;
end
C_all_edges = 1-C_all_edges;


return;


