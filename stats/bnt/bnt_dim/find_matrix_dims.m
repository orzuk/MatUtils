% Plot all graphs of size 4/5 along with their dimension
N = 6;

ug = enumerate_unlabeled_graphs(N); %

num_graphs = size(ug, 1)

dim_vec = zeros(1, num_graphs);

% Set the dims manually
dim_vec = zeros(1,num_graphs);
num_edges_vec = zeros(1,num_graphs);
num_unconnected_vec = zeros(1,num_graphs); 

EMats = zeros( N,N, num_graphs);
i=1;
for j=1:N
    for k=j+1:N
        %   for s=1:
        EMats(j,k,:) = ug(:,i);
        EMats(k,j,:) = ug(:,i);  % Symmetric
        i=i+1;
    end
end





figure; hold on;

if(N==4)
    grid_x=3; grid_y=4;
end
if(N==5) 
    grid_x=6; grid_y=6;
end
if(N==6)
    grid_x=12; grid_y=13;
end

for i=1:num_graphs

    % Build the matrix of 0's and 1's indicating
    % which potential appears in which probability
    cur_num_edges = sum(sum(EMats(:,:,i)))/2; num_edges_vec(i) = cur_num_edges; 
    num_unconnected_vec(i) = sum(sum(EMats(:,:,i)) == 0);
    A = zeros(4*cur_num_edges,2^N);

    % Loop over potentials
    ind=1;
    for  j=1:N
        j_bit = mod(floor([0:2^N-1] ./ 2^(j-1)), 2);
        for k=j+1:N
            k_bit = mod(floor([0:2^N-1] ./ 2^(k-1)), 2);
            jk_4 = 2*j_bit + k_bit;
            if(EMats(j,k,i)) % There is a potential

                A(ind,find(jk_4 == 0)) = 1;
                A(ind+1,find(jk_4 == 1)) = 1;
                A(ind+2,find(jk_4 == 2)) = 1;
                A(ind+3,find(jk_4 == 3)) = 1;

                ind=ind+4; % Each time step 4 times !!!
            end
        end
    end

    dim_vec(i) = rank(A);

    subplot(grid_x,grid_y,i); hold on;
    graph_fixed_draw( EMats(:,:,i)); title(['DIM is ' num2str(dim_vec(i)+num_unconnected_vec(i)-1)]);
end


figure; plot(num_edges_vec, dim_vec+num_unconnected_vec, '+');
