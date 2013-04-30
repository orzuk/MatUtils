% Plot all graphs of size 4/5 along with their dimension
N = 5;

ug = enumerate_unlabeled_graphs(N); %

num_graphs = size(ug, 1)

dim_vec = zeros(1, num_graphs);

% Set the dims manually
dim_vec = [1 2];


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

for i=1:num_graphs
    subplot(6,6,i); hold on;
    graph_fixed_draw( EMats(:,:,i));
end


