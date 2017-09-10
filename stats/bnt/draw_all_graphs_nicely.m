% draw_all_graphs_nicely


N=4; num_graphs = 11;


dims_vec = [4 5 6 7 7 6 7 8 8 9 10]; % The last one 10 is a ganble


p=0.3; XY = [p 1-p p 1-p; 1-p 1-p p p];

figure; subplot(4,3,1); hold on;

subplot(4,3,2); hold on; title('Dimensions of all 4-nodes non-isomorphic undirected graphs on binary r.v.s');

Emats = {};

for i=1:num_graphs
    Emats{i} = zeros(4);
end

Emats{2}(1,2)=1;

Emats{3}(1,2)=1; Emats{3}(2,3)=1;

Emats{4}(1,2)=1; Emats{4}(2,3)=1; Emats{4}(2,4)=1;

Emats{5}(1,2)=1; Emats{5}(2,3)=1; Emats{5}(1,3)=1;

Emats{6}(1,4)=1; Emats{6}(2,3)=1; 

Emats{7}(1,2)=1; Emats{7}(1,4)=1; Emats{7}(2,3)=1;

Emats{8}(1,2)=1; Emats{8}(1,3)=1; Emats{8}(1,4)=1; Emats{8}(2,3)=1;

Emats{9}(1,3)=1; Emats{9}(1,4)=1; Emats{9}(2,3)=1; Emats{9}(2,4)=1;

Emats{10}(1,2)=1; Emats{10}(1,3)=1; Emats{10}(1,4)=1; Emats{10}(2,3)=1; Emats{10}(2,4)=1;

Emats{11}(1,2)=1; Emats{11}(1,3)=1; Emats{11}(1,4)=1; Emats{11}(2,3)=1; Emats{11}(2,4)=1; Emats{11}(3,4)=1;

for i=1:num_graphs
    Emats{i} = Emats{i} + Emats{i}';
end

for i=1:num_graphs
    subplot(4,3,i); hold on;
    % Now show the graphs
    graph_draw_given_layout(Emats{i}, XY);
    text(0.1, 0.9, ['G #' num2str(i)]);
    text(0.1, 0.1, ['dim :' num2str(dims_vec(i))]);
end

