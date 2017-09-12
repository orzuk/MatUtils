%function tea_cycle_riddle()

n = 20; % number of people
E = generate_random_graph(0, n, 0.2); % generate adj matrix
G = E; 
V = rand(n,1) > 0.5; % randomize initial conditions 

Z = ones(n,1);
iters = 10; full_figure; 
for i=1:iters
       G = set_diag(G, V); 
   %figure; 
   graph_draw(G); pause(3);
   W = sign(E*(2*V-Z));
   V(W==1) = 1; V(W==-1) = 0; % take the majority 

end
