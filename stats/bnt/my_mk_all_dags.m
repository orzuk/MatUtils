function Gs = my_mk_all_dags(N, only_connected)
% my_mk_all_dags generate all DAGs on N variables
% G = my_mk_all_dags(N)
%
% G = mk_all_dags(N, order) only generates DAGs in which node i has parents from 
% nodes in order(1:i-1). Default: order=[] (no constraints).
%
% G{i} is the i'th dag
%
% Note: the number of DAGs is super-exponential in N, so don't call this with N > 4.

if nargin < 2, order = []; end

use_file = 0;

global BNT_HOME
fname = sprintf('%s/DAGS%d.mat', BNT_HOME, N);
if use_file & exist(fname, 'file')
  S = load(fname, '-mat');
  fprintf('loading %s\n', fname);
  Gs = S.Gs;
  return;
end



% the dag
dag = zeros(N);

indy = (1:N*N);
up = find(mod(indy-1, N)+1 < floor( (indy-1)/N)+1);


% possible digraphs
m = 2^(N*(N-1)/2);

% all perms of size n
p = perms(1:N);


directed = 1; undirected = 0;

Gs = {};
%jjj = 1;

GG = zeros(m, N*N);

TOT = zeros(1, N*N);

ind =  ind2subv(2*ones(1,N*(N-1)/2), 1:m);
%enumerate all permutations on N 
for i = 1:length(p)
    % enumerate all of the graphs
    for j=1:m
        % Assign the graph
        dag(up) = ind(j,:)-1;
        
        % Check if acyclic
        if acyclic(dag, directed)
            % Check if connected
            % we need to permute ...
                GG(j,:) = reshape( dag(p(i,:), p(i,:)),1, N*N ) ; %jjj = jjj+1;
        end
    end
    
    % Here we must do unique!!!
    GG = unique(GG, 'rows');
    %sprintf('do perm %d out of %d\n', i, length(p)) 
    TOT = union(TOT, GG, 'rows');
    
end



% Remove all non connected graphs
if(only_connected)
    % Remove all non connected graphs
    j = 1;
    for i=1:length(TOT)
        if(connected_graph( reshape(TOT(i,:), N, N), undirected ))
            Gs{j} = reshape(TOT(i,:), N, N);
            j = j+1;
        end
    end
else
    j = 1;
    for i=1:length(TOT)

        Gs{j} = reshape(TOT(i,:), N, N);
        j = j+1;
    end
end

        
           
    

