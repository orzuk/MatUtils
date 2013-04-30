% My tests for the BNT , just prelimenary runs (Or)

% Make Some Sinthetic Data
N = 4; % Number Of Points
dag = zeros(N);  % Adjancy Matrix
C = 1; S = 2; R = 3; W = 4;  % Give Names : Cloudy 1 Sprinkler 2 Rain 3 WetGrass 4
dag(C,[R S]) = 1;
dag(R, W) = 1;
dag(S, W) = 1;

discrete_nodes = 1:N;  % All nodes are discrete (this is default ... )
node_sizes = 2*ones(1,N);   % All nodes are binary


% Construct The Bayesian Network
% bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);
bnet = mk_bnet(dag, node_sizes);  % discrete is default


% Parameters : Determine the CPD (Conditional Probability Distribution  )
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);


% rand('state', seed);   % This gives a random CPD
% We can also use sampling according to dirichlet using sample_dirichlet






% Inference 
%%%%%%%%%%%

% Choose inference engine
engine = jtree_inf_engine(bnet);  % Junction Tree


% Create evidence (obsreved r.v.)
evidence = cell(1,N);
evidence{W} = 2;  % This means W = 2

% Add evidence to the engine
[engine, loglik] = enter_evidence(engine, evidence);  % This depends on the algorithm


% Now, compute p = Pr(S=2/W=2) : 
marg = marginal_nodes(engine, S);
marg.T % Display The Marginal Distribution
p = marg.T(2); % Should be p = 0.4298







%%%%%%%%%%%%%%
% Printing Utilities 
%
% graph_to_dot(bnet.dag, 'filename', 'foo.dot'); dos('dot -Tps foo.dot -o foo.ps'); dos('gsview32 foo.ps')
%
%%%%%%%%%%%%%%











