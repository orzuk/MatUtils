% addpath 'C:\Weizmann\Research\ipod';
% addpath 'C:\weizmann\Research\BayesianNetworks\BNsoftware\KevinMurphy\FullBNT_2005_01_31\FullBNT\graph';
% addpath 'C:\weizmann\Research\BayesianNetworks\BNsoftware\KevinMurphy\FullBNT_2005_01_31\FullBNT\KPMtools';
addpath  'E:\Research\networks\BNsoftware\KevinMurphy\FullBNT\graph';
addpath  'E:\Research\ipod';

N=5; p=0.5;

G1 = create_random_DAG(N,p/2); G2 = create_random_DAG(N,p);

graph_draw(G1); figure; graph_draw(G2);

is_sub = my_is_submodel(G1, G2)


% Now also check for graph G all of their Cij's
 [Cij_dags C_all_edges] = get_all_Cij_for_dag(G1);
figure; hold on; graph_draw(G1); title('G graph');
figure; hold on; graph_draw(C_all_edges); title('All Edges');

[x, y, labels] = draw_dot_two_graphs(G1, C_all_edges); %, labels);

