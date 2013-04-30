% Get source nodes of a directed graph 
function s = get_graph_sinks(A)

s = find(sum(A - diag(diag(A)),2)==0);

