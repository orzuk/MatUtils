% Get source nodes of a directed graph 
function s = get_graph_sources(A)

s = find(sum(A - diag(diag(A)))==0);

