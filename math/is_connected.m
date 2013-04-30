% Check if an undirected graph is conneted
function IC = is_connected(G)

IC=isempty(find(bfs(G,1) < 0));

