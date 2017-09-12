% My plot graph function for the BNT , uses the graphViz freeware package (Or)
function ret = my_plot_graph(dag)



% Transfer the DAG matrix into a 'dot' format file
graph_to_dot(dag, 'filename', 'foo.dot');



% Convert the 'dot' file into a postscript file using the GraphViz package
dos('dot -Tps foo.dot -o foo.ps & ');

%  Open the postscript file containing the graph using Ghostview,
dos('D:\Ghostgum\gsview\gsview32 foo.ps &');


ret = 0; % return 0 





%%%%%%%%%%%%%%
% Printing Utilities 
%
% graph_to_dot(bnet.dag, 'filename', 'foo.dot'); dos('dot -Tps foo.dot -o foo.ps'); dos('gsview32 foo.ps')
%
%%%%%%%%%%%%%%











