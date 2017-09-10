% Draw two graphs one on top of the other.
% [x, y, labels] = draw_dot(adj1, adj2, lables)  
% The first graph determines the layout !
% by Or Zuk May 2006.
%
function [x, y, labels] = draw_dot_two_graphs(adj1, adj2, labels);

adj1_only = max(adj1 - adj2,0);
adj2_only = max(adj2 - adj1,0);
and_adj = adj1 .* adj2; % This graph contains and edge only if both graphs contain it!
or_adj = min(adj1 + adj2,1); adj_sym = min(or_adj+or_adj',1); % make directed edges undirected. Use the union of graphs!

[n,m] = size(adj_sym);
if n ~= m, warning('not a square adjacency matrix!'); end
if isequal(triu(adj_sym,1),tril(adj_sym,-1)'), directed = 0; else, directed = 1; end
adj_sym = double(adj_sym > 0);    % make sure it is a binary matrix cast to double type
% to be platform independant no use of directories in temporary filenames
tmpDOTfile = '_GtDout.dot';           tmpLAYOUT  = '_LAYout.dot';
graph_to_dot(adj_sym, 'directed', directed, 'filename', tmpDOTfile); % save in file
if ispc, shell = 'dos'; else, shell = 'unix'; end                %  Which OS ?
%cmnd = strcat(shell,'(''neato -V'')');    % request version to check NEATO is there
%status = eval(cmnd);
%if status == 1,  warning('DOT/NEATO not accessible'); end
%  put all your favorite  NEATO attributes  here
neato = '(''neato -Tdot  -Gmaxiter=25000 -Gregular'; % -Gstart="regular" -Gregular
neato = strcat([neato '-Gminlen=5 -Goverlap=false ']);   % minimal edge length, no overlap
if n > 100   % some extra NEATO options for over populated graphs
    neato = strcat([neato '-x']);
end



cmnd = strcat([shell neato ' -o' tmpLAYOUT ' ' tmpDOTfile ''')']);    % -x compact
status = eval(cmnd);                 %  get NEATO to layout

[trash, names, x, y] = dot_to_graph(tmpLAYOUT);  % load NEATO layout
num_names = str2num(char(names))';
nam_len = length(names);
if nam_len < n  % plot singletons without coordinates all together in a lower left
    num_names(nam_len+1:n) = my_setdiff(1:n, num_names);
    x(nam_len+1:n) = 0.05*ones(1,n-nam_len);
    y(nam_len+1:n) = 0.05*ones(1,n-nam_len);
end
[ignore,lbl_ndx] = sort(num_names);  % recover from dot_to_graph node_ID permutation
x = x(lbl_ndx); y = y(lbl_ndx);
if nargin == 2                                   % if no labels were provided
    labels = names(lbl_ndx);
end
% now pick a healthy font size and plot
if n > 40, fontsz = 7; elseif n < 12, fontsz = 12; else fontsz = 9; end
figure; clf; axis square      %  now plot

% Draw the first graph
[x, y, h] = graph_draw(adj1_only, 'node_labels', labels, 'fontsize', fontsz, ...  % Draw the original (possibly directed) graph
    'node_shapes', zeros(size(x,2),1), 'X', x, 'Y', y, 'linecolor', 'g');
% Draw on top of it the second graph
[x, y, h] = graph_draw(adj2_only, 'node_labels', labels, 'fontsize', fontsz, ...  % Draw the original (possibly directed) graph
    'node_shapes', zeros(size(x,2),1), 'X', x, 'Y', y, 'linecolor', 'b');
% Draw on top of it the union graph
[x, y, h] = graph_draw(and_adj, 'node_labels', labels, 'fontsize', fontsz, ...  % Draw the original (possibly directed) graph
    'node_shapes', zeros(size(x,2),1), 'X', x, 'Y', y, 'linecolor', 'r', 'linewidth', 3);
%legend('G1 only', 'G2 only', 'G1 and G2');
xlabel('green - G1, blue - G2, thick red - Both graphs');
delete(tmpLAYOUT); delete(tmpDOTfile);     % clean up temporary files
