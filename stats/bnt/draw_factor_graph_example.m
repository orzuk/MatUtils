% Draw   an example of a factor graph

N=25+3;

E=zeros(N); 

E(26,26)=1; E(27,27)=1; E(28,28)=1; % Color the factors

E(5:5:25,26)=1; % add row constrain
E(21:25,27)=1;  % add column constrain
E(11,28)=1; E(19,28)=1; % add edge constrain


E=E+E'; % Make symmetric

% Prepare x and y vecs
p=0.12;
y_vec = repmat(p:p:5*p, 1, 5); y_vec(26:N)= [15*p/2, 5*p, 4*p/2];
x_vec = reshape(repmat(p:p:5*p,5,1),1,25); x_vec(26:N)= 1-p;


labels = int2str( [15:-1:11, 25:-1:21, 35:-1:31,45:-1:41,55:-1:51]'); %  str2num('row'), str2num('col'); str2num('edge')]

labels = cellstr(labels)   %  labels = cellstr(char(zeros(N,1)+double('+')));

labels{26} = 'row'; labels{27} = 'column'; labels{28} = 'edge';

% Now plot with the x,y coordinates
graph_fixed_draw(E, x_vec', y_vec', labels);