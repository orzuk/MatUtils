% Generate a graph representing the possible 
% transitions between binary sequences with n-1
% bits overlap
function [E]=generate_seqs_graph(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the transitions sequences with one additional 'overall' state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E = zeros(2^N);

for i=[0:2^N-1]
   E(i+1, bitshift(i,1,N)+1)=1; E(i+1, bitshift(i,1,N)+1+1)=1; % To whom can we pass    
end

size(E)

%E(end,:)=1;

%E=min(E+E', 1);