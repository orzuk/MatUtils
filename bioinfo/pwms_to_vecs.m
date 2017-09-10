% Convert 4*L matrices to vectors. 
% The format of the vector is: 
% pos1: A C G T  , pos2: A C G T, ....
% 
% Input: 
% P - a 3-dim array of pwms, of size 4*L*M (unless working in cell-mode, then output&input are cell-arrays)
%
% The output: 
% V - a matrix of size 4L*m
function V = pwms_to_vecs(P)

if(iscell(P)) % here we handle also a cell-array of pwms
    V = cell(1,length(P));
    for i=1:length(P)
        V{i} = pwms_to_vecs(P{i});
    end
    return;
end
L = size(P,2); 
V = reshape(P, L*4, size(P,3));


