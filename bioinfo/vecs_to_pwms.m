% Convert 4L vectors to 4*L matrices. 
% The format of the vector is: 
% (pos1: A C G T  , pos2: A C G T, ....)
% 
% The input: 
% V - a matrix of size 4L*m (unless working in cell-mode, then output&input are cell-arrays)
%
% The output: 
% P - a 3-dim array of pwms, of size 4*L*M
%
function P = vecs_to_pwms(V)

if(iscell(V)) % here we handle also a cell-array of pwms
    P = cell(1,length(V));
    for i=1:length(V)
        P{i} = vecs_to_pwms(V{i});
    end
    return;
end
L = size(V,1) / 4;
P = reshape(V, 4, L, size(V,2));


