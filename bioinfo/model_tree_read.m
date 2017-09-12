% Reads a phylogenetic model: this include both tree and rate matrices
% For some reason the phytreeread part doesn't work for some files in pc
% (only in unix) and throws matlab 
%
function [Q PI tree] = model_tree_read(model_file_name)

Assign24MammalsGlobalConstants(); 

R = textread(model_file_name, '%s', 'delimiter', '\n'); 

PI = str2nums(R{1}); 
Q = zeros(4); 
for i=1:4
    Q(i,:) = str2nums(R{i+2});
end

R{end} = strdiff(R{end}, 'TREE: ');

if(machine == UNIX)
    tree = phytreeread(R{end});
else % a hack to deal with PC's problem with reading strings with phytreeread
    savecellfile(R(end), 'tmp_tree.txt');
    tree = phytreeread('tmp_tree.txt');
    delete('tmp_tree.txt');
end


