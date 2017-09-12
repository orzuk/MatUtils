% Draws a phylogenetic tree with weights on the branches, 
% where the weights are not the branch length but give some extra
% information like rate of evolution for a specific set etc.
%
% Input: 
% phytree_file - input file with phylogenetic tree
% weighted_flag - whether to plot a weighted tree (by colors)
% w - the weights themselves (optional)
%
function dummy = draw_phytree(phytree_file, weighted_flag, w, varargin)

% read file (assume in Newick format)
if(isa(phytree_file, 'char'))
    Tree = phytreeread(phytree_file);
else % enable giving a tree as input 
    Tree = phytree_file;
end

if(~exist('weighted_flag', 'var'))
    weighted_flag = 1; % flag saying if to draw the tree in colors
end
if(weighted_flag)
    if(~exist('w', 'var') || isempty(w)) % set weights inside function 
        w = weights(Tree); w = [w' w(1:end-1)']; n = length(w);
        w = randn(n,1).*0.2+1; w(1) = 19;  
    end
    %w(2) = 1.61234; w(3) = 1.42423; w(4) = 1.644;
    %w(5) = 1.6331; w(22) = 1.323; w(23) = 1.5123; w(24) = 1.443;
%    w(:) = 1;
    H = phytree_plot(Tree, 'weights', w, 'orientation', 'left', 'type', 'square'); % use my colored edges plot
    title('HAR1 - 118bp');
else
    H = plot(Tree); % just use matlab bioinfo-toolbox plot (un-colored)

end

dummy = 0;
