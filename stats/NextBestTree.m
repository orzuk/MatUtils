% Computes the 'next best' spanning tree
% after the optimal one. This is useful when we want to
% generate a distribution from a tree, and see how close it is
% to some 'wrong' tree
% This is a more efficient version. It follows Gabow's algorithm:
% Gabow, H.N., Two algorithms for generating weighted spanning trees in order.
% SIAM Journal on Computing, Vol. 6, No. 1. (1977), pp. 139-150
%
% Input:
% Mat - matrix of pairwise distances
% Tree - current spanning tree
%
% Output:
% NextTree - the next best spanning tree
% ScoreDiff - difference in score between input and output trees
%
function [NextTree ScoreDiff] = NextBestTree(Mat, Tree, algorithm_used, varargin)


n = size(Tree,1); % # of nodes
Tree = triu(Tree);
CompTree = triu(1-Tree-eye(n)); % complement tree - used Gabow's algorithm!

TreeEdges = find(Tree);
CompTreeEdges = find(CompTree);

[SortedTreeEdgesVals SortedTreeEdgesInds] = sort(Mat(TreeEdges)); % Sort edges by weight
[SortedCompTreeEdgesVals SortedCompTreeEdgesInds] = sort(Mat(CompTreeEdges),'descend');

DiffMat = repmat(SortedTreeEdgesVals, 1, length(CompTreeEdges)) - ...
    repmat(SortedCompTreeEdgesVals, 1, length(TreeEdges))';  % Get the payments for each pair of edges
CheckedVec = sum(DiffMat < 0); % We disregard the negatives


if(~exist('algorithm_used', 'var') || isempty(algorithm_used))
    algorithm_used = 0; % 0 - new way, 1 Gabow's algorithm (dosen't work), 2 slow implementation
end
switch algorithm_used
    case 0  % sort differences and try one by one
        [DiffMatSorted DiffMatPerm] = sort(-DiffMat(:), 'descend'); % , 'descend'); % just sort differences
        for checked_ind = 1:(n-1)^2*(n-2)/2
            I = mod_max(DiffMatPerm(checked_ind), (n-1)); % edge to remove
            J = div(DiffMatPerm(checked_ind)-1, (n-1))+1; % edge to add
            
            TempTree = Tree; TempTree(TreeEdges(SortedTreeEdgesInds(I))) = 0;
            TempTree(CompTreeEdges(SortedCompTreeEdgesInds(J))) = 1;
            if(acyclic(TempTree+TempTree', 0)) % test if we have a tree
                NextTree = TempTree; NextTree = NextTree + NextTree'; % Make undirected
                ScoreDiff = -DiffMatSorted(checked_ind);
                return;
            end
        end
    case 1  % Gabow's algorithm. Still a bug in implementation
        for checked_ind = 1:(n-1)^2*(n-2)/2 % indices: n-1 edges in tree x (n-1) over 2 in complement
            %    ind_is = checked_ind
            cur_js = [1 find(CheckedVec(1:end-1) > CheckedVec(2:end))+1];
            %     if(CheckedVec(cur_js(end)) == 0) % don't try to replace last one
            %         cur_js = cur_js(1:end-1);
            %     end
            CurDiffs = SortedTreeEdgesVals(CheckedVec(cur_js)+1) - SortedCompTreeEdgesVals(cur_js);
            [MinDiff next_j] = min(CurDiffs); next_j = cur_js(next_j); % find best difference
            
            TempTree = Tree; TempTree(TreeEdges(SortedTreeEdgesInds(CheckedVec(next_j)+1))) = 0;
            TempTree(CompTreeEdges(SortedCompTreeEdgesInds(next_j))) = 1;
            if(acyclic(TempTree+TempTree', 0)) % test if we have a tree
                NextTree = TempTree; NextTree = NextTree + NextTree'; % Make undirected
                ScoreDiff = MinDiff;
                return;
            end
            CheckedVec(next_j) = CheckedVec(next_j)+1;   % Update chekced counts (why addition?)
        end
        
    case 2 % old slow implementation
        DiffMat = repmat(Mat(TreeEdges), 1, length(CompTreeEdges)) - ...
            repmat(Mat(CompTreeEdges), 1, length(TreeEdges))';  % Get the payments for each pair of edges
        
        CompTreeEdges_I = mod(CompTreeEdges-1,n)+1; CompTreeEdges_J = ceil(CompTreeEdges./n); % get edges of tree complement
        
        % Try other way (slower but correct)
        CyclicMat = zeros(size(DiffMat,1), size(DiffMat,2));
        for j=1:length(TreeEdges) % loop on edges
            TempTree = Tree; TempTree(TreeEdges(j)) = 0; TempTree = TempTree+TempTree';
            for i=1:length(CompTreeEdges) % loop over pairs of possible problematic edges
                TempTree2 = TempTree; TempTree2(CompTreeEdges(i)) = 1;
                TempTree2((CompTreeEdges_I(i)-1)*n + CompTreeEdges_J(i)) = 1; % Make transpose also 1
                CyclicMat(j,i) = 1-acyclic(TempTree2, 0);
            end
        end
        DiffMat = DiffMat + 99999999999 .* CyclicMat; % Take only valid edges removals
        [ScoreDiff I] = min(DiffMat); [ScoreDiff J] = min(ScoreDiff); I = I(J); % Get minimal values of the matrix
        I = TreeEdges(I); J = CompTreeEdges(J);
        
        NextTree = Tree;
        NextTree(I) = 0; NextTree(J) = 1; % replace the edges
        NextTree = NextTree + NextTree'; % Make undirected
end


