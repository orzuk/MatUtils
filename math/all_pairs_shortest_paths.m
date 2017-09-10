% Compute shortest paths for all pairs of points in a graph
%
% Input: 
% E - graph adjacency matrix
% 
% Output: 
% apsp - the shortest paths between each pair of vertices
% 
function apsp = all_pairs_shortest_paths(E)

N = length(E);

apsdNew = -ones(N) + eye(N) + 2.*E; apsd = apsdNew;
apspNew = -ones(N);
for i=1:N
    for j=1:N
        if(E(i,j))
            apspNew(i,j)=1;
        end
    end
end
apsp = apspNew; 
counted = N + sum(sum(E)); % Start with num. edges + N
while (counted < N^2)
    for i=1:N
        for j=setdiff([1:N],i)
            for k=setdiff([1:N],[i,j])
                if((apsd(i,j)==(-1)) && E(i,k) && (apsd(j,j) >=  0))
                    apsdNew(i,j) = apsd(k,j)+1;
                    apspNew(i,j) = k;
                    counted = counted + 1;
                end
            end
        end
    end
    apsd = apsdNew; apsp = apspNew;  
end



% Dan's python code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% def allPairsShortestPaths(E):
%     N = len(E)
% 
%     apsdNew = -ones(E.shape)+identity(N)+2*E
%     apsd = array(apsdNew)
% 
%     apspNew = -ones(E.shape)
%     for i in xrange(N):
%         for j in xrange(N):
%             if E[i,j]:
%                 apspNew[i,j] = i
%     apsp = array(apspNew)
%     counted = N+sum(sum(E))
%     while counted < N**2:
%         for i in xrange(N):
%             for j in xrange(N):
%                 if i==j: continue
%                 for k in xrange(N):
%                     if i==k: continue
%                     if j==k: continue
%                     if apsd[i,j]==-1 and E[k,j] and apsd[i,k]>=0:
%                         apsdNew[i,j] = apsd[i,k]+1
%                         apspNew[i,j] = k
%                         counted += 1
%         apsd = array(apsdNew)
%         apsp = array(apspNew)
%     return apsp
% 
% def printPath(apsp,i,j):



