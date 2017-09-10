% Compute the Normal Fan of a 2-dimensional Polytope
function NF = NormalFan(NP)

edges = diff(NP);

NF =  [edges(:,2), - edges(:,1)];
NF = NF./ repmat(sqrt( edges(:,2).^2 + edges(:,1).^2), 1, 2); % compute edges of the normal fan
% eq_inds = find(abs(NP(:,1)) == abs(NP(:,2)));
% if(~isempty(eq_inds))
%     NF(eq_inds,1) = sign(NP(eq_inds,2)) ./sqrt(2);
%     NF(eq_inds,2) = sign(NP(eq_inds,1)) ./ sqrt(2);
% end


