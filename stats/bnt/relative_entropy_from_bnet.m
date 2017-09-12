% Compute the relative entropy between dag G and distribution P.
% This is difficult. How do we minimize?
% Currently the function computes the other (easier) direction
function KL = relative_entropy_from_bnet(G, P)


Q = project_on_bnet(P, G); % first project

KL = sum(P .* log(P)) - sum(P .* log(Q));  % Now compute KL