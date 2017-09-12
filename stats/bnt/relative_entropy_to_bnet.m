function KL = relative_entropy_to_bnet(P, G)
% Compute the relative entropy between distribution P and dag G
Q = project_on_bnet(P, G); % first project
KL = sum(P .* log2(P)) - sum(P .* log2(Q));  % Now compute KL
