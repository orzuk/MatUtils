% Add a relative derichlet correction
% Input:
% V - input probability/counts vector (with possible zeros)
% alpha - the relative correction factor
% M - the number of elements in V
%
% Output:
% W - the corrected probability vector
%
function W = derich_correct(V, alpha, M, varargin)

if(alpha == 0) % small check to see that we don't correct by zero
    W=V;
    return;
end
if(iscell(V)) % work on cell array 
    n = length(V); 
    W = cell(n,1);
    for i=1:n
        W{i} = derich_correct(V{i},alpha);
    end
else
    if(nargin == 2) % here the function finds M
        M = size(V, 1); % the dimension (four for pwms) 
        if(M == 1) % special case when V is a vector
            M = size(V, 1);
        end
    end
    W = (V + alpha) ./ (1 + M*alpha);
end


