% Compute information content for pwms vectors
%
% Input: pwms - a set of position weight matrices 
% 
% Output: ic - information content (in bits) of each pwm
% 
function ic = pwms_information_content(pwms)

V = pwms_to_vecs(pwms);
if(iscell(V))
    n = length(V);
    ic = zeros(n,1);
    for i=1:n
        ic(i) = 0.5*length(V{i}) - entropy(V{i});
    end
else
    ic = 2*size(pwms, 2) - entropy(V);
end

