% Get the reverse-complement of a pwm
function P_rc = pwmrcomplement(P)
if(iscell(P))
    n = length(P);
    P_rc = cell(n,1);
    for i=1:n
        P_rc{i} = pwmrcomplement(P{i});
    end
else
    P_rc = P(end:-1:1,end:-1:1);
end
