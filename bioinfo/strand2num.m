% Convert a strand vector marked by '+' and '-' to numeric
function s_num = strand2num(s)

Assign24MammalsGlobalConstants;

n = length(s); 
s_num = zeros(n,1) + POS_STRAND; 
s_num(strmatch('-', s)) = REV_STRAND;






