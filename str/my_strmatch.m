% Strmatch for cell, deals with empty cells (from Tal Shay)
function A = my_strmatch(str, strarray, flag)

I = find(cellfun(@ischar, strarray));  % New: Changed! (don't check if empty. Check if char)
if (nargin == 2)
    B = strmatch(str, strarray(I));
else
    B = strmatch(str, strarray(I), flag);
end
A = I(B);
