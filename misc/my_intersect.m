% Intersect cells, deal with empty cells (from Tal Shay) 
function [C IA IB]  = my_intersect(A,B)

I_fullA = find(~cellfun('isempty', A));
fullA = A(I_fullA);

I_fullB = find(~cellfun('isempty', B));
fullB = B(I_fullB);
[C C_fullA C_fullB] = intersect(fullA, fullB);

IA = I_fullA(C_fullA);
IB = I_fullB(C_fullB);
