%function [C, AI, BI] = intersect_order_by_first_gr(A, B)
function [C, AI, BI] = intersect_order_by_first_gr(A, B)

[c, ai, bi] = intersect(A, B);

if(length(c) == length(A))
    C = A;
    AI = [1:length(A)];
    BI = AI;
    BI(ai) = bi;
else
    t = 'intersect_order_by_first_gr: not same size'
    [AI, ind_sorted] = sort(ai);
    C = A(AI);
    BI = bi(ind_sorted);
end