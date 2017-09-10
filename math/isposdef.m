% Check if a square matrix is positive definite
function f = isposdef(A)
f = min(eig(A)) > 0;

