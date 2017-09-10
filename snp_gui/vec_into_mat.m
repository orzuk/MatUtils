%function mat = vec_into_mat(vec, m)
function mat = vec_into_mat(vec, m)

n = (length(vec))/m;

mat = (reshape(vec, m, n));
