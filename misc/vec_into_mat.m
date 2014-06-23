% Convert a 1-d vector into a matrix (written by Libi Hertzberg) 
function mat = vec_into_mat(vec, m)

n = (length(vec))/m;

mat = (reshape(vec, m, n));
