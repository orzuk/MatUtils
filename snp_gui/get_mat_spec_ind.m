%function vals = get_mat_spec_ind(mat, I,J)
function vals = get_mat_spec_ind(mat, I,J)

num_rows = size(mat, 1);
num_columns = size(mat, 2);
vec = reshape(mat,1, num_rows*num_columns);

ind = num_rows*(J-1) + I;
vals = vec(ind);