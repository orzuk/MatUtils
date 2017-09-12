% Check if a vector is a column vector
function a = is_column(v)
a = (size(v,2) == 1);