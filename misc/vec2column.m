% Force a vector to be a column vector
function v_col = vec2column(v)
if(size(v,2) > max(1, size(v,1)))
    v_col = v';
else
    v_col = v;
end
    

