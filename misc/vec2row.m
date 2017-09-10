% Force a vector to be a row vector
function v_row = vec2row(v)
if(size(v,1) > max(1,size(v,2)))
    v_row = v';
else
    v_row = v;
end
