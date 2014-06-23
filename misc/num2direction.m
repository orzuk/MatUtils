% Convert a numeric directions vector to marked by 'LEFT' and 'RIGHT' 
function d = num2direction(d_num)


AssignGeneralConstants;

n = length(d_num);
d = mat2cell(repmat('LEFT', n,1), ones(n,1));
right_inds = find(d_num == RIGHT);
for i=vec2row(right_inds)
    d{i} = 'RIGHT';
end

