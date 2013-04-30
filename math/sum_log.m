% compute sum(log(x_i))
function s = sum_log(x_vec, y_vec)

if(~exist('y_vec', 'var')) % sum a vector 
    m = max(x_vec); % maximum log 
    s = log(sum(exp(x_vec-m))) + m;
else % sum two vector element-wise. 
    m_vec = max(x_vec, y_vec);
    s = log(exp(x_vec-m_vec) + exp(y_vec-m_vec)) + m_vec;
end
