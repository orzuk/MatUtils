% Compute entropy for histogram 
%
% Input:
% x1 - x values for first histogram
% p1 - probablilites for first histogram
% x2 - x values for second histogram
% p2 - probablilites for second histogram
% scale_factor - resolution increase in output histogram
% KL_mode - whether to compute relative entropy or just log-likelihood part
% pseudo_count - new! allow pseudo counts for empty cells, to avoid log(0)
%
% Output:
% KL - relative entropy: D(p1 || p2) 
% 
function KL = relative_entropy_hist(x1, p1, x2, p2, KL_mode, pseudo_count)

if(~exist('KL_mode', 'var') || isempty(KL_mode)) % set default normalization
    KL_mode = 1;
end
if(~exist('pseudo_count', 'var') || isempty(pseudo_count)) % Default: no pseudo counts 
    pseudo_count = [0 0];    
end
if(isscalar(pseudo_count))
    pseudo_count = repmat(pseudo_count, 2, 1);
end

x1 = x1(p1>0); p1=p1(p1>0); % leave only non-zero elements 
x2 = x2(p2>0); p2=p2(p2>0); % leave only non-zero elements 


if(max(pseudo_count)>0) % add to x: 
    x = union(x1, x2);
    [~, I1, J1] = intersect(x1, x); p1_new = zeros(size(x)); p1_new(J1) = p1(I1); p1=p1_new + pseudo_count(1); x1=x;
    [~, I2, J2] = intersect(x2, x); p2_new = zeros(size(x)); p2_new(J2) = p2(I2); p2=p2_new + pseudo_count(2); x2=x;
end

p1=p1./sum(p1); p2=p2./sum(p2); % normalize 
[x_inter, I1, I2] = intersect(x1, x2); 
if(length(x_inter) < length(x1)) % there is an element not appearing in x2 !! 
    KL = inf;
    error('Error !!! reference distribution has zeros where first distribution is non-zero!!'); 
end
if(KL_mode) % compute KL
    KL = sum( vec2row(p1(I1)) .* vec2row(log(p1(I1)) - log(p2(I2))) ); 
else % compute just P*log(Q) part (also switch sign!!!)
    KL = sum( vec2row(p1(I1)) .* vec2row(log(p2(I2))) );     
end
