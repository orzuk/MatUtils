% Compute the hamming weight of a vector or a matrix
function H = hamming_weight(x)

% H = zeros(size(x,1), size(x,2));
% L = max(ceil(log2(x(:)+1)));
% for i=0:L
%     H = H + mod(  floor(x/2^i), 2); 
% end

if(isa(x, 'vpi'))
    H = sum(vpi2bin(x) == '1')
else
    H = sum(dec2bin(x) == '1')
end

%H=sum(mod(floor(x./2.^[0:31]), 2));