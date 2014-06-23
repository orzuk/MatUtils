% Transfer a number in a decimal format to binary format
% 
function [b]=dec_to_bin(a)


n = max(ceil(log2(a+1))); b=zeros(1,length(a));

for i=0:n-1
    b = b + (mod(floor(a./2.^i), 2)).*10.^i;
end

