% Transfer a number in a decimal format to a different base binary format
% 
% Input: 
% a - number in decimal base
% base - what base do we want to represent a in
% num_digits - ensure a certain number of digits (allow for leading zeros)
%
% Output: 
% b - vector of digits representing the number in new base
%
function b=my_dec2base(a, base, num_digits)

if(~exist('num_digits', 'var'))
    num_digits = max(ceil(log2(a+1)./log2(base)));
end
b=zeros(num_digits,length(a));

for i=0:num_digits-1
    b(i+1,:) = mod(floor(a./base.^i), base); % num_digits-i
end
b = b';

