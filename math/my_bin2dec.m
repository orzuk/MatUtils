% Transfer a number in a binary format to decimal format
% 
% Input: 
% b - vector of digits representing the number in binary
%
% Output: 
% a - number in decimal base
%
function a=my_bin2dec(b)

num_digits = size(b, 2); 
a = sum(b .* repmat(2.^(0:num_digits-1), size(b,1), 1), 2);

