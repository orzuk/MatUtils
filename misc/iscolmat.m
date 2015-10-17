% Check if matrix is longer than wider
% Example: 
% XXXX
% XXXX
% XXXX
% XXXX
% XXXX
% XXXX
% 
function t = iscolmat(x)
t = (size(x,1) >= size(x,2));

