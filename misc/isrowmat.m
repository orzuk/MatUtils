% Check if matrix is wider than longer 
% Example: 
% XXXXXXXXXXXX
% XXXXXXXXXXXX
% XXXXXXXXXXXX
% 
function t = isrowmat(x)
t = (size(x,1) <= size(x,2));

