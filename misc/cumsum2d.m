% Perform 2-dimensional cumsum on an array 
% In the resulting array a_cum we have: 
% a_cum(i,j) = \sum_{i' <= i,j' <= j} a(i', j')
% 
% Input: 
% a - original array
% 
% Output: 
% a_cum - 2-dimensional cumulative sum 
function a_cum = cumsum2d(a) 

a_cum = cumsum(cumsum(a),2); %  + cumsum(a,2); 
