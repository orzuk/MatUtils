% Check if array is sorted with tolerance. 
% We check each line seperately 
% Input: 
% x - array
% tol - tlerance for checking sorted
% direction - increasing (default) or decreasing 
%
% Output: 
% ret - 
% 
function ret = issortedtol(x, tol, direction)

if(~exist('tol', 'var') || isempty(tol))
    tol = 10e-10;
end
if(~exist('direction', 'var') || isempty(direction))
    direction = 'increasing';
end
if(strcmp(direction, 'decreasing'))
    ret = issortedtol(-x, tol); 
else
    ret = 1;
    for i=1:size(x, 1)
        if(min(diff(x(i,:))) < -tol) % find unsorted line 
            bad_i = i
            bad_tol = min(diff(x(i,:)))
            ret = 0; return; 
        end
    end    
end
    
    
