% Perform num2str but choose any delimiter rather than space 
% (e.g. one can use tab, '-', new lines etc.)
% 
% Input: 
% x - number 
% delimiter - delimiter instead of spice 
% precision - how many digits to keep 
%
function s = num2str_delim(x, delimiter, precision, varargin)

if(exist('precision', 'var'))
    s = num2str(x, precision);
else
    s = num2str(x);
end
del_inds = strfind(s, ' '); % remove consecutive occurences of delimiter (we want just one delimiter each time)
del_inds = del_inds(find(del_inds(2:end) == del_inds(1:end-1)+1)+1);

s = strrep(s(setdiff(1:length(s), del_inds)), ' ', delimiter); % replace spaces by delimiters

