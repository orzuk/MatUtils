% Convert a sting to all the numbers (doubles) contained in it,
% assumed to be seperated by anything other then '.'
%
% Input:
% s - string with numbers
% ind - take only a certain number out of all (optional)
% special_parse_flag - enable also converting percantage to numbers, and e-6 to 10^(-6) (default: false)
% do_arithmatics - compute simple numeric operations automatically (default: yes)
%
% Output:
% nums - numbers appearing in string
%
function nums = str2nums(s, ind, special_parse_flag, do_arithmatics, varargin)

if(isempty(s)) % return empty array
    nums = []; return;
end
tmp = regexp(s, '[0-9]');
if(isempty(tmp))
    nums = []; return;
end
if(~exist('do_arithmatics', 'var') || isempty(do_arithmatics))
    do_arithmatics = 1; 
end
if(do_arithmatics)
    dots = regexp(s, '[\/.*+]'); % enable all arithmatic operations except minus (minus sign can be used to set intervals)
else  % don't do arithmatics
    dots = []; 
end
if(~exist('special_parse_flag', 'var') || isempty(special_parse_flag))
    special_parse_flag = 0;
end

seperators = tmp(diff(tmp) > 1)+1; % these seperate between different numbers
dot_seperators = intersect( tmp(diff(tmp) == 2) + 1, dots); % don't seperate over dots
if(special_parse_flag)
    exps = strfind(lower(s), 'e'); % numbers represented in exponential form
    exps_seperators = intersect( tmp(diff(tmp) == 2) + 1, exps);
    exps = strfind(lower(s), 'e-'); % numbers represented in exponential form (negative exponent)
    neg_exp_seperators = intersect( tmp(diff(tmp) == 3) + 1, exps);
    exps_seperators = union(exps_seperators, [neg_exp_seperators neg_exp_seperators+1]); %intersect( tmp(diff(tmp) == 3) + 1, exps));
    %    dot_seperators = union(dot_seperators, exps_seperators);
    seperators = setdiff(seperators, exps_seperators);
    tmp = union(tmp, exps_seperators);
end
seperators = setdiff(seperators, dot_seperators);
seperators = [seperators length(s)+1];
tmp = union(tmp, dot_seperators);

n = length(seperators);
nums = zeros(1,n);  prev_sep = 0;
[good_dots, good_dots_inds] = intersect(dots+1, tmp); % find symbols like '.' which are within an number
tmp = union(tmp, good_dots-1); % add these to tmp 

for i=1:n
    inds = find( (tmp < seperators(i)) & (tmp >= prev_sep) );
    prev_sep = seperators(i);
    if(isempty(inds))
        nums(i) = [];
    else
        if(ismember(s(tmp(inds(1))), '/\*')) % avoid division signs at the beginning
            inds = inds(2:end); 
        end
        nums(i) = str2num(s(tmp(inds)));
    end
    if(s(max(1,tmp(inds(1))-1)) == '-') % enable negative numbers
        nums(i) = -nums(i);
    end
end


if(special_parse_flag) % Deal with special percentage and exponents
    percent_inds = strfind(s, '%');
    for i=1:length(percent_inds)
        cur_ind = find(percent_inds(i) == [seperators(1:end-1), max(tmp)+1], 1); % adjust to catch also last separator
        if(~isempty(cur_ind))
            nums(cur_ind) = nums(cur_ind)/100; % convert percentage to 0-1 number
        end
    end
end

if(exist('ind', 'var') && (~isempty(ind)))
    nums = nums(ind);
end






