% Find for each element of a the closest element of b
%
% Input:
% a - first vector (query)
% b - second vector (how close to it)
% in_c_flag - use fast c implementation (default)
%
% Output:
% closest_ind - inds of closest elements in b to a
% closest_dist - distance of closest elements in b to a
%
function [closest_ind closest_dist] = match_closest_vals(a, b, in_c_flag, varargin)

if(isempty(b))
    closest_ind = []; closest_dist = []; return; % nothing to compare to
end
[a_sorted sort_perm_a] = sort(a); n = length(a);
[b_sorted sort_perm_b] = sort(b); m = length(b);
inv_sort_perm_a = inv_perm(sort_perm_a);
% inv_sort_perm_b = inv_perm(sort_perm_b);
if(~exist('in_c_flag', 'var')) % default is using the fast C code 
    in_c_flag = 1;
end
if(in_c_flag) % call function in c (fast)
%     ab = [a_sorted b_sorted]; % new implementation: sort everything 
%     [ab_sorted sort_perm_ab] = sort(ab); % sort the joint vector 
%     a_inds = sort_perm_ab(1:n); b_inds = sort_perm_ab(n+1:end); 
%     is_a_vec = zeros(m+n,1); is_a_vec(a_inds) = 1; 
    [closest_ind, closest_dist] = match_closest_vals_c(double(a_sorted), double(b_sorted)); % Convert to double (in the future might want to enable single support in c)
else% do everything in matlab (slow)
    closest_ind = zeros(n,1); closest_ind = zeros(n,1);
    ctr = 1;
    for i=1:n
        f = find(a_sorted(i) < b_sorted(ctr:end), 1);
        if(isempty(f)) % this means that this value of a is bigger than all values of b
            closest_ind(i:n) = m;
            closest_dist(i:n) = a_sorted(i:n) - b_sorted(m);
            break; % get out of loop
        else
            if(f == 1)
                closest_ind(i) = ctr;
            else
                [dummy closest_ind(i)] = min(abs(b_sorted(ctr+f-2:ctr+f-1) - a_sorted(i)));
                closest_ind(i) =  closest_ind(i) + f-2 + ctr-1;
            end
            closest_dist(i) = a_sorted(i) - b_sorted(closest_ind(i));
            ctr = closest_ind(i);
        end
        f_is = f;
        ctr_is = ctr;
    end
end
    
% apply inverse permutations to return to original
closest_ind = sort_perm_b(closest_ind(inv_sort_perm_a));
closest_dist = closest_dist(inv_sort_perm_a);


    
