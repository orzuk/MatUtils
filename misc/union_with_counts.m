% Perform union of two lists, take into account also the counts of the lists.
%
% Input: 
% inds1 - the indices of the first set
% counts1 - the counts of the first set
% inds2 - the indices of the second set
% counts2 - the counts of the second set
%
% Output:
% union_inds - the indices of the union
% union_counts - the counts of the union
%
function [union_inds union_counts] = union_with_counts(inds1, counts1, inds2, counts2)

if(isempty(inds1)) % if first is empty -> just take the second
    union_inds = inds2; union_counts = counts2;
else
    union_inds = inds1; union_counts = counts1;
    if(~isempty(inds2)) % nothing to do if the second one is empty
        if(size(inds1, 2) == 1)
            [intersect_inds I1 I2] = intersect(inds1, inds2); % common inds
        else
            [intersect_inds I1 I2] = intersect(inds1, inds2, 'rows'); % common inds
        end
        if(~isempty(intersect_inds))
            if(size(counts1,2) == 1)
                union_counts(I1) = counts1(I1) + counts2(I2);
            else
                union_counts(I1,:) = counts1(I1,:) + counts2(I2,:);
            end
        end
        if(size(inds2, 1) == 1)
            [new_inds I2] = setdiff(inds2, inds1); % new inds
        else
            [new_inds I2] = setdiff(inds2, inds1, 'rows'); % new inds
        end
        if(~isempty(new_inds))
            union_inds(end+1:end+size(new_inds,1),:) = new_inds;
            if(size(counts2,2) == 1)
                union_counts(end+1:end+size(new_inds,1)) = counts2(I2);
            else
                union_counts(end+1:end+size(new_inds,1),:) = counts2(I2,:);
            end
        end
    end
end

if(size(union_counts,2) == 1)   % sort by counts again 
    [union_counts sort_perm] = sort(union_counts, 'descend'); 
else % sort by the first counts positions
    [dummy sort_perm] = sort(union_counts(:,1), 'descend');  union_counts = union_counts(sort_perm,:);
end
union_inds = union_inds(sort_perm,:);




