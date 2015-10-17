% Perform union of the elements of two structures
% by one field (so if we have many times the same elements in another field
% it's kept duplicated)
% 
% Input: 
% s1 - the first structure
% s2 - the second structure
% field_name - name of subfield the unique is based on
% 
% Output: 
% union_s - the structure, with all fields as union
%
function union_s = union_struct_by_field(s1, s2, field_name)

eval(['len1 = length(s1.' field_name ');']); 
eval(['len2 = length(s2.' field_name ');']);  % get two lengths 

eval(['[union_vec union_inds1 union_inds2] = union(s1.' field_name ', s2.', field_name, ',  ''rows'');']);
fields = fieldnames(s1);
%union_inds1 = union_inds(union_inds <= len1);    % seperate inds one from inds two
%union_inds2 = union_inds(union_inds > len1) - len1;

union_s = {};
for i=1:length(fields)
    if(length(eval(['s1.' fields{i}])) == len1)
        eval(['s1.' fields{i} ' = vec2row(s1.' fields{i} ');']);
        eval(['s2.' fields{i} ' = vec2row(s2.' fields{i} ');']);
        eval(['union_s.' fields{i} ' = vec2column([ s1.' fields{i} ...
            '(:,union_inds1)  s2.' fields{i} '(:,union_inds2) ]);']);
    else
        if(iscell(eval(['s1.' fields{i}])))
            if(length(eval(['s1.' fields{i} '{1}'])) == len1)
                eval(['union_s.' fields{i} '{1} = vec2column([ vec2row(s1.' fields{i} ...
                    '{1}(union_inds1,:))  vec2row(s2.' fields{i} '{1}(union_inds2,:)) ]);']); % always return a column vector
            else
                eval(['union_s.' fields{i} ' = s1.' fields{i} ';']); % just copy one of them 
            end
        else
                eval(['union_s.' fields{i} ' = s1.' fields{i} ';']); % just copy one of them
        end
    end
    %    eval(['unique_s.' fields{i} ' = s.' fields{i} '(unique_inds);']);
end



