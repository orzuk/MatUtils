% Sort all the elements of a structure by one field
%
% Input: 
% s - struct name
% field_name - field name 
% sort_str - ascend or descend (default ascend)
% top_k - (optional) keep only the top k elements after sorting 
%
% Output: 
% sorted_s - the struct sorted by the field 'field_nane'
%
function sorted_s = sort_struct_by_field(s, field_name, sort_str, top_k, varargin)

if(~exist('sort_str', 'var') || isempty(sort_str))
    sort_str = 'ascend';
end
eval([' cell_flag = iscell(s.' field_name ');']); % check if it's a cell array
eval(['n = size(s.' field_name ', 1);']); % get length 
if(cell_flag)
    eval(['[sorted_vec sort_perm] = sort(s.' field_name ');']);
else
    eval(['[sorted_vec sort_perm] = sort(s.' field_name ', ''' sort_str, ''');']);
end
fields = fieldnames(s); len = length(sort_perm);
sorted_s = {};
for i=1:length(fields)
    if(iscell(eval(['s.' fields{i}]))) % current field is cell array
        if(length(eval(['s.' fields{i} '{1}'])) == len)
            eval(['sorted_s.' fields{i} '{1} = vec2column(s.' fields{i} '{1});']);
            eval(['sorted_s.' fields{i} '{1} = sorted_s.' fields{i} '{1}(sort_perm,:);']);
            if(exist('top_k', 'var'))
                eval(['sorted_s.' fields{i} '{1} = sorted_s.' fields{i} '{1}(1:' num2str(top_k) ',:);']); % keep the top k with lowest scores
            end
        else
            cell_len = length(eval(['s.' fields{i}]));
            if((cell_flag) || (cell_len == len))
                eval(['sorted_s.' fields{i} ' = s.' fields{i} '(sort_perm,:);']);
            else
                eval(['sorted_s.' fields{i} ' = s.' fields{i} ';']);
            end
        end
    else % current field is a vector 
        if(length(eval(['s.' fields{i}])) == len)
            eval(['sorted_s.' fields{i} ' = vec2column(s.' fields{i} ');']);
            eval(['sorted_s.' fields{i} ' = sorted_s.' fields{i} '(sort_perm,:);']);
            if(exist('top_k', 'var'))
                eval(['sorted_s.' fields{i} ' = sorted_s.' fields{i} '(1:' num2str(top_k) ',:);']); % keep the top k with lowest scores
            end
        else
            eval(['sorted_s.' fields{i} ' = s.' fields{i} ';']);
        end
    end
end



