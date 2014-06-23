% Perform unique of all elements of a structure
% by one field
% 
% Input: 
% s - the structure
% field_name - name of subfield the unique is based on
% 
% Output: 
% unique_s - the structure, with all fields as unique
%
function unique_s = unique_struct_by_field(s, field_name)

eval(['[unique_vec unique_inds] = unique(s. ' field_name ', ''rows'');']);
eval(['n = size(s.' field_name ', 1);']); % get length 
fields = fieldnames(s);
unique_s = {}; 
for i=1:length(fields) 
    eval(['field_len = size(s.' fields{i} ', 1);']); % get current fields length
    if(n == field_len)
        eval(['unique_s.' fields{i} ' = s.' fields{i} '(unique_inds,:);']);
    else
        eval(['unique_s.' fields{i} ' = s.' fields{i} ';']);
    end
end



