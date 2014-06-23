% Take only certain inds of each field in a struct
%
% Input:
% s - the structure
% inds_vec - indices to extract for each field
%
% Output:
% new_s - the structure, with all fields taken only specific inds
%
function new_s = struct_by_inds(s, inds_vec)

%eval(['[unique_vec unique_inds] = unique(s. ' field_name ', ''rows'');']);
% eval(['n = size(s.' field_name ', 1);']); % get length
field_names = fieldnames(s);
new_s = {};
for i=1:length(field_names)
    i_is = i
    field_is = field_names{i}
    %    eval(['field_len = size(s.' fields{i} ', 1);']); % get current fields length
    %    if(n == field_len)
    eval(['new_s.' field_names{i} ' = s.' field_names{i} '(inds_vec,:);']);
    %    else
    %        eval(['unique_s.' fields{i} ' = s.' fields{i} ';']);
    %    end
end


