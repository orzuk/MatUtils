% Get indices of important class of alleles: missense, synonymous and stop
% Input: 
% S - structure 
% exome_struct - structure with exome information
% 
% Output: 
% good_allele_inds - indices of good alleles
%
function good_allele_inds = get_good_allele_inds(S, exome_struct)

good_allele_inds = cell(4,1);
switch exome_struct.data_str
    case 'ESP'
        good_allele_inds{1} = strfind_cell(lower(S.allele_types), 'syno');
        good_allele_inds{2} = strfind_cell(lower(S.allele_types), 'missen');
        good_allele_inds{3} = strfind_cell(lower(S.allele_types), 'stop-gained');
        for i=1:3
            good_allele_inds{i} = sort(setdiff(good_allele_inds{i}, strfind_cell(lower(S.allele_types), 'splice')));
            good_allele_inds{4} = sort(union(good_allele_inds{4}, good_allele_inds{i}));
        end
        
        % %         good_allele_inds = union(strfind_cell(lower(S.allele_types), 'syno'), strfind_cell(lower(S.allele_types), 'missen'));
        % %         good_allele_inds = union(good_allele_inds, strfind_cell(lower(S.allele_types), 'stop-gained'));
        % %         % S.good_allele_inds = union(good_allele_inds, strfind_cell(lower(allele_types), 'coding-notmod3')); % NEW! Add frameshifts!!!
        % %         good_allele_inds = sort(setdiff(good_allele_inds, strfind_cell(lower(S.allele_types), 'splice'))); % set which types of alleles to plot
        
    case 'ExAC'
        good_allele_inds{1} = strfind_cell(lower(S.allele_types), 'syno');
        good_allele_inds{2} = strfind_cell(lower(S.allele_types), 'missen');
        good_allele_inds{3} = strfind_cell(str2title(lower(S.allele_types)), 'stop-gained');
        for i=1:3
            good_allele_inds{i} = sort(setdiff(good_allele_inds{i}, strfind_cell(lower(S.allele_types), 'splice')));
            good_allele_inds{4} = sort(union(good_allele_inds{4}, good_allele_inds{i}));
        end
        
        % %         good_allele_inds = union(strfind_cell(lower(S.allele_types), 'syno'), strfind_cell(lower(S.allele_types), 'missen'));
        % %         good_allele_inds = union(good_allele_inds, strfind_cell(str2title(lower(S.allele_types)), 'stop-gained'));
        % %         % S.good_allele_inds = union(good_allele_inds, strfind_cell(lower(allele_types), 'coding-notmod3')); % NEW! Add frameshifts!!!
        % %         good_allele_inds = sort(setdiff(good_allele_inds, strfind_cell(lower(S.allele_types), 'splice'))); % set which types of alleles to plot
end

