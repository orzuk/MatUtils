% Take the string of an architecture pathway from the string for the entire architecture
% 
% Input: 
% 
function one_pathway_str = architecture_str_to_one_pathway_str(architecture_str)

switch architecture_str
    case {'and-of-ors', 'sum-of-ors', 'CNF'}
        one_pathway_str = 'or';
    case {'or-of-ands', 'sum-of-ands', 'DNF'}
        one_pathway_str = 'and';
    case {'and-of-sigmoids', 'or-of-sigmoids', 'sigmoid'}
        one_pathway_str = 'sigmoid';
    otherwise
        delim_pos = find(architecture_str == '-');
        if(~isempty(delim_pos))
            one_pathway_str = architecture_str(delim_pos(2)+1:end);
        else
            one_pathway_str = architecture_str;
        end
end



