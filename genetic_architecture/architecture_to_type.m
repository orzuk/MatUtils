% Determine the type of an architecture (e.g. binary, discrete etc.)
function arch_type = architecture_to_type(architecture_str)

AssignStatsConstants;

switch architecture_str
    
    case {'and', 'or', 'xor',  'sigmoid-additive', 'multiplicative', 'logistic', ...
            'sigmoid', 'and-of-sigmoids', 'or-of-sigmoids', 'liability', 'and-of-k-of-n', 'and-of-k_or_more_of_n', ...
            'k-of-n', 'k_or_more_of_n', ...
            'CNF', 'DNF', 'sum-of-ands', 'sum-of-ors', 'dominant', 'recessive',  'circuit'}
        arch_type = BINARY;
        
    case {'additive',  'monotone', 'random', 'given', 'linear'}
        arch_type = CONTINUOUS;
end


