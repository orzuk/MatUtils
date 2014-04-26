% Display population names nicely
% 
% Input: 
% input_pop_str - short population name 
% 
% Output: 
% nice_pop_str - full population strings to display on figure 
% 
function nice_pop_str = get_nice_population_names(input_pop_str)

if(iscell(input_pop_str))
   
    nice_pop_str = cell(size(input_pop_str));
    for i=1:length(input_pop_str)
        nice_pop_str{i} = get_nice_population_names(input_pop_str{i});
    end
    return;    
end

switch input_pop_str
    case '2phase'
        nice_pop_str = 'Expansion2';
    case 'expan1'
        nice_pop_str = 'Expansion1';
    case 'expan2'
        nice_pop_str = 'expansion-2'; 
    case {'equil', 'equilibrium'}
        nice_pop_str = 'Equilibrium'; 
    case 'europ'
        nice_pop_str = 'Europe';
    case {'ice', 'ice3'} % new simulations from Steven 
        nice_pop_str = 'Iceland';
    case {'finn1', 'finn3'}
        nice_pop_str = 'Finland';         
    case 'finn2'
        nice_pop_str = 'finland2';         
    case {'varsel1', 'varsel2'} 
        nice_pop_str = 'variable-selection';
    case 'ice-bneck'
        nice_pop_str = 'Iceland-Bottleneck';
    case 'ice-nobneck'
        nice_pop_str = 'Iceland-NoBottleneck';        
    otherwise
        nice_pop_str = 'not-assigned'; 
        
end
