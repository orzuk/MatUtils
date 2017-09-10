% Compute a string representing the architecutre formula. 
% If formula is too long we split it to differenct lines 
% 
% Input: 
% architecture_str - string representing architecture
% params_struct - structure with architecture parameters
% split_formula_to_lines - whether to split a long formula to different lines
% 
% Output: 
% architecture_formula - string with architecture formula
% 
function architecture_formula = ...
    get_architecture_formula(architecture_str, params_struct, split_formula_to_lines)

if(~exist('split_formula_to_lines', 'var') || isempty(split_formula_to_lines))
    split_formula_to_lines = 0;
end
iters = length(params_struct.z_std);
if(iters > 1)
    architecture_formula = cell(iters,1); 
    for i=1:iters
        cur_params_struct = params_struct;
        cur_params_struct.z_std = params_struct.z_std(i);
        if(isfield(cur_params_struct, 'min_freq'))
            cur_params_struct.min_freq = params_struct.min_freq(i);
            cur_params_struct.max_freq = params_struct.max_freq(i);
        end
        architecture_formula{i} = get_architecture_formula(architecture_str, ...
            cur_params_struct, split_formula_to_lines);
    end
    return; 
end

precision = 3;
% q_z  = std_to_q_binary(params_struct.z_std); % (1-sqrt(1-4*params_struct.z_std^2))/2; % perform convulution with a binary variable
architecture_formula = '\be Pr(z=1) = ';

if(isfield(params_struct, 'min_freq') && ((params_struct.min_freq > 0) || (params_struct.max_freq < 1)))
    architecture_formula = [architecture_formula num2str(params_struct.min_freq, precision) ' + ' ...
        num2str( (params_struct.max_freq-params_struct.min_freq), precision) '\Big\{ '];
else % still add parentethis 
    architecture_formula = [architecture_formula '\Big\{ '];
end

% if(params_struct.z_std > 0) % convolution: x*(1-q) + (1-x)*q
%     architecture_formula = [architecture_formula num2str(q_z, precision) ' + ' ...
%         num2str(1-2*q_z, precision) '\Big\{ '];
% end
if(isfield(params_struct, 'N'))
    N = params_struct.N; 
else
    if(isfield(params_struct, 'linear_coef_vec'))
        N = length(params_struct.linear_coef_vec);
    else
        N = -1 % We don't know what's the number of loci!!!! 
    end
end

switch architecture_str % choose architecture to use
    case {'additive', 'linear'} % simple additive model
        for i=1:N
            architecture_formula = [architecture_formula '+' ...
                num2str(params_struct.linear_coef_vec(i), precision) i_to_x_str(i)]; % 'x_{' num2str(i) '}'];
        end
    case 'random' % choose each value randomly
    case 'monotone' % a random monotone function (how to force it?). Here z must be full! (doesn't work with sampling)
    case 'xor' % take exclusive or of all genotypes
    case 'and' % take and gate of all genotypes
    case 'or' % take or gate of all genotypes
    case {'sigmoid', 'sigmoid-additive'} % sigmoid applied to a linear sum
        architecture_formula = [architecture_formula ...
            '0.5 \{1 + \tanh[' num2str(params_struct.a, precision) '('];
        for i=1:N
            if(i > 1)
                plus_str = '+';
            else
                plus_str = '';
            end
            if(params_struct.linear_coef_vec(i) ~= 1)
                num_str = num2str(params_struct.linear_coef_vec(i), precision);
            else
                num_str = '';
            end
            architecture_formula = ...
                [architecture_formula plus_str num_str i_to_x_str(i)]; % 'x_{' num2str(i) '}'];
        end
        architecture_formula = ...
            [architecture_formula ' - ' num2str(params_struct.b, precision) ')]\}'];
        
    case {'and-of-sigmoids', 'or-of-sigmoids'}
%        architecture_formula = '';
        for i=1:params_struct.num_clauses % loop on clauses
            cur_params_struct = params_struct;
            cur_params_struct.linear_coef_vec = params_struct.linear_coef_vec( ...
                (i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause );
            cur_params_struct.min_freq(:) = 0; cur_params_struct.max_freq(:) = 1; 
            cur_architecture_formula = ...
                get_architecture_formula('sigmoid', cur_params_struct);
            if(i > 1)
                switch architecture_str
                    case 'and-of-sigmoids'
                        architecture_formula = [architecture_formula ' \wedge '];
                    case 'or-of-sigmoids'
                        architecture_formula = [architecture_formula ' \vee '];
                end
            else
%                architecture_formula = [strdiff(architecture_formula, 'Pr(z=1) ') ...
%                    strdiff(strdiff(cur_architecture_formula(1:14), '\be')];
            end
            architecture_formula = [architecture_formula ...
                ' \bigg\{ ' cur_architecture_formula(15:end-3) ' \bigg\} '];
        end
%        architecture_formula = [architecture_formula ' \ee']; % this is done once in the and
                
    case {'CNF', 'DNF', 'sum-of-ands', 'sum-of-ors'} % all combinatorial polynomials
        for i=1:params_struct.num_clauses % loop on clauses
            if(i == 1)
                architecture_formula = [architecture_formula ' ('];
            else
                switch architecture_str
                    case 'DNF',
                        architecture_formula = [architecture_formula ' \vee ('];
                    case  'CNF'
                        architecture_formula = [architecture_formula ' \wedge ('];
                    case {'sum-of-ands',  'sum-of-ors'}
                        architecture_formula = [architecture_formula ' + \frac{1}{' ...
                            num2str(params_struct.num_clauses) '} ('];
                end
            end
            switch architecture_str
                case {'DNF', 'sum-of-ands'}
                    for j=1:params_struct.k_in_clause % loop on variables in the clause
                        if(j < params_struct.k_in_clause)
                            wedge_str = ' \wedge ';
                        else
                            wedge_str = '';
                        end
                        architecture_formula = ...
                            [architecture_formula i_to_x_str(params_struct.clause_mat(i,j)) ...
                            wedge_str]; % wedge=AND, vee=OR
                        %                            'x_{' num2str(params_struct.clause_mat(i,j)) '}'
                    end
                case {'CNF', 'sum-of-ors'}
                    for j=1:params_struct.k_in_clause % loop on variables in the clause
                        if(j < params_struct.k_in_clause)
                            vee_str = ' \vee ';
                        else
                            vee_str = '';
                        end
                        architecture_formula = ...
                            [architecture_formula i_to_x_str(params_struct.clause_mat(i,j)) ...
                            vee_str]; % wedge=AND, vee=OR
%                            'x_{' num2str(params_struct.clause_mat(i,j)) '}' 
                    end
            end
            architecture_formula = [architecture_formula ')'];
        end % loop on clauses
    case 'given' % just use the input parameters to fix the architecture
end
% if(params_struct.z_std > 0) % we don't use this one anymore ... 
     architecture_formula = [architecture_formula ' \Big\} '];
% end
architecture_formula = [architecture_formula  ' \ee'];

if(split_formula_to_lines)
    max_chr_in_line = 100; % Split to lines
    len = length(architecture_formula);
    new_line_pos = max_chr_in_line:max_chr_in_line:len;
    for i=1:length(new_line_pos)
        tmp = find(architecture_formula(new_line_pos(i):end) == ' ', 1) + new_line_pos(i)-1;
        if(isempty(tmp))
            new_line_pos = new_line_pos(1:i-1);
            break;
        else
            new_line_pos(i) = tmp;
        end
    end
    if(~isempty(new_line_pos))
        tmp_formula =  ['$$ ' architecture_formula(4:new_line_pos(1))];
        for i=2:length(new_line_pos)
            tmp_formula = [tmp_formula ' $$ $$ ' architecture_formula(new_line_pos(i-1)+1:new_line_pos(i))];
        end
        tmp_formula = [tmp_formula ' $$ \be ' architecture_formula(new_line_pos(end)+1:end)];
        architecture_formula = tmp_formula;
    end
end


function x_str = i_to_x_str(i)

if(mod(i,2) == 1) % maternal SNP
    x_str = ['x_{' num2str((i+1)/2) '}'];
else % paternal SNP
    x_str = ['y_{' num2str(i/2) '}'];
end
