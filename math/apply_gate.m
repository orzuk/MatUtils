% Apply a gate (e.g. boolean)
% We try to keep intermediate and output values as probabilities (not just binary)
%
% Input:
% x - input values
% gate_type - gate function
% gate_params - function parameters
%
% Output:
% y - gate output value
%
function y = apply_gate(x, gate_type, gate_params)

    

AssignGeneralConstants;


m = size(x,1); % number of vectors to work on
N = size(x,2); % number of loci in a vector

if(~isfield(gate_params, 'a'))
    gate_params.a = 1;
end
if(~isfield(gate_params, 'b'))
    gate_params.b = 0;
end
switch gate_type
    case AND
        y = prod(x);
    case OR
        y = 1 - prod(1-x);
    case XOR
        % y = mod(sum(x,2), 2); % binary version 
        y = ((-1).^(N-1).*prod(2.*x - 1,2) + 1) ./ 2; % prob version
    case SUM
        y = sum(x,2);
    case NOT
        y = 1-x;
    case MAX
        y = max(x);
    case MIN
        y = min(x);
    case AFFINE
        y = x .* gate_params.a + gate_params.b;
    case EXP
        y = gate_params.a .^ x;
    case MAJORITY
        y = majority(x);
    case RANDOM % output a random number (gaussian) 
        y = randn(m,1) .* gate_params.a + gate_params.b; 
    case THRESHOLD % used for liability threshold model
        y = (sign(x- gate_params.b) + 1) ./ 2; % see if larger than threshold 
    case SIGMOID
        y = sum(x .* gate_params.linear_coef_vec); % take weighted sum
        y = 0.5 * (1 + tanh(gate_params.a .* (y - gate_params.b))); % new! use tanh!!!
    case K_OF_N
        y = (sum(x,2) == gate_params.k);
    case K_OR_MORE_OF_N
        y = (sum(x,2) >= gate_params.k);
    case DOMINANT % genetic: works if at least one is on
        gate_params.k = 1;
        y = apply_gate(x, K_OR_MORE_OF_N, gate_params);
    case RECESSIVE % genetic: works if all are on
        gate_params.k = N;
        y = apply_gate(x, K_OR_MORE_OF_N, gate_params);
    case MULTIPLICATIVE
        y = apply_gate(x, ADDITIVE, gate_params); y = exp(y);
    case ADDITIVE
        y = sum(x .* repmat(vec2row(gate_params.linear_coef_vec) ,m, 1),2) + ...
            repmat(gate_params.b, m, 1);  % take weighted sum
    case LOGISTIC
        y = apply_gate(x, MULTIPLICATIVE, gate_params); y = y ./ (1+y);
    case LIABILITY % the liability threshold model
%         y = sign( sum(x .* repmat(vec2row(gate_params.linear_coef_vec) ,m, 1),2) - ...
%             gate_params.b); % b is the threshold 
%         y = (y+1)./2; % set sign to be zero/one  % BINARY version 
 
        y = sum(x .* repmat(vec2row(gate_params.linear_coef_vec) ,m, 1),2); % compute genetic part 
        h = gate_params.a; threshold = gate_params.b; % heritability (part of genetic component)
        y = 1 - normcdf((threshold - y) ./ sqrt(1-h)); % PROBABILITY Version 

end

