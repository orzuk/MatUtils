% Simulate a genetic architecture: compute trait given genotypes
%
% Input:
% x - genotypes vector. For some architectures this works also with probabilities
% architecture_str - what kind of architecture to use (e.g. additive, random etc.)
% params_struct - parameters for architecture (e.g. linear coefficients for an additive model) % z_std - extra parameter giving the variance
% num_outputs - how many outputs to simulate for EACH input
% output_type - probabilities or binary values. Optional, only for binary values
%
% Output:
% z - output trait (could be binary, continuous ...). We now use z as the prob. of z being one
% z - output clean (without final noise)
% architecture_formula - a formula of the architecture (latex format)
%
function [z z_clean architecture_formula] = ...
    genetic_architecture(x, architecture_str, ...
    params_struct, num_outputs, output_type, varargin)

AssignGeneralConstants;
AssignStatsConstants; 

N = size(x, 2); % number of loci
M = size(x, 1); % number of different x vectors
x = repmat(x, num_outputs, 1); % easy (but costly) implementation: just duplicate x

if(~isfield(params_struct, 'z_std') || isempty(params_struct.z_std))
    params_struct.z_std = 0;
end
if(~exist('output_type', 'var') || isempty(output_type))
    output_type = PROBS; 
end
architecture_formula = get_architecture_formula(architecture_str, params_struct,1);
q_z  = std_to_q_binary(params_struct.z_std); % perform convulution with a binary variable

switch architecture_str % choose architecture to use
    case {'additive', 'linear'} % simple additive model
        z = sum(x .* repmat(vec2row(params_struct.linear_coef_vec(1:N)), M*num_outputs, 1), 2); % take weighted sum
        if(isfield(params_struct, 'b'))
            z = z + params_struct.b; % add baseline affine
        end
    case 'random' % choose each value randomly
        %     params_struct.z_std = params_struct.linear_coef_vec(1); % take variance
        z_base_vec = rand(2^N, 1); % generate the 'baseline' output for each input genotypes vector
        index_vec = sum(repmat(2.^[0:N-1], M, 1) .* x(1:M,:), 2) + 1;
        z = vec2column(mat_into_vec(z_base_vec(repmat(index_vec, 1, num_outputs))));
    case 'monotone' % a random monotone function (how to force it?). Here z must be full! (doesn't work with sampling)
        z = zeros(2^N,1);
        z(1) = rand(1);
        x_vec =  dec2bin(0:2^N-1) - '0'; % take all 2^N possibilities
        ham_x_vec = sum(x_vec, 2);
        for ham = 1:N % loop on hamming distance
            ind_vec = find(ham_x_vec == ham); % find all sites with current hamming vec
            for i=1:length(ind_vec) % loop on invdividual genotype vectors x - heavy part
                x_ones = find(x_vec(ind_vec(i),end:-1:1));
                for j=1:ham % loop over all individual indices and remove them
                    z(ind_vec(i)) = max(z(ind_vec(i)), z(1+bitxor(ind_vec(i)-1, 2^(x_ones(j)-1))));
                end
            end
            z(ind_vec) = z(ind_vec) + rand(length(ind_vec), 1); % add a random component
        end
        if(size(x,1) ~= 2^N) % just use coordinates of the full function
            ind_vec = sum(repmat(2.^[0:N-1], size(x, 1), 1) .* x, 2)+1;
            z = z(ind_vec);
        end
        
    case 'xor' % take exclusive or of all genotypes
        z = mod(sum(x,2), 2);            
    case 'and' % take AND gate of all genotypes
        z = x(:,1);
        for i=2:size(x,2)
            z = z .* x(:,i); % bitand(z, x(:,i)); % AND gate
        end
    case 'or' % take OR gate of all genotypes
        z = x(:,1);
        for i=2:size(x,2)
            z =  z + x(:,i) - z .* x(:,i); %  bitor(z, x(:,i)); % OR gate
        end
    case 'multiplicative' % (generalizes 'and')
        z = exp(sum(x .* repmat(vec2row(params_struct.linear_coef_vec(1:N)), M*num_outputs, 1), 2)); % take weighted sum and multiply 
        if(isfield(params_struct, 'b'))
            z = z .* exp(params_struct.b); % mutiply by baseline affine
        end
    case {'sigmoid', 'sigmoid-additive'} % sigmoid applied to a linear sum
        z = sum(x .* repmat(vec2row(params_struct.linear_coef_vec(1:N)), M*num_outputs, 1), 2); % take weighted sum
        z = 0.5 * (1 + tanh(params_struct.a .* (z - params_struct.b))); % new! use tanh!!!   %  1 ./ (1+a.*exp(-b.*(z-c))); % apply sigmoid function
        %        a = params_vec(end-2); b = params_vec(end-1); c = params_vec(end); % set sigmoid parameters
        %        draw_sigmoid();
        
    case {'logistic', 'logit'} % logistic regression model 
        z = exp(sum(x .* repmat(vec2row(params_struct.linear_coef_vec(1:N)), M*num_outputs, 1), 2)); % take weighted sum and multiply 
        if(isfield(params_struct, 'b'))
            z = z .* exp(params_struct.b); % mutiply by baseline affine
        end
        z = z ./ (1+z); 
        
    case {'liability', 'liability-threshold', 'probit'} % threshold function 
        z = exp(sum(x .* repmat(vec2row(params_struct.linear_coef_vec(1:N)), M*num_outputs, 1), 2)); % take weighted sum and multiply
        z = 1 - phi(params_struct.b - z); % assume noise is with st.d. one! here already take environmental effect into account 
    
    case {'and-of-sigmoids', 'or-of-sigmoids', 'and-of-k-of-n', 'and-of-k_or_more_of_n'}
        arch_first_gate = str2word('-', architecture_str, 1);
        pathway_architecture_str = architecture_str_to_one_pathway_str(architecture_str);
        switch arch_first_gate
            case 'and'
                z = ones(M*num_outputs,1);
            case 'or'
                z = zeros(M*num_outputs,1); % take ands
        end
        for i=1:params_struct.num_clauses % loop on clauses
            cur_params_struct = params_struct;
            cur_params_struct.linear_coef_vec = ...
                params_struct.linear_coef_vec((i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause);
            cur_params_struct.min_freq(:) = 0; cur_params_struct.max_freq(:) = 1; % apply affine only at the end
            w = genetic_architecture( ...
                x(:,(i-1)*params_struct.k_in_clause+1:i*params_struct.k_in_clause), ...
                pathway_architecture_str, cur_params_struct, num_outputs);
            switch arch_first_gate
                case 'and'
                    z = z .* w; % AND gate for probabilities
                case 'or'
                    z = 1 - (1-z) .* (1-w); % OR gate for probabilities
            end
        end
        
        
    case {'CNF', 'DNF', 'sum-of-ands', 'sum-of-ors'} % all combinatorial polynomials
        switch architecture_str
            case 'CNF'
                z = ones(M*num_outputs,1);
            otherwise
                z = zeros(M*num_outputs,1); % take ands
        end
        for i=1:params_struct.num_clauses % loop on clauses
            switch architecture_str
                case {'DNF', 'sum-of-ands'}
                    w = ones(M*num_outputs,1); % take ands
                    for j=1:params_struct.k_in_clause % loop on variables in the clause
                        w = w .* x(:,params_struct.clause_mat(i,j)); % * replaces bitand
                    end
                case {'CNF', 'sum-of-ors'}
                    w = zeros(M*num_outputs,1); % take ands
                    for j=1:params_struct.k_in_clause % loop on variables in the clause
                        w = w + x(:,params_struct.clause_mat(i,j)) - ...
                            w .* x(:,params_struct.clause_mat(i,j)); %   bitor(w, x(:,params_struct.clause_mat(i,j))); % should replace bitor
                    end
            end
            switch architecture_str
                case 'DNF',
                    z = z + w - z .* w; % bitor(z,w); % OR gate
                case  'CNF'
                    z = z .* w; % bitand(z,w); % AND gate
                case {'sum-of-ands',  'sum-of-ors'}
                    z = z + w; % sum
                    if(i == params_struct.num_clauses) % normalize to get a probability
                        z = z ./ max(z);
                    end
            end
        end % loop on clauses

    case {'k-of-n', 'k_or_more_of_n'} % binomial architecture (k of n)
        gate_params.k = params_struct.a;
        z = apply_gate(x, K_OR_MORE_OF_N, gate_params);
        
    case 'given' % just use the input parameters to fix the architecture
        z = vec2column(params_vec);
        
        
    case 'circuit' % New: here the architecture is given by a circuit with all possible
        z = zeros(M,1); 
        for i=1:M % circuit not vectorized yet :-(
            z(i) = apply_circuit(x(i,:), params_struct.circuit, ...
                params_struct.circuit_gate_params, num_outputs);
        end
        % gates (most general form!!!)
end

z_clean = z; % noiseless version
arch_type = architecture_to_type(architecture_str);
switch arch_type
    case BINARY % in the binary case the output is the PROB of one!!
        if(isfield(params_struct, 'min_freq'))
            z = z .* (params_struct.max_freq - params_struct.min_freq) + ...
                params_struct.min_freq; % apply general affine transfomration
        else
            z = z .* (1-q_z) + (1-z) .* q_z; % apply convolution with q
        end
        if(output_type == BINARY) % here report actual values
           z = z > rand(size(z)); 
        end
    case CONTINUOUS
        z = z + randn(M*num_outputs, 1) .* params_struct.z_std; % add extra variance with Gaussian noise
end


% function  draw_sigmoid()
%                %         if(size(x,1) == 2^N)
%             x_vec = [-10:0.01:10]; y_vec = 0.5 * (1 + tanh(a .* (x_vec - b))); % 1 ./ (a+b.*exp(-c.*x_vec));
% %            y_diff_vec = b*c.* exp(-c.*x_vec) ./ (a+b.*exp(-c.*x_vec)).^2;
%             %            figure; plot(x_vec, y_diff_vec); title('y derivative');
%             figure; hold on; plot(x_vec, y_vec);
%             plot(sum(x .* repmat(vec2row(params_vec(1:N)), M, 1), 2), z, 'r*'); % a.*z
%
% %             my_saveas(gcf, ['../../common_disease_model/figures/sigmoid_example_N_' ...
% %                 num2str(N) '_a_' num2str(a, 2) '_b_' num2str(b, 2)], ... % no need for c  '_c_' num2str(c, 2)], ...
% %                 {'jpg', 'epsc', 'fig'});
%         end
%        z = z .* a; % make it a probability
%        % z = double(rand(size(z,1),1) > a.*z); % last stage: randomize z according to probabilities - wrong !!!!




