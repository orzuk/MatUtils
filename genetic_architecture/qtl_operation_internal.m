% Internal function for transferring genotype (data) to phenotype (qtl_vec)
% 
% Input: 
% data - set of individual liabilities (Gaussians) 
% operation - type of function  
% operation_param - function parameters
% N_vec - number(s) of different QTLs 
%
% Output: 
% qtl_vec - phenotype vector
% inds_vec - indices picked (if architecture is MAX) 
% 
function [qtl_vec inds_vec] = qtl_operation_internal(data, operation, operation_param, N_vec)

iters = size(data,1); N = size(data,2);
%N_vec = res:res:N;
qtl_vec = zeros(length(N_vec), iters);
inds_vec = ones(length(N_vec), iters);

for i=1:length(N_vec)
    cur_data = data(:,1:N_vec(i));
    switch upper(operation) % compute qtl vector. We also do the inverse operation when possible !!! 
        case 'MAX'
            [qtl_vec(i,:) inds_vec(i,:)] = max(cur_data,[],2); % take maximum
        case 'MIN'
            [qtl_vec(i,:) inds_vec(i,:)] = min(cur_data,[],2); % take minimum
        case 'SUM'
            qtl_vec(i,:) = sum(cur_data,2); % take sum
        case 'PROD'
            qtl_vec(i,:) = prod(cur_data,2); % take product
        case {'DIV', 'RATIO'}
            qtl_vec(i,:) = (cur_data(:,1)+operation_param(1)) ./ ...
                (cur_data(:,2)+operation_param(1)); % take ratio
        case 'DIFF' % take difference square (larger difference -> larger QTL phenotype)
            if(~exist('operation_param', 'var') || isempty(operation_param))
                operation_param = 2; 
            end
            qtl_vec(i,:) = (abs(cur_data(:,1) - cur_data(:,2))).^operation_param(1);             
        case 'INTERACTION' % Specify the strength of pairwise interaction 
            qtl_vec(i,:) = cur_data(:,1) + cur_data(:,1) + ...
                N_vec(i)*operation_param(1) .* cur_data(:,1).*cur_data(:,2);
        
        case 'CHI'
            qtl_vec(i,:) = sqrt(sum(cur_data.^2,2)); % take sum of squares
        case 'LP' % Take Lp norm but keep signs !!!! !
            qtl_vec(i,:) = sum(sign(cur_data).*abs(cur_data).^operation_param(1),2);
            qtl_vec(i,:) = sign(qtl_vec(i,:)) .* abs(qtl_vec(i,:).^(1/operation_param(1)));
%            qtl_vec(i,:) = (abs(sum(cur_data.^operation_param(1),2))).^(1/operation_param(1)); % take sum of squares
            
        case 'GAUSSMULTIVARIATE' % just the distance from zero 
            qtl_vec(i,:) = exp(-sum(cur_data.^2,2));             
        case 'EXP-SUM'
            qtl_vec(i,:) = log(sum(exp(operation_param(1).*cur_data),2)); % take sum of exponents
        case 'ILOGIT'
            qtl_vec(i,:) = sum(ilogit((cur_data-operation_param(1))./operation_param(2)),2); % take sum of inverse logistic        
        case 'GXE'
            env_vec = rand(size(cur_data,1),1); % environmental variance
            alpha = 5; % level of GxE interaction
            
            
            qtl_vec(i,:) = max(cur_data,[],2); % take maximum
            %        qtl_vec = max(qtl_vec, env_vec); % interaction with environment
            qtl_vec(i,:) = qtl_vec(i,:) + alpha .* qtl_vec(i,:) .* env_vec; % interaction with environment
        case 'RAND'
            qtl_vec(i,:) = rand(size(cur_data,1),1);  % take product            
    end
end
