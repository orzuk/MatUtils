% Estimate IBD sharing for each pair of individuals from SNP data
%
% Input:
% snp_mat - matrix of snp values for all individuals
% f_vec - vector of allele frequencies for all snps
% constant_flag - flag choosing the esitmator type (how to combine evidence from different SNPs) 
% 
% Ouptut:
% IBD_mat - matrix of IBD sharing between each pair of individuals.
%
function IBD_mat = estimate_IBD_sharing_from_snps(snp_mat, f_vec, constant_flag)


if(~exist('constant_flag', 'var') || isempty(constant_flag))
    constant_flag=0;
end
% There is no closed form solution :( - we solve numerically
% which is rather slow for N^2 people

[n, m] = size(snp_mat); % m - number of individuals, n - number of SNPs

compact_flag = 1; 
if(compact_flag) 
    [unique_f_vec, f_inds, unique_f_counts] = unique_with_inds(f_vec);
%     unique_f_counts = length_cell(f_inds)
%     [unique_f_vec, unique_f_counts] = unique_with_counts(f_vec);
    num_unique_f = length(unique_f_vec);
    xy00_nums = zeros(1, num_unique_f); 
    xy11_nums = zeros(1, num_unique_f); 
end
    
IBD_mat = zeros(m);
for i=1:m
    if(mod(i,20)==0)
        run_i = i
    end
    for j=(i+1):m % loop on pairs of individuals
        xy_vec = snp_mat(:,i)+snp_mat(:,j);
        %        IBD_mat(i,j) = fzero(@(k) LL_IBD_derviative_internal(k, xy_vec, vec2column(f_vec)), 0.5);
        if(compact_flag) 
            xy01_nums = sum(xy_vec == 1);             
            for k=1:num_unique_f
                xy00_nums(k) = sum(xy_vec(f_inds{k}) == 0);
                xy11_nums(k) = sum(xy_vec(f_inds{k}) == 2);
            end
%            xy11_nums = unique_f_counts - xy00_nums;
        end
        switch constant_flag
            case 0 % MLE
                if(~compact_flag)
                    IBD_mat(i,j) = fminbnd(@(k) -LL_IBD_internal(k, xy_vec, vec2column(f_vec)), 0, 1);
                else
                    IBD_mat(i,j) = fminbnd(@(k) -LL_IBD_internal_compact(k, ...
                        xy00_nums, xy01_nums, xy11_nums, vec2row(unique_f_vec)), 0, 1);
                end
                    
            case 1 % here all f_vecs are the same
                n00 = sum(xy_vec==0);
                n01 = sum(xy_vec==1);
                n11 = sum(xy_vec==2);
                f=f_vec(1);
                A = n*f*(1-f);
                B = n00*(f^2-f*(1-f)) + n01*(f^2+(1-f)^2) + n11*((1-f)^2-f*(1-f));
                C = n00*(-f^2) + n01*f*(1-f) - n11*(1-f)^2;
                IBD_mat(i,j) = (-B+sqrt(B^2-4*A*C))/ (2*A);
            case 2 % Eric's estimator (weighted average)
                p_equal_vec = f_vec.^2+(1-f_vec).^2;
                equal_vec = 1-mod(xy_vec,2);
                V_vec = 1-2.*f_vec .*(1-f_vec);
                w_vec = 1 ./ V_vec; w_vec = w_vec ./ sum(w_vec);
                IBD_mat(i,j) = sum(w_vec .* (equal_vec - p_equal_vec) ./ (1-p_equal_vec));
                %                IBD_mat(i,j) = mean( (equal_vec - p_equal_vec) ./ (1-p_equal_vec));
                
            case 3 % Eric's estimator (simple average)
                p_equal_vec = f_vec.^2+(1-f_vec).^2;
                equal_vec = 1-mod(xy_vec,2);
                V_vec = 1-2.*f_vec .*(1-f_vec);
                w_vec = 1 ./ V_vec; w_vec = w_vec ./ sum(w_vec);
                %                IBD_mat(i,j) = sum(w_vec .* (equal_vec - p_equal_vec) ./ (1-p_equal_vec));
                IBD_mat(i,j) = mean( (equal_vec - p_equal_vec) ./ (1-p_equal_vec));
            case 4 % new estimator - what's that? 
                n00_vec = (xy_vec==0);
                n01_vec = (xy_vec==1);
                n11_vec = (xy_vec==2);
                
                IBD_mat(i,j) = mean( (n00_vec - (1-f_vec).^2) ./ (f_vec.*(1-f_vec)) - ...
                    (n01_vec - 2.*f_vec.*(1-f_vec)) ./ (2.*f_vec.*(1-f_vec)) + ...
                    (n11_vec - f_vec.^2) ./ (f_vec.*(1-f_vec)) ) / 3;
                
            case 5 % Visscher genetic relationship matrix
                G = genetic_relationship_matrix(snp_mat, 'binary', f_vec);
                IBD_mat(i,j) = G(1,2);
        end % switch constant
        
    end % loop on j
end % loop on i
IBD_mat = IBD_mat+IBD_mat'+eye(m); % -diag(diag(IBD_mat));



% Loglikelihood of data given k
%
% Input:
% k - the IBD sharing parameter
% xy_vec - vector denoting jointly x and y (actually their sum)
% f_vec - allele frequency of SNPs
%
% Output:
% ret - the log-likelihood at k
%

function ret = LL_IBD_internal(k, xy_vec, f_vec)

ret = sum(log(k.*f_vec+1-f_vec) .* (xy_vec == 0)) + ... % 00
    sum(log(1-k) .* (xy_vec == 1)) + ... % 01 or 10
    sum(log(-k.*f_vec+k+f_vec) .* (xy_vec == 2)); % 11


% Compute likelihood in a compact way
function ret = LL_IBD_internal_compact(k, xy00_nums, xy01_nums, xy11_nums, unique_f_vec)

ret = xy01_nums .* log(1-k) + ...
    sum(xy00_nums .* log(k.*unique_f_vec+1-unique_f_vec)) +  ...
    sum(xy11_nums .* log(k.*(1-unique_f_vec)+unique_f_vec)); 



% Derivative of log-likelihood dLL(x,y; k, f) / dk
% Mainly used to solve dLL(x,y; k, f) / dk = 0 and get k maximizing LL
%
% Input:
% k - the IBD sharing parameter
% xy_vec - vector denoting jointly x and y (actually their sum)
% f_vec - allele frequency of SNPs
%
% Output:
% ret - the derivative at k
%
function ret = LL_IBD_derviative_internal(k, xy_vec, f_vec)

ret = sum((f_vec ./ (k.*f_vec+1-f_vec)) .* (xy_vec == 0)) + ... % 00
    sum((1 ./ (k-1)) .* (xy_vec == 1)) + ... % 01 or 10
    sum(((1-f_vec) ./ (-k.*f_vec+k+f_vec)) .* (xy_vec == 2)); % 11



