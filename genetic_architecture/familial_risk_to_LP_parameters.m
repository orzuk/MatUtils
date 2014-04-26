% Fit parameters of the LP model to observed familial risk/correlations
%
% Input:
% lambda_or_qtl_R - relative risk for relative of degree R
% N - number of liabilities
% k_R - how close are relatives
% trait_type - disease or quantitative
% mu - disease prevalence
% c_R_vec - vector of coefficients for shared environment (of the total c_R)
%
% Output:
% h_x_one_liab - heritability of each liability
% c_R - shared environment
%
function [h_x_one_liab c_R] = ...
    familial_risk_to_LP_parameters(lambda_or_qtl_R, N,  k_R, trait_type, mu, c_R_vec)

m = length(lambda_or_qtl_R); % number of parameters
options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance


%if(m==1) % here only lambda_mz
twin_str = {'MZ', 'DZ'};
for i=1:m
    switch trait_type
        case 'quantitative' % match r_MZ
            h_x_one_liab(i) = fminbnd(@(x) ...
                abs(qtl_familial_correlations_internal(N, N, k_R, 'MAX', ...
                [], x, 0, 'numeric', [], 'MZ')-lambda_or_qtl_R(i)), ...
                lambda_or_qtl_R(i), 1, options); % we know that h_z < h_x < 1 but h_pop can be MORE than h_x !!!!
        case 'binary'
            mu_l = 1-(1-mu).^(1/N);
            h_x_one_liab(i) = fminbnd(@(x) ...
                abs(compute_lambda_R_LP_internal(N, 1, x, 0, k_R(i), ...
                mu_l, mu)-lambda_or_qtl_R(i)), ...
                0, 1/k_R(i), options);
    end
end % loop on m
if(m==1)
    c_R=0; h_x_one_liab = min(1, h_x_one_liab);
else
    r_R = h_x_one_liab .* k_R; % compute familial correlations
    
    if(~exist('c_R_vec', 'var') || isempty(c_R_vec)) % here let all relationships have shared environment
        h_x_one_liab= (r_R(1)-r_R(2)) / (k_R(1)-k_R(2));
        h_x_one_liab = min(1, h_x_one_liab);
        c_R = r_R(2) - k_R(2) * h_x_one_liab;
        if(c_R<0) % can't have this !!!
            h_x_one_liab =  h_x_one_liab+c_R.*2; % (1-h_x_one_liab)*2; % soak up c_R into h_x_liab
            c_R=0;
        end
        if(h_x_one_liab==1)
            c_R=1; % assign 'full power' to shared environment
        else
            c_R = c_R / (1-h_x_one_liab); % get relative fraction of environment shared 
            c_R = min(c_R,1);
        end
        
    else % here allow different shared environment for different familial relationship
        alpha = c_R_vec(2)/c_R_vec(1);
        c_R = (k_R(2)*r_R(1)-k_R(1)*r_R(2)) / (k_R(2)-k_R(1)*alpha);
        if(c_R<0) % can't have this !!!
            c_R=0;
        end
        h_x_one_liab = (r_R(1)-c_R) / k_R(1);
        h_x_one_liab = min(1, h_x_one_liab);
        if(h_x_one_liab==1)
            c_R=1; % assign 'full power' to shared environment
        else
            c_R = min(c_R,1);
        end    
    end    
    %     c_R = 2*h_x_one_liab(2)-h_x_one_liab(1);
    %     h_x_one_liab=2*h_x_one_liab(1)-2*h_x_one_liab(2);
end

% if m==1
%
%    else % use two parameters to fit both h_x and c_R


