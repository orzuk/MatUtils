% Simulate linear regression model and quasi regression
%
% Input:
% n_vec - vector of number of samples
% d_vec - vector of number of variables
% iters - # iterations to estimate variance components
% sigma_e - noise level 
% model_str - (optional) what is the model used 
%
% Output:
% mean_V - quasi-regression mean estimator for linear part of the variance
% std_V -  quasi-regression std estimator for linear part of the variance
% bias_V1 -  quasi-regression bias estimator  for linear part of the variance
% estimated_bias_V - 
% true_std_V1 - analytic computation of st.d. of quasi-regression estimator of variance
% true_std_V2 -  analytic computation of st.d. of quasi-regression estimator of variance of quadratic part
% 
function [mean_V std_V true_bias_V estimated_bias_V true_std_V1 true_std_V2 V_empirical_ave V_theoretical_ave] = ...
    simulate_quasi_regression(n_vec, d_vec, iters, sigma_e, model_str)

if(~exist('sigma_e', 'var') || isempty(sigma_e))
    sigma_e = 0;
end
if(~exist('model_str', 'var') || isempty(model_str))
    model_str = 'linear';
end
sigma_g = sqrt(1-sigma_e^2); 

num_n = length(n_vec); num_d = length(d_vec); max_n = max(n_vec); max_d = max(d_vec);
mean_V = zeros(num_n, num_d); std_V = zeros(num_n, num_d);
true_std_V1 = zeros(num_n, num_d);

n_start_vec = [1 n_vec(1:end-1)+1];
d_start_vec = [1 d_vec(1:end-1)+1];


switch model_str
    case 'linear'
        beta_vec = repmat(1 ./ sqrt(max_d), max_d, 1);  % This assumes no noise at all ... 
    case 'quadratic'
        beta_vec = repmat(1 ./ sqrt(max_d*(max_d-1)), max_d*(max_d-1), 1);
end
beta_vec = beta_vec .* sigma_g; 

quasi_reg_var_mat = zeros(iters, num_d, num_n);
quasi_reg_bias_mat = zeros(iters, num_d, num_n);


for k=1:num_n
    for j=1:num_d
        true_std_V1(k, j) = ( 2*d_vec(j)^2 + d_vec(j)* (10*n_vec(k)-12) + 2*(n_vec(k)-1)*(4*n_vec(k)-9) - ...
            8*(n_vec(k)-1)^2 / d_vec(j) ) / n_vec(k)^3; % this is the formula for variance for the linear case
        true_std_V2(k, j) = (26*d_vec(j)^4 - 27*d_vec(j)^3+(24*n_vec(k)-151)*d_vec(j)^2+(-3*n_vec(k)+156)*d_vec(j)+ ...
            2*(n_vec(k)-1)*(31*n_vec(k)-61)+ 8*(n_vec(k)-1)*(30*d_vec(j)^2 + 6*d_vec(j) - 12*n_vec(k)) * d_vec(j)*(d_vec(j) - 1)) / n_vec(k)^3; % this is the formula for variance for the quadratic case
    end
end
true_std_V1 = sqrt(true_std_V1); % take st.d.

for i=1:iters % loop on iterations
    run_iter = i
    X = 2*(rand(max_n, max_d) > 0.5) - 1; % Example: use Binary r.v.s. !!! +/- 1/2. Otherwise approximations don't work ..
    switch model_str
        case 'quadratic'
            %            X_prod = X * X'; % matrix holding x_j x_k for any j and k
    end
    y = zeros(max_n, num_d);
    for j=1:num_d % Generate observed variable (no noise currently)
        switch model_str
            case 'linear'
                y(:,j) = X(:,1:d_vec(j)) * repmat(1 ./ sqrt(d_vec(j)), d_vec(j), 1); % beta_vec(1:d_vec(j));
            case 'quadratic'
                beta = 1 / sqrt(d_vec(j)*(d_vec(j)-1));
                for k=1:d_vec(j) % can we somehow employ this together?
                    y(:,j) = y(:,j) + (X(:,1:d_vec(j)) .* repmat(X(:,k), 1, d_vec(j))) * repmat(beta, d_vec(j), 1);
                end
        end % switch model
        y(:,j) = y(:,j) + randn(max_n, 1) .* sigma_e; % Add noise 
    end
    
    
    
    for j=1:num_d % compute quasi-regression estimator. loop on dimension. Each time take only 1..d variables
        switch model_str
            case 'linear'
                tmp_quasi_reg_beta = zeros(max_d,1); tmp_quasi_reg_bias_vec = zeros(max_d,1);
                for k=1:num_n % loop on sample size
                    tmp_quasi_reg_bias_vec = tmp_quasi_reg_bias_vec + ...
                        (X(n_start_vec(k):n_vec(k),:).^2)' * y(n_start_vec(k):n_vec(k),j).^2;
                    tmp_quasi_reg_beta = tmp_quasi_reg_beta + ...
                        X(n_start_vec(k):n_vec(k),:)' * y(n_start_vec(k):n_vec(k),j);
                    quasi_reg_var_mat(i,j,k) = sum((tmp_quasi_reg_beta(1:d_vec(j)) ./ n_vec(k)).^2); % adjust beta as function of d to keep total variance constant
                    quasi_reg_bias_mat(i,j,k) = sum((tmp_quasi_reg_bias_vec(1:d_vec(j)) .* ((max_d / d_vec(j)) / n_vec(k)^2))); % adjust beta as function of d to keep total variance constant
                end
            case 'quadratic'
                tmp_quasi_reg_beta = zeros(max_d*max_d,1); tmp_quasi_reg_bias_vec = zeros(max_d*max_d,1);
                for k=1:num_n % loop on sample size
                    for j2=1:d_vec(j) % loop on one of the two x_j variables
                        tmp_quasi_reg_bias_vec(((j2-1)*d_vec(j)+1):j2*d_vec(j)) = tmp_quasi_reg_bias_vec(((j2-1)*d_vec(j)+1):j2*d_vec(j)) + ...
                            ((X(n_start_vec(k):n_vec(k),1:d_vec(j)).^2) .* repmat(X(n_start_vec(k):n_vec(k), j2).^2, 1, d_vec(j)))' * ...
                            y(n_start_vec(k):n_vec(k),j).^2;
                        tmp_quasi_reg_beta(((j2-1)*d_vec(j)+1):j2*d_vec(j)) = tmp_quasi_reg_beta(((j2-1)*d_vec(j)+1):j2*d_vec(j)) + ...
                            (X(n_start_vec(k):n_vec(k),1:d_vec(j)) .* repmat(X(n_start_vec(k):n_vec(k), j2), 1, d_vec(j)))' * ...
                            y(n_start_vec(k):n_vec(k),j);
                    end
                    quasi_reg_var_mat(i,j,k) = sum((tmp_quasi_reg_beta(1:d_vec(j)) ./ n_vec(k)).^2); % adjust beta as function of d to keep total variance constant
                    quasi_reg_bias_mat(i,j,k) = sum((tmp_quasi_reg_bias_vec(1:d_vec(j)) .* ((max_d / d_vec(j)) / n_vec(k)^2))); % adjust beta as function of d to keep total variance constant
                end
        end
    end
    
    %     for j=1:num_d % we don't need this
    %         [V_empirical{j}(i,:) V_theoretical{j}(i,:) ...
    %             V_diag_theoretical{j}(i) V_off_diag_theoretical{j}(i) ...
    %             V_diag{j}(i) V_off_diag{j}(i)] = compute_seven_parts(X(:,1:d_vec(j)), y(:,j));
    %     end
end


V = sum(beta_vec.^2);
mean_V = reshape(mean(quasi_reg_var_mat), num_d, num_n)';
std_V = reshape(std(quasi_reg_var_mat), num_d, num_n)';
true_bias_V = mean_V - V; % compute bias
estimated_bias_V = reshape(mean(quasi_reg_bias_mat), num_d, num_n)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j=1:num_d % we don't need this
%     V_empirical_ave{j} = mean(V_empirical{j})
%     V_theoretical_ave{j} = mean(V_theoretical{j})
%     should_be_zero = sum(V_theoretical_ave{j}) - true_std_V1(end,j)
%     V_diag_theoretical_ave{j} = mean(V_diag_theoretical{j})
%     V_diag_empirical_ave{j} = mean(V_diag{j})
%     V_off_diag_theoretical_ave{j} = mean(V_off_diag_theoretical{j})
%     V_off_diag_empirical_ave{j} = mean(V_off_diag{j})
%     should_be_all_ones =  V_empirical_ave{j} ./ V_theoretical_ave{j}
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% true_std_V1 = 2.*(n_vec -1).*(2.*n_vec-3) .* V1; % Compute true variance according to formula with seven terms
% true_std_V1 = true_std_V1 ./ (n_vec.^3); % normalize by n^3
% mu_2 =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Internal function for estimating the different variance componenets (goal
% is to compare empirical with theoretical)
function [V_empirical V_theoretical ...
    V_diag_theoretical V_off_diag_theoretical ...
    V_diag V_off_diag] = compute_seven_parts(X, y, model_str)

if(~exist('model_str', 'var') || isempty(model_str))
    model_str = 'linear';
end
[n d] = size(X);

V_diag_theoretical = d; % 2*(n-1)*d;
V_off_diag_theoretical = (d-1)*4/d; % 2*(n-1)*..

S2 = sum(X.^2, 2); % \sum_j x_j^2
mu2 = sum((X.^2)'*y.^2)/n; % sum_i x_i
X_times_y_sqr_mat = X .* repmat(y.^2, 1, d);


switch model_str
    case 'linear'
        V_theoretical(1) = 2*(n-1)* (d + 4 - 4/d); %   (3*d-2);
        V_theoretical(2) = d*(3*d-2); %    (3*d^3-11*d^2+26*d-16)/(2*d);
        V_theoretical(3) = 4*(n-1)*(3*d-2);
        V_theoretical(4) = 4*(n-1)*(n-2)*(3*d-2)/d;
        V_theoretical(5) = -d^2;
        V_theoretical(6) = -4*(n-1)*d;
        V_theoretical(7) =- 2*(n-1)*(2*n-3);
        V_theoretical = V_theoretical ./ (n^3);
        
        
        V_empirical(1) = 0;
        V_empirical_mat = (X_times_y_sqr_mat' * X) .^ 2 ./ n^2;
        V_empirical(1) = 2*(n-1)*sum(V_empirical_mat(:));
        V_empirical(1) = V_empirical(1) * n / (n-1) - 2*d *(3*d-2); %  / (n-1); % correct for bias !
        V_diag = trace(V_empirical_mat); V_off_diag = sum(V_empirical_mat(:))-trace(V_empirical_mat);
        V_empirical(2) = sum(y.^4 .* S2.^2)/n;
        V_empirical(3) = 4*(n-1)*sum(S2 .* y.^4)/n;
        V_empirical(4) = 4*(n-1)*(n-2) * sum(y.^4)/n;
        V_empirical(5) = -mu2^2;
        V_empirical(6) = -4*(n-1)*var(y) * mu2;
        V_empirical(7) = - 2*(n-1)*(2*n-3) * var(y);
        V_empirical = V_empirical ./ (n^3);
        
    case 'quadratic'
        V_theoretical(1) = -9999999999999; % ??? need to fill! %   (3*d-2);
        V_theoretical(2) = -9999999999999; %    (3*d^3-11*d^2+26*d-16)/(2*d);
        V_theoretical(3) = -9999999999999;
        V_theoretical(4) = -9999999999999;
        V_theoretical(5) = -9999999999999;
        V_theoretical(6) = -9999999999999;
        V_theoretical(7) = -9999999999999;
        V_theoretical = V_theoretical ./ (n^3);
        
        V_empirical(1) = -9999999999999;
        V_empirical(2) = -9999999999999;
        V_empirical(3) = -9999999999999;
        V_empirical(4) = -9999999999999;
        V_empirical(5) = -9999999999999;
        V_empirical(6) = -9999999999999;
        V_empirical(7) = -9999999999999;
        V_empirical = V_empirical ./ (n^3);
        
        
        
        
end

