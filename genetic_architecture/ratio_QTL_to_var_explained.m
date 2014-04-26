% Compute variance explaiend for a rare allele when sequencing to extremes of a distribution
% Find the most likely values of beta and f (using maximal likelihood).
% We need to compute this numerivally since we've got a sum of two gaussians
% (no analytic solution)
%
% Input:
% n1 - # of individuals sequenced in left tail
% n2 - # of individuals sequenced in right tail
% r1 - # of occurences of rare allele in individuals in the left tail
% r2 - # of occurences of rare allele in individuals in the right tail
% t1 - fraction of population in left tail
% t2 - one minus fraction of population in right tail
% gamma - fraction of functional rare alleles
%
% Output:
% beta - inferred effect size
% f - inferred aggregate allele frequency of all rare alleles
% V - inferred variance explained
% min_gamma - lower bound on gamma, the fraction of functional rare alleles
%
function  [beta f V min_gamma] = ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, ...
    gamma, input_f) % Use shamil's formula

if(~exist('input_f', 'var') || isempty(input_f))
    input_f = 0.01; % temp. Fix the allele frequency f (this is aggregate of ALL rare variants, not just functional)
    search_f_flag = 1;
else
    search_f_flag = 0;
end

AssignGeneralConstants;
fancy_figures_flag = 0; % do many checks and figures
R = (r2/n2)/(r1/n1); % compute enrichment of right vs. left tail

if(nargout > 3) % output also minimal gamma
    min_gamma = ratio_to_gamma_lowerbound_internal(R, [t1 t2], input_f); % [s1 s2],
%    min_gamma_bimodal = ratio_to_gamma_lowerbound_internal(R, [t1 t2], input_f, 1); % [s1 s2],
    gamma = max(gamma, min_gamma + 0.001); % don't allow gamma values lower than possible
end

%x_vec = linspace(0, 1, 2000); x_vec = x_vec(2:end-1);
%beta_grid = linspace(-20, 20, 1000); % possible beta values


if(~search_f_flag) % here we get f as input. Perform 1d optimization
    f = input_f; min_beta = -10; max_beta = 10; % set beta possible range
    beta_MLE = fminbnd(@(x) -loglikelihood_QTL_internal( ...
        f, x, [n1 n2], [r1 r2], [t1 t2], gamma), min_beta, max_beta); % perform 1-d optimization using likelihood. n counts individuals
    beta2 = fminbnd(@(x) abs(ratio_QTL_internal(x, R, [t1 t2], gamma, f)), ...
        -10, 10); % a different approach: solve for beta  (not used) 
    both_beta_should_be_equal = beta2 - beta_MLE
    beta = beta_MLE % take the MLE estimator 
else % search the 2d space of f and beta
    beta = 0.1; % initilize beta for search
    f = ((r1/(2*n1)) + (r2/(2*n2)))/2; % Temp! take average of two frequencies. Divide by another two due to two chromosomes (diplod)
    f_beta = fminsearch(@(x) -loglikelihood_QTL_2d_input_internal( ...
        x, [n1 n2], [r1 r2], [t1 t2], gamma), [f beta]); % perform 2-d optimization
    f = f_beta(1); beta = f_beta(2);
end



R_substituted = beta_to_ratio_internal(beta, gamma, [t1 t2], f); % Check R computation in two ways (doesn't affect results) 
R_diff_should_be_zero = R -  R_substituted
if(abs(R_diff_should_be_zero) > 0.01)
    two_Rs_dont_agree_something_is_wrong = 234324
    should_be_zero = R_substituted - R;
end


% Compute variance explained - normalize beta by the trait's variance 
sigma_Z = sqrt(1 + beta^2*gamma*f*(1-gamma*f)); % mu_Z = gamma*f*beta; % what is this 
beta = beta ./ sigma_Z; % normalize beta to st.d. units
V = beta_to_heritability(beta, gamma*f, 1, 'diploid'); % compute var explained 

if(fancy_figures_flag) % what figures do we show here? need to include parameters 
    draw_rare_variants_figures_internal()
end % draw figures (currently not used) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old Stuff
%
% for i=1:2
%     max_f2(i) = fminbnd(@(x) -loglikelihood_QTL_internal(x, max_beta(i), [n1 n2], [r1 r2], [t1 t2], gamma), ...
%         0, 1)/2;
% end
%
% ll_vec = loglikelihood_QTL_internal(x_vec, beta, [n1 n2], [r1 r2], [t1 t2], gamma);
% figure; plot(x_vec, ll_vec, '.'); title('log-likelihood'); xlabel('RAF'); ylabel('LL');
%
% [max_val max_ind] = max(ll_vec); max_xxx = x_vec(max_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Internal function which needs to be zero. x represents beta
function ret = ratio_QTL_internal(x, R, t, gamma, f)

ret = beta_to_ratio_internal(x, gamma, t, f) - R;

%ret = gamma*(normcdf(t2-x) + R*normcdf(t1-x)) + ... % this needs to be solved for x (beta)
%    (1-gamma)*(normcdf(t2) - 1 + R*normcdf(t1)) - gamma; % moved terms around in formula from Shamil


% Compute log-likelihood when beta and f are given seperately
%
% Input:
% f_beta - vector of size 2 [f beta]
% n - number of individuals
% r - number of rare alleles
% t - quantiles for two tails
% gamma - fraction of functional rare alleles
%
% Output:
% LL - log likelihood
%
function [LL LL1 LL2] = loglikelihood_QTL_2d_input_internal(f_beta, n, r, t, gamma)

f = f_beta(1:end-1); beta = f_beta(end);
[LL LL1 LL2] = loglikelihood_QTL_internal(f, beta, n, r, t, gamma);

% Internal function for computing the log-likelihood of allele frequency given two tails
%
% Input:
% f - allele frequency
% beta - effect size
% n - number of individuals
% r - number of rare alleles
% t - quantiles for two tails
% gamma - fraction of functional rare alleles
%
% Output:
% LL - log likelihood
% LL1 - first part of likelihood
% LL2 - second part of likelihood 
%
function [LL LL1 LL2] = loglikelihood_QTL_internal(f, beta, n, r, t, gamma)

s(1,:) = quantile_to_threshold_local(t(1), gamma, f, beta); % Compute thresholds. Currently slow
s(2,:) = quantile_to_threshold_local(t(2), gamma, f, beta); 

f_beta_t1_t2_alpha = [f beta t(1) t(2) gamma]
s_is = s 


% Compute probabilities of having a rare allele
p = (f./t(1)).*( gamma.*normcdf(s(1,:)-beta)+(1-gamma).*normcdf(s(1,:)) ); % Pr(x=1 | z < s_1)
p = min(p, 0.99999999); % avoid rounding/underflow errors
LL1 = r(1) .* log(p) + (2*n(1)-r(1)).*log(1-p); % 2n is number of chromosomes
p = (f./(1-t(2))).*( gamma.*(1-normcdf(s(2,:)-beta))+(1-gamma).*(1-normcdf(s(2,:))) ); % Pr(x=1 | z > s_2)
p = min(p, 0.99999999); % avoid rounding/underflow errors
LL2 = r(2) .* log(p) + (2*n(2)-r(2)).*log(1-p); % 2n is number of chromosomes
LL = LL1 + LL2;


% simple model from Shamil
%    p = (f2*normcdf(t(i)-beta) + (f-f2)*normcdf(t(i))) ./ ...
%        (f2*normcdf(t(i)-beta) + (1-f2)*normcdf(t(i)));
% simple model from Shamil
%    p = (f2*(1-normcdf(t(i)-beta)) + (f-f2)*(1-normcdf(t(i)))) ./ ...
%        (f2*(1-normcdf(t(i)-beta)) + (1-f2)*(1-normcdf(t(i))));
%    LL = LL + r(i) * log(f2*normcdf(t(i)-beta)) + (n(i)-r(i))*log((1-f2)*normcdf(t(i))) - ...
%        n(i)*log(f2*normcdf(t(i)-beta)+(1-f2)*normcdf(t(i)));
%    LL = LL + r(i) * log(f2*(1-normcdf(t(i)-beta))) + (n(i)-r(i))*log((1-f2)*(1-normcdf(t(i)))) - ...
%        n(i)*log(f2*(1-normcdf(t(i)-beta))+(1-f2)*(1-normcdf(t(i))));


% Formula for R. Modified from Shamil
% Now R depends also on f (not just gamma and beta. We cannot decouple
% these ...)
% Input:
% beta - effect size (additive)
% gamma - fraction of functional rare variants
% t - quantiles
% f - allele frequency of rare alleles
%
% Output:
% R - ratio of rare alleles in right to left tail
% R2 - same
% s - thresholds corresponding to quantiles
%
function [R R2 s] = beta_to_ratio_internal(beta, gamma, t, f)

s = zeros(2, length(beta));
for i=1:2
    s(i,:) = quantile_to_threshold_local(t(i), gamma, f, beta);
end
%mu_Z = gamma*f*beta;
%sigma_Z = sqrt(1 + beta^2*gamma*f*(1-gamma*f));
%s = norminv(t) * sigma_Z + mu_Z;
% s(2) = norminv(1-t(2)) * sigma_Z + mu_Z;

R = (gamma.*(1-normcdf(s(2,:)-beta)) + (1-gamma).*(1-normcdf(s(2,:)))) .* t(1) ./ ...
    ((gamma.*normcdf(s(1,:)-beta) + (1-gamma).*normcdf(s(1,:))) .* (1-t(2))); % scale by two tails areas !!!

R2 = -t(1) .* ((f-1).*normcdf(s(2,:)) + t(2) - f) ./ ...
    ((1-t(2)) .* ((f-1).*normcdf(s(1,:)) + t(1))); % Should give the same value

R_err = abs(R2-R)./max(1,R); % compute relative error
if(max(R_err) > 0.00000001)
    error('BAD R COMPUTATION!!!!');
end


% Simplified formula (from Shamil. Doesn't depend on f)
%R = (gamma.*(1-normcdf(t(2)-beta)) + (1-gamma).*(1-normcdf(t(2)))) .* normcdf(t(1)) ./ ...
%    ((gamma.*normcdf(t(1)-beta) + (1-gamma).*normcdf(t(1))) .* (1-normcdf(t(2)))); % scale by two tails areas !!!




% Compute the minimal fraction of rare alleles given an enrichment ratio in Shamil's model
%
% Input:
% R - ratio of rare alleles in right vs left tails
% t - two quantile thresholds
% f - allele frequency of rare alleles
% bimodal_flag - here we require that the phenotypic distribution is not bimodal
%
% Output:
% gamma_min - minimum possible value of gamma
%
function gamma_min = ratio_to_gamma_lowerbound_internal(R, t, f, bimodal_flag)

if(~exist('bimodal_flag', 'var') || isempty(bimodal_flag))
    bimodal_flag = 0;
end

% New: take case 3 when \beta -> infty and case 1 when \beta -> -infty
gamma_min = max( t(1).*(1-R) ./ (t(1).*(1-R) + R.*(1-f)), ...
    (1-t(2)).*(R-1) ./ ((1-t(2)).*R+t(2)-f) );

gamma_min = max(gamma_min, 0); % sometimes bound gets negative



if(bimodal_flag)
    gamma_max = 1; % (gamma_min was already set)
    is_bimodal_flag = 0;
    while(gamma_max-gamma_min > 0.00001)
        gamma_mid = (gamma_min + gamma_max)/2;
        beta = fminbnd(@(x) abs(ratio_QTL_internal(x, R, t, gamma_mid, f)), ...
            -10, 10); % solve for beta
        x_vec = min(-5,-5+beta):0.0001:max(5,5+beta); % joint x-vec for all distributions
        z_vec = normpdf(x_vec); % standard gaussian
        
        mix_vec = z_vec .* 0.5 + ... % (1-gamma_mid.*f) + ...
            normpdf(x_vec-beta) .* 0.5; % gamma_mid.*f;
        is_bimodal_flag = isbimodal_hist(x_vec, mix_vec, f*gamma_mid*0.1);
        %        is_bimodal_flag = (abs(beta) > 3); % crude way! allow maximum effect size of 3 st.d.s.
        if(is_bimodal_flag)
            gamma_min = gamma_mid;
        else
            gamma_max = gamma_mid;
        end
    end
    
    xxx = gamma_mid;
    %    [x p] = 0;% generate MoG distribution with the desired parameters
    
end


% New MoG model:
%gamma_min = max(1 - t(1) / (normcdf(s(1))*(t(1)+R*(1-t(2)))), ...
%    1 - R*(1-t(2)) / (R*(1-normcdf(s(1)))*(1-t(2))+ (1-normcdf(s(2)))*t(1))  );

% Old simple Shamil's model
%gamma_min = max( normcdf(t(1))*(1-R) / (normcdf(t(1))*(1-R) + R), ...
%    (1-normcdf(t(2))) * (R-1) / (  (1-normcdf(t(2)))*(R-1)+1 ) );

% Determine the maximum/minimum possible variance explained
function [V_min V_max] = ratio_to_var_explained_min_max_internal(R, t, f)

beta = fminbnd(@(x) abs(ratio_QTL_internal(x, R, [t(1) t(2)], 1, f)), ...
    -10, 10); % solve for beta, with gamma=1
tmp = beta.^2.*f.*(1-f);
V_min = 2.*tmp./(1+tmp);

tmp = (1-t(2)).*(R-1); % take beta -> infinity
V_max = 2 .* (tmp + 1-f).^2 ./ (tmp.*(tmp+1));

tmp = t(1).*(R-1); % take beta -> -infinity
V_max = max(V_max, 2.*(tmp + R.*(1-f)).^2 ./ (tmp.*(tmp+R)));

% compute MixtureOfGaussians threshold numerically
% Input:
% t - desired quantile
% gamma - fraction of functional rare alleles
% f - combined allele frequency of rare alleles
% beta - effect size of rare alleles
%
% Output:
% s - the threshold giving the quantile t
%
function s = quantile_to_threshold_local(t, gamma, f, beta)

num_opt = max(length(beta), length(f)); % enable multiple beta to optimize 
if(length(f) == 1)
    f = repmat(f, num_opt, 1);
end
if(length(beta) == 1)
    beta = repmat(beta, num_opt, 1);
end
s = zeros(num_opt, 1); % size(beta));
s(1) = norminv(t) + gamma*beta(1)*f(1);
for i=1:num_opt % this is the slow loop: how to make this more efficient???
    s(i) = fzero(@(x) gamma*f(i)*normcdf(x-beta(i)) + (1-gamma*f(i))*normcdf(x) - t, ...
        s(max(1,i-1))); % provide approximate starting point from the previous s
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figures
function   draw_rare_variants_figures_internal()


[y_vec y_vec1 y_vec2] = loglikelihood_QTL_internal( ...
    f, beta_grid, [n1 n2], [r1 r2], [t1 t2], gamma);
[MAX_LL MAX_I] = max(y_vec);
plot_inds = max(1,MAX_I-250):min(length(y_vec),MAX_I+100);

figure; hold on;
plot(beta_grid(plot_inds), y_vec(plot_inds), '.');
plot(beta_grid(plot_inds), y_vec1(plot_inds), 'c.');
plot(beta_grid(plot_inds), y_vec2(plot_inds), 'g.'); %plot(beta_grid, y_vec1+y_vec2, 'm.');
plot(beta_grid(MAX_I), y_vec(MAX_I), 'rx', 'linewidth', 3);
title(['Log-Likelihood for different effect sizes, \gamma=' num2str(gamma,2)]);
xlabel('\beta'); ylabel('LL');
legend('LL', 'LL-left', 'LL-right', 'MLE');



gamma_vec = linspace(min_gamma+0.000000001, 1, 100); % start with small vector
for i=1:length(gamma_vec)
    run_gamma_i = i
    beta_vec(i) = fminbnd(@(x) abs(ratio_QTL_internal(x, R, [t1 t2], gamma_vec(i), input_f)), ...
        -10, 10); % solve for beta
    f_vec(i) = fminbnd(@(x) -loglikelihood_QTL_internal(x, beta_vec(i), ...
        [n1 n2], [r1 r2], [t1 t2], gamma), 0, 1)/2;
end
%    mu_Z = gamma.*f.*beta_vec;
sigma_Z = sqrt(1 + beta_vec.^2.*gamma.*f.*(1-gamma.*f));
beta_vec = beta_vec ./ sigma_Z; % normalize beta to st.d. units
%end
%if(fancy_figures_flag)

[V_min V_max] = ratio_to_var_explained_min_max_internal(R, [t1 t2], f); % get minimum and maximum values
V_explained_vec = beta_to_heritability(beta_vec, gamma_vec.*f_vec, 1, 'diploid'); % compute var explained (what the hell is the allele frequency?)
figure; hold on; % plot ...
subplot(2,1,1);
plot(gamma_vec, V_explained_vec, '.'); line([min_gamma min_gamma], [0 max(V_explained_vec)], 'color', 'r');
%    line([min_gamma, 1], [V_min V_min], 'color', 'g');
%    line([min_gamma, 1], [V_max V_max], 'color', 'g');
title('Var. explained as a function of \gamma');
xlabel('\gamma'); ylabel('var. explained'); % New: fraction of functional variants is gamma
ylim([0 0.1]); % why cut this

subplot(2,1,2);
plot(gamma_vec, abs(beta_vec), '.');line([min_gamma min_gamma], [0 max(abs(beta_vec))], 'color', 'r');
title('Effect size (\beta) as a function of \gamma');
xlabel('\gamma'); ylabel('\beta');
ylim([0 4]);

%    my_saveas(gcf, '../../common_disease_model/docs/figs/QTL_two_tails_heritability_vs_gamma', format_fig_vec);

subplot(2,2,3);
plot(gamma_vec, f_vec, '.'); line([min_gamma min_gamma], [0 max(f_vec)], 'color', 'r');
title('Risk Allele Freq. (f) as a function of \gamma');
xlabel('\gamma'); ylabel('RAF');

subplot(2,2,4);
plot_x_vec = linspace(-4,4,1000);
plot_y_vec = normpdf(plot_x_vec);
plot(plot_x_vec, plot_y_vec, '.');
line([s1 s1], [0 0.4], 'color', 'k');
line([s2 s2], [0 0.4], 'color', 'k');
xlabel('Z'); ylabel('freq.');
title(['r_1=' num2str(r1) ', n_1=' num2str(n1) ', r_2=' num2str(r2) ', n_2=' num2str(n2) ...
    ', R=' num2str(R, 3) ', t_1=' num2str(100*t1,3) '%, t_2=' num2str(100*t2,3) '%']);
%
% beta = fminbnd(@(x) abs(ratio_QTL_internal(x, R, [t1 t2], gamma)), -10, 10); % solve for beta
%x_vec =  linspace(-100, 10, 2000);
%figure; plot(x_vec, abs(ratio_QTL_internal(x_vec, R, [t1 t2], gamma)), '.');

%f = 0.5; % Temp!!! should take some weighted average of two allele frequencies at extremes
% min_x = min(r1*0.5/n1, r2*0.5/n2);
% max_x = max(r1*0.5/n1, r2*0.5/n2);
% f = fminbnd(@(x) -loglikelihood_QTL_internal(x, beta, [n1 n2], [r1 r2], [t1 t2], gamma), ...
%     0*min_x, 1+0*max_x)/2;

if(1 < 0)
    %       if(fancy_figures_flag) % Plot log-likelihood 2-D figure (heavy part !!)
    for i=1:length(beta_grid)
        run_i = i
        LL_mat(i,:) = loglikelihood_QTL_internal(x_vec, beta_grid(i), ...
            [n1 n2], [r1 r2], [t1 t2], gamma); % This is the heavy loop !!!!
    end
    figure; imagesc_with_labels(real(LL_mat), x_vec, beta_grid, [], [], 100); % , skip, varargin)
    %imagesc(real(LL_mat)); colorbar;
    [A B] = max(real(LL_mat(:))); [I J] = ind2sub(size(LL_mat), B);
    %max_beta(1) = beta; max_f(1) = f;
    max_beta = beta_grid(I), max_f = x_vec(J)/2;
    %       end
end



% Simulate data in two tails
function simulate_rare_variants_tails(beta, f, gamma, t1, t2, n1, n2)

%beta = 1;
% figure; hold on;
% [data X] = MixtureOfGaussiansSimulateData([1-gamma*f gamma*f], [0 beta], [1 1], 1000); X=X-1;
% MixtureofGaussiansDraw1dGaussians(data', [1-gamma*f gamma*f], [0 beta], [1 1], {'null', 'shifted'}, ...
%      {'null', 'shifted', 'sum'}, 'gr:', [], 100, [], 0);
% t_draw(1) = t1 + beta*gamma*f;
% t_draw(2) = t2 + beta*gamma*f;
%
% f_empirical(1) = mean(X(find(data<t_draw(1))))
% f_empirical(2) = mean(X(find(data>t_draw(2))))
% f_observed(1) = r1/(2*n1)
% f_observed(2) = r2/(2*n2)
% for i=1:2
%    line([t_draw(i) t_draw(i)], [0 0.2]);
% end
%  xx = 99
