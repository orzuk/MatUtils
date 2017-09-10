% Compute the difference between two distributions.
% The distributions do not have to be normalized (to sum to one)
% and the x values (bins) do not have to be equi-distanct. 
% Thus, they can represent both histograms (in which case they sum to one)
% or densities (in which case \sum(bin_size * p) = 1
%
% Input:
% x1 - the x vector of first distribution
% p1 - probabilities of first distribution
% x2 - the x vector of second distribution
% p2 - probabilities of second distribution
% diff_flag - how do we measure the difference. We currently implement three ways:
%                                   1 - FDR: Compute only the fraction detected at a certain FDR level q
%                                   2 - Conservative: take max x ...
%                                   3 - Liberal: take min x ....
%                                   4 - Intermidiate: take  mixture-of-gaussians
%                                   5 - 'minimal difference' - assume that one of the dists is at least epsilon away from the other
% diff_direction - which side to look at the difference (meaningful only for some of the methods)
% diff_param - additional method-specific parameter (e.g. the q level for fdr)
% shift - value we move ditributions with respect to one another to be conservative (default is 1)
% do_smooth (optional) - perform Gaussian smoothing on curves
% plot_flag - if to plot the two (scaled) distributions, along with the difference
% middle_val - value representing distribution mean/median/whatever distinguishes right/left
%
% Output:
% diff_frac - maximum difference fraction (dist. 1 has large area than dist. 2)
% diff_cut - where is the best cut seperating curves
% diff_x - x vector for the difference distribution 
% diff_dist - the difference distribution
% diff_frac_one_sided - NEW! sometimes the difference is two sided, in
% which case we also give the difference of each side (left%right)
%
function [diff_frac diff_cut diff_x diff_dist diff_frac_one_sided] = ...
    distributions_diff(x1, p1, x2, p2, ...
    diff_flag, diff_direction, diff_param, ...
    shift, do_smooth, plot_flag, legends_vec, middle_val, varargin)

AssignStatsConstants(); % set statistical constants
delta = 0.01; % don't try to estimate stuff below this threshold 
LEFT = 0; RIGHT = 1;

diff_frac = 0; diff_cut =0; % diff_dist = p1; % set defaults

if(any(isnan(p2))) % deal with a case where the 2nd distribution is not 'proper'
    diff_frac = 1.0;  diff_dist = p1;
    if(diff_direction == RIGHT)
        diff_cut = min(x1);
    else
        diff_cut = max(x1);
    end
    diff_x = x1; diff_dist = p1;
    return;
end

if(~exist('diff_flag', 'var') || isempty(diff_flag)) % default is conservative
    diff_flag = 1;
end
if(~exist('shift', 'var') || isempty(shift))
    shift=1;
end
if(~exist('do_smooth', 'var') || isempty(do_smooth))
    do_smooth = 0;
end
if(~exist('middle_val', 'var') || isempty(middle_val))
    middle_val = mean_hist(x2, p2);
end


sum_p1 = sum(p1); sum_p2 = sum(p2);
p1 = p1 ./ sum(p1); p2 = p2 ./ sum(p2); % normalize distributions (they could come from density functions)


if(do_smooth)% perform Gaussian smoothing BEFORE doing interpolation 
    chop = 1;
    if(do_smooth == 1) % default sigma
        sigma = GAUSS_SMOOTH_SIGMA;
    else
        sigma = do_smooth;
    end
    p1 = gauss_smooth(p1, sigma, chop);
    p2 = gauss_smooth(p2, sigma, chop);
end    

old_working = 0; % old method: refinement by union. New: try uniform refinement!
if(old_working)
    diff_x = union(x1, x2);
else
    diff_x = union(x1, x2); diff_min = min(diff_x); diff_max = max(diff_x); scale_factor=2;
    new_n = round(length(diff_x) * scale_factor);
    new_step = (diff_max - diff_min) / (new_n-1); 
    diff_x = diff_min:new_step:diff_max;
end
n = length(diff_x); % First generate a refinement of the distributions and bring them both to the same x vec
p1 = interp1(x1, p1, diff_x, 'linear', 0); p1 = normalize_hist(diff_x, p1);
p2 = interp1(x2, p2, diff_x, 'linear', 0); p2 = normalize_hist(diff_x, p2); % p2 ./ sum(p2);

tmp_diff_direction = 0; 
if( (diff_direction == LEFT) && (diff_flag == FDR) )
    p1 = p1(end:-1:1);
    p2 = p2(end:-1:1);
    tmp_diff_direction = 1;
    diff_direction = RIGHT;
end

% p_tmp = zeros(n,1); [inter_x I1 J1] = intersect(diff_x, x1); p_tmp(I1) = p1(J1); p1 = p_tmp; 
% p_tmp = zeros(n,1); [inter_x I2 J2] = intersect(diff_x, x2); p_tmp(I2) = p2(J2); p2 = p_tmp;

cum_p1 = cumsum(p1); cum_p2 = cumsum(p2); Z = cum_p1(end); % Z is normalization constant

switch diff_flag
    case FDR % use FDR (here enable many different cutoffs)
        fdr_q = diff_param; n = length(fdr_q); % just change name for better readability
        if(diff_direction == RIGHT)
            diff_vec = (Z - cum_p2) ./ ( (Z-cum_p2) + (Z-cum_p1) );
        else
            diff_vec = cum_p2 ./ (cum_p2 + cum_p1);
        end
        diff_frac = zeros(n,1); diff_cut = zeros(n,1); 
        for i=1:n
            tmp_diff = find(diff_vec < fdr_q(i), 1);
            if(isempty(tmp_diff))
                diff_frac(i) = 0;
            else
                diff_frac(i) = tmp_diff;
            end
            if(diff_frac(i) > 0)
                diff_cut(i) = diff_x(diff_frac(i));
                diff_frac(i) = Z - cum_p1(diff_frac(i));
            else
                if(diff_direction == RIGHT)
                    diff_cut(i) = max(diff_x);
                else
                    diff_cut(i) = min(diff_x);
                end
            end
        end
        diff_frac = diff_frac ./ Z; % normalize to get fraction 
    case FDR_INV % use FDR, but with input as cutoff (not FDR level)
        diff_cut = diff_param; % just change name for better readability
        diff_frac = find(diff_x >= diff_cut, 1);
        diff_frac = 1 - cum_p1(diff_frac)/Z;
    case APARAM_INTERSECT % methods for separating curves  % use conservative estimate: (1-F1) <= (1-F2) right most point
        cum_p1 = cum_p1 ./ cum_p1(end); cum_p2 = cum_p2 ./ cum_p2(end); % normalize again
        if(diff_direction == RIGHT) % right-tail difference between two distributions (this is the Kolmogorov-Smirnov test statistics)
            f_start = max( find(cum_p1 > 0.5, 1), find(cum_p2 > 0.5, 1) ); % get a reasonable cutoff
            [diff_frac diff_cut] = max(cum_p2(f_start:max(f_start,end-shift)) - ...
                cum_p1(min(f_start+shift,n):end)); diff_cut = diff_cut + f_start - 1; diff_cut = diff_x(diff_cut);
        else % left-tail difference between two distributions (this is the Kolmogorov-Smirnov test statistics)
            f_start = min( find(cum_p1 < 0.5, 1, 'last'), find(cum_p2 < 0.5, 1, 'last') ); % get a reasonable cutoff
            if(isempty(f_start)) % deal with case where first bin is already very large
                f_start=1;
            end
            [diff_frac diff_cut]= max(cum_p1(1:max(1,f_start-shift)) - ...
                cum_p2(min(1+shift,f_start):f_start)); diff_cut = diff_x(diff_cut);
        end
    case  APARAM_DIFF % use liberal estimate:min_x (p2(diff_x) / p1(diff_x)). We need to perform Gaussian smoothing?
        min_alpha = 0.02; 
        p1_for_diff = p1; cumsum_p1_diff = cumsum(p1_for_diff) ./ sum(p1_for_diff);
        min_p1_ind = find(cumsum_p1_diff > min_alpha, 1);
        max_p1_ind = find(cumsum_p1_diff < 1-min_alpha, 1, 'last');
        p2_for_diff = p2; cumsum_p2_diff = cumsum(p2_for_diff) ./ sum(p2_for_diff);
        min_p2_ind = find(cumsum_p2_diff > min_alpha, 1);
        max_p2_ind = find(cumsum_p2_diff < 1-min_alpha, 1, 'last');
        min_p_ind = max(min_p1_ind, min_p2_ind);
        max_p_ind = min(max_p1_ind, max_p2_ind);
        
%%%        save_p2 = p2; p2(767:end) = p1(767:end); p2= normalize_hist(diff_x, p2); 
        diff_frac = 1 - min(max(p1(min_p_ind:max_p_ind), delta) ./ ...
            max(p2(min_p_ind:max_p_ind), delta));  % Get the upper-bound
        [dummy diff_cut] = min(abs(cum_p1 - cum_p2 - diff_frac)); diff_cut = diff_x(diff_cut);

    case APARAM_CUM_DIFF % similar to APARAM_DIFF but on the cumulative distributions (from Siepel's paper)
        cumsum_p1 = cumsum(p1); cumsum_p2 = cumsum(p2); 
        diff_frac = 1 - min(max(cumsum_p1(min_p_ind:max_p_ind), delta) ./ ...
            max(cumsum_p1(min_p_ind:max_p_ind), delta));  % Get the upper-bound
        [dummy diff_cut] = min(abs(cum_p1 - cum_p2 - diff_frac)); diff_cut = diff_x(diff_cut);
        
        
    case  MIX_GAUSS % use mixture-of-gaussians model, where we know one Gaussian (not supported yet)
        n=100000; % number of data points
        data = distrnd(diff_x, p1, n, 1); iters=100;
        prior = [0.5 0.5]; mu = [mean_hist(diff_x,vec2row(p2)), mean_hist(diff_x,vec2row(p1))];
        sigma = [std_hist(diff_x,vec2row(p2)), std_hist(diff_x,vec2row(p1))];
        [P_F M_F S_F LogLike_F] = MixtureOfGaussiansGivenInit(data,2,iters, prior, mu, sigma, ...
            [0 0], [1 0], [1 0]); % we know and freeze mu_1 and sigma_1
        diff_frac = P_F(2);
    case UPPER_BOUND_EPS  % use minimal difference epsilon . The formula is: p0 = (mu_0 - mu + eps) / mu_0
        eps = diff_param;
        mu1 = mean_hist(diff_x, vec2row(p1)); mu2 = mean_hist(diff_x, vec2row(p2));
        if(diff_direction == RIGHT)
            diff_frac = 1 - eps / (mu1 - mu2 + eps);
        else
            diff_frac = 1 - eps / (mu2 - mu1 + eps);
        end
        diff_frac = max(diff_frac, 0); % might be negative if means are in reverse order !!! 
end % switch

diff_dist = ((p1 - (1-diff_frac(1)) .* p2) ./ diff_frac(1)) .* sum_p1; % Compute difference distribution
if(~exist('diff_cut', 'var'))  % here take the intersection point where which dist. is maximal
    if(diff_direction == RIGHT)
        [dummy diff_cut] = min(  diff_dist >= p1 ); % get the index
    else
        [dummy diff_cut] = max(  diff_dist >= p1 ); % get the index
    end
    diff_cut = diff_x(diff_cut); % move from index to value
end

if(tmp_diff_direction)
    diff_cut = max(diff_x) - (diff_cut - min(diff_x));
end
    
if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 0;
end
if(plot_flag)% Plot two distributions and their difference    
    figure; subplot(2,2,1); hold on; 
    plot(diff_x, p1, '.'); plot(diff_x, p2, 'r.'); title('two dists. unscaled');
    legend('dist. 1', 'dist. 2');
    subplot(2,2,2); hold on; 
    plot(diff_x, max(p1, delta) ./ max(p2, delta), '.'); title('two dists. ratio p_1 / p_2');
    legend('p_1 / p_2'); line( [min(diff_x) max(diff_x)], [1-diff_frac(1) 1-diff_frac(1)], 'color', 'red'); 
    subplot(2,2,3); hold on;
    plot(diff_x, max(p2, delta) ./ max(p1, delta), '.'); title('two dists. ratio p_2 / p_1');
    legend('p_2 / p_1');
    
    if(~exist('legends_vec', 'var'))
        legends_vec = {'dist. 1', 'dist. 2', 'diff'};
    end
    subplot(2,2,4); hold on; title(['two dists. and their diff. (scaled) ' num2str(100*diff_frac(1)) ' % ']);
    plot(diff_x, p1, '.'); plot(diff_x, (1-diff_frac(1)) .* p2, 'r.');
    plot(diff_x, diff_frac(1) .* normalize_hist(diff_x, diff_dist), 'g.');    
    legend(legends_vec);
end

left_inds = find(diff_x<= middle_val); 
right_inds = find(diff_x > middle_val);
if(nargout > 4) % assume just one output
    diff_frac_one_sided(1) = diff_frac(1) * sum(diff_dist(left_inds)) / sum(diff_dist); % right side
    diff_frac_one_sided(2) = diff_frac(1) * sum(diff_dist(right_inds)) / sum(diff_dist); % right side
end

