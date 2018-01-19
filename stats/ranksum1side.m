% Wilcoxon rank sum test but with one-sided p-value
% Note: this function is not needed anymore due to the 'tail' option in Matlab's ranksum function
%
function [p, h, stats] = ranksum1side(x,y,varargin)
%RANKSUM Wilcoxon rank sum test for equal medians.
%   P = RANKSUM(X,Y) performs a one-sided rank sum test of the hypothesis
%   that two independent samples, in the vectors X and Y, come from
%   distributions with equal medians, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("medians are equal") is
%   true.  Small values of P cast doubt on the validity of the null
%   hypothesis.  The two sets of data are assumed to come from continuous
%   distributions that are identical except possibly for a location shift,
%   but are otherwise arbitrary.  X and Y can be different lengths.
%   The two-sided p-value is computed by doubling the most significant
%   one-sided value.
%
%   The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%
%   [P,H] = RANKSUM(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("medians are equal") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = RANKSUM(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = RANKSUM(...,'method',M) computes the p-value exactly if M is
%   'exact', or uses a normal approximation if M is 'approximate'.  If you
%   omit this argument, RANKSUM uses the exact method for small samples and
%   the approximate method for larger samples.
%
%   [P,H,STATS] = RANKSUM(...) returns STATS, a structure with one or two
%   fields.  The field 'ranksum' contains the value of the rank sum
%   statistic.  For the 'approximate' method, the field 'zval' contains the
%   value of the normal (Z) statistic.
%
%   See also SIGNTEST, SIGNRANK, KRUSKALWALLIS, TTEST2.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.14.4.5 $

[p, h, stats] = ranksum(x,y);
if(stats.zval < 0) % this means x < y
    p = 1-p/2;
else
    p = p/2; % move from two sided to one sided test
end

