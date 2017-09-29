% Builds a histogram by binning the elements of val into containers, one
% for each unique value of val.  Each element of val is weighted by the
% value in the corresponding column of weights.
%
% Input:
% vals = a 1xn matrix of positive integer values.  The (weighted) frequency of these values
% will be made into the output histogram. (why must they be integers?)
% weights  -  a 1xn matrix of weights (default is ones(1,n))
% bins - number of bins, or the bins themselves
% min_val = the minimum integer bin to be included on the output histogram
% (default is 0)
%
% Output:
% H = a row vector of weighted counts, such that H(x) is the weighted
% number of times min_val+(x-1) appears in vals
%
function [H, bins] =  weighted_hist(vals, weights, bins, min_val)

if ~exist('weights','var')
    weights = ones(1, size(vals,2));
end
if exist('min_val','var')
    vals = vals - min_val+1;
end

epsilon = 0.00000000000001;
if(length(bins) == 1) % input is number of bins
    [~, bins] = hist(vals, bins);
    bins = [min(vals)-epsilon bins max(vals)+epsilon]; % make sure all are represented 
end
[~, bin_inds] = histc(vals,bins);
missing_bins = setdiff((1:length(bins)-1), unique(bin_inds));
if(~isempty(missing_bins)) % here there are missing bins 
    bin_inds = [vec2row(bin_inds) missing_bins]; 
    weights = [vec2row(weights) zeros(1, length(missing_bins))]; 
end
H = (accumarray(vec2column(bin_inds), vec2column(weights)))';
if(length(bins) > length(H))
    bins = 0.5 .* (bins(1:end-1) + bins(2:end));
end

