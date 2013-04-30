function imputed = knnimpute(data,K,varargin)
%KNNIMPUTE imputes missing data using the nearest-neighbor method
%
%   KNNIMPUTE(DATA) replaces NaNs in DATA with the corresponding value from
%   the nearest-neighbor column using Euclidean distance. If the nearest-
%   neighbor column also contains a NaN value, then the next nearest column
%   is used.
%
%   KNNIMPUTE(DATA, K) replaces NaNs in DATA with a weighted mean of the K
%   nearest-neighbor columns. The weights are inversely proportional to the
%   distances from the neighboring columns.
%
%   KNNIMPUTE(..., 'DISTANCE',DISTFUN) computes nearest neighbors using
%   distance metric DISTFUN.  Choices are
%
%       'euclidean'   Euclidean distance
%       'seuclidean'  Standardized Euclidean distance, in which each
%                     coordinate in the sum of squares is inverse weighted
%                     by the sample variance of that coordinate
%       'cityblock'   City block distance
%       'mahalanobis' Mahalanobis distance
%       'minkowski'   Minkowski distance with exponent 2
%       'cosine'      One minus the cosine of the included angle
%                     between observations (treated as vectors)
%       'correlation' One minus the sample correlation between
%                     observations (treated as sequences of values)
%       'hamming'     Hamming distance, percentage of coordinates
%                     that differ
%       'jaccard'     One minus the Jaccard coefficient, the
%                     percentage of nonzero coordinates that differ
%       'chebychev'   Chebychev distance (maximum coordinate difference)
%       function      A distance function specified using @, for
%                     example @DISTFUN
%
%   See PDIST for more details.
%
%   KNNIMPUTE(..., 'DISTARGS',ARGS) passes the arguments ARGS to the
%   function DISTFUN. ARGS can be a single value or a cell array of
%   values.
%
%   KNNIMPUTE(...,'WEIGHTS',W) allows you to specify the weights used in
%   the weighted mean calculation. W should be a vector of length(K).
%
%   KNNIMPUTE(...,'MEDIAN',true) uses the median of the K nearest neighbors
%   instead of the weighted mean.
%
%   Examples:
%
%       load yeastdata
%       % Remove data for empty spots
%       emptySpots = strcmp('EMPTY',genes);
%       yeastvalues(emptySpots,:) = [];
%       genes(emptySpots) = [];
%       % Impute missing values
%       imputedValues = knnimpute(yeastvalues);
%
%   See also ISNAN, KNNCLASSIFY, NANMEAN, NANMEDIAN, PDIST.

% Reference: Speed, T. Statistical Analysis of Gene Expression
%            Microarray Data (2003), Chapman & Hall

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.1.12.1 $   $Date: 2004/12/24 20:42:47 $

if nargin < 2
    K = 1;
end
K = min(K,size(data,2));

metric = 'euclidean';
uWeights = [];
userWeights = false;
distargs = {};
useWMean = true;

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'distance','distargs','median','weights'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        karg = strmatch(lower(pname), okargs); %#ok
        if isempty(karg)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(karg)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(karg)
                case 1  % distance
                    metric = pval;

                case 2 %distfun args
                    if ~iscell(pval)
                        distargs = {pval};
                    else
                        distargs = pval;
                    end

                case 3 % median flag
                    medianFlag = opttf(pval);
                    if isempty(medianFlag)
                        error('Bioinfo:InputOptionNotLogical','%s must be a logical value, true or false.',...
                            upper(char(okargs(karg))));
                    end
                    if medianFlag
                        useWMean = false;
                    end
                case 4  % weights
                    uWeights = pval(:);
                    userWeights = true;
            end
        end
    end
end

% create a copy of data for output
imputed = data;

% identify missing vals
nanVals = isnan(data);
% use rows without nans for calculation of nearest neighbors
noNans = sum(nanVals,2) == 0;

dataNoNans = data(noNans,:);
if isempty(dataNoNans)
    error('Bioinfo:KNNIMPUTE:AllNans',...
        'All rows of the input data contains missing values. Unable to impute missing values.');
end

% check weights
if userWeights
    if numel(uWeights) ~= K || ~isnumeric(uWeights)
        error('Bioinfo:KNNIMPUTE:BadWeights',...
            'WEIGHTS must be a numeric vector of length K.');
    end
    if ~useWMean
        warning('Bioinfo:KNNIMPUTE:MedianWeights',...
            'WEIGHTS are ignored when using MEDIAN.');
    end
end
% calculate pairwise distances between columns
distances = pdist(dataNoNans',metric,distargs{:});

% sort and get indices of nearest columns
SqF = squareform(distances);
% force the diagonals to be negative so that they automatically sort to the
% top.
[dists, ndx] = sort(SqF - eye(size(SqF)));
dists(1,:) = 0;
equalDists = [zeros(1,size(dists,2));(diff(dists)== 0);];
% map ndx back to columns in data

% get rows and columns of missing values
[rows,cols] = find(nanVals);

% for each missing value find the nearest column without a corresponding
% NaN.

% it is possible that we could get divide by zero warnings in here.
ws = warning('off','MATLAB:divideByZero');


rowWarn = false(size(data,1),1);
for count = 1:numel(rows)
    % check that we don't have a row of nans
    if all(isnan(data(rows(count),:)))
        if rowWarn(rows(count)) == false
            warning('Bioinfo:KNNIMPUTE:RowAllNans',...
                'Row %d contains all NaNs.',rows(count));
            rowWarn(rows(count)) = true;
        end
        continue
    end
    % start at 2 as we know that each col is the closest element
    for nearest = 2:size(ndx,1)-K+1
        % Look for exactly equal columns and add any of these into the
        % calculation.
        L = max(find(equalDists(nearest+K-1:end,cols(count))==0,1)-2,0);
        dataVals = data(rows(count),ndx(nearest:nearest+K-1+L,cols(count)));
        if useWMean
            if ~userWeights
                weights = 1./dists(2:K+L+1,cols(count));
            else
                if L > 0
                    weights = [uWeights(1:end-1); repmat(uWeights(end)/L,L,1)];
                else
                    weights = uWeights;
                end
            end
            val = wnanmean(dataVals,weights);
        else
            val = nanmedian(dataVals);
        end
        if ~isnan(val)
            imputed(rows(count),cols(count)) = val;
            break
        end
    end
end
% reset the warning state
warning(ws);


function m = wnanmean(x,weights)
%WNANMEAN Weighted Mean value, ignoring NaNs, infs are special

% Find NaNs and set them to zero
x = x(:); weights = weights(:);
nans = isnan(x);
infs = isinf(weights);
if sum(~nans) == 0
    m = NaN;
    return
end
if sum(infs) > 0
    m = nanmean(x(infs));
    return 
end
% normalize the weights
weights = weights./sum(weights);
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = nansum(x.*weights);
