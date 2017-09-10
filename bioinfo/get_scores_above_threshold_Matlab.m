% Get all the binding sites for a given motif which are
% above a certain threshold. It uses as its input the output of the function 'calc_all_sites_scores'
% It support both kinds of thresholding : Absolute (given value) or fractional
% (given fraction)
% We also output the threshold since it may be used for other blocks
% This function doesn't do any actual work, but serves as an interface to
% the c function (mex file) with the same name
%
% Input:
% loc_scores - matrix/cell-array of location scores
% loc_scores_length - vector of lengths of the different regions
% threshold - the threshold above which we consider a 'hit'
% L - the length of the Binding Site
% smart_threshold - the function can determine the threshold herself.
%   The options for smart threshold are the following:
%       -1 - pick the singel best instance of loc_scores (in each member of the cell array)
%       0 - pick the input threshold given from outside : 'threshold'
%       1 - choose threshold in a smart way. Treat 'threshold' as the top fraction above which we keep sites
%       2 - choose threshold in a smart way, such that exactly k scores pass, where k is the 'threshold' parameter
% diff_lengths_flag - saying if the regions are of different lengths
%
% Output :
% BS_genes - indices to genes containing the thresholds
% BS_places - places where the BS were found
% BS_scores - scores of the found BS
% thresh_output - the 'optimal' threshold as computed by the function (if it was required to do so)
%
function [BS_genes BS_places BS_scores thresh_output] = ...
    get_scores_above_threshold_Matlab(loc_scores, loc_scores_lengths, threshold, L, smart_threshold, diff_lengths_flag)


switch smart_threshold
    case -1 % this means that we pick only the best score in each promoter
        n = length(loc_scores_lengths);
        BS_genes = [1:n];
        if(iscell(loc_scores))
            [BS_scores BS_places] = max_cell(loc_scores);
        else
            [BS_scores BS_places] = max(loc_scores);
        end
        thresh_output = min(BS_scores); % give the 'min-max' score
        return; % here we don't need to call the c function
    case 0 % 'no real threshold' - just read the input threshold
        thresh_output = threshold;
    case 2  % pass only 'threshold' scores (so here 'threshold' is a positive integer)
        if(iscell(loc_scores))
            thresh_output = my_quantile(cell2vec(loc_scores), 1 - threshold ./ sum(loc_scores_lengths));
        else
            thresh_output = my_quantile(loc_scores(:), 1 - threshold ./ loc_scores_lengths(1));
        end
    otherwise % use smart threshold. This includes case 1
        if(iscell(loc_scores))
            thresh_output = my_quantile(cell2vec(loc_scores), 1-threshold); % get_gene_scores_quantile_threshold(loc_scores,loc_scores_lengths, threshold, diff_lengths_flag);
        else
            thresh_output = my_quantile(loc_scores(:), 1-threshold); % get_gene_scores_quantile_threshold(loc_scores,loc_scores_lengths, threshold, diff_lengths_flag);
        end
end

[BS_genes BS_places BS_scores] = get_scores_above_threshold_mex(loc_scores, loc_scores_lengths, ...
    thresh_output, L, smart_threshold, diff_lengths_flag); % call the .mex function interface to decide which c function to run
if(smart_threshold == 2) % take top k
    n = length(BS_genes);
    BS_genes = BS_genes(1:min(n,threshold));
    BS_places = BS_places(1:min(n,threshold));
    BS_scores = BS_scores(1:min(n,threshold));
end




% An helper function that helps to coordinate the calling to the c function
% It decides based on the type of the input if to calls the mex c function
% with singles or doubles, and also if the input is a cell array it splits
% it into its components
function [BS_genes BS_places BS_scores] = ...
    get_scores_above_threshold_mex(loc_scores, loc_scores_lengths, threshold, L, smart_threshold, diff_lengths_flag)

if(iscell(loc_scores))
    n = length(loc_scores);
    BS_genes = cell(n,1); BS_places = cell(n,1); BS_scores = cell(n,1);
    if(length(threshold) == 1)
        threshold = repmat(threshold, n, 1);
    end
    for i=1:n
        [BS_genes{i} BS_places{i} BS_scores{i}] = ...
            get_scores_above_threshold_mex(loc_scores{i}, loc_scores_lengths(i), threshold(i), L, smart_threshold, diff_lengths_flag);
        if(~isempty(BS_genes{i})) % set the genes indices according to the cell locations
            BS_genes{i}(:) = i;
        end
    end
    BS_genes = cell2vec(BS_genes); BS_places = cell2vec(BS_places); BS_scores = cell2vec(BS_scores);
else
    %    loc_scores = double(loc_scores); % Temp!!! convert to doubles because we don't yet have a 'singles' implementation !!!
    is_double = double(isa(loc_scores, 'double')); % check type of scores
    if(is_double)
        threshold = double(threshold);
    else
        threshold = single(threshold);
    end
    if(diff_lengths_flag == 0)
        [BS_genes BS_places BS_scores] = ...
            get_scores_above_threshold(loc_scores, threshold, L, smart_threshold, is_double);  % Note: the function doesn't actually use the smart_threshold flag
    else
        [BS_genes BS_places BS_scores] = ...
            get_scores_above_threshold_diff_lengths(loc_scores, loc_scores_lengths, ...
            threshold, L, smart_threshold, is_double);  % 0 means 'single'. Note: the function doesn't actually use the smart_threshold flag
%          [BS_genes2 BS_places2 BS_scores2] = ...
%             get_scores_above_threshold_diff_lengths(double(loc_scores), loc_scores_lengths, ...
%             double(threshold), L, smart_threshold, is_double+1);  % 0 means 'single'. Note: the function doesn't actually use the smart_threshold flag
   end
end

