% Normalization - Subtracts the mean of each row and divide by standard deviation (from Tal Shay)
%
% Input: 
% data - matrix
% orientation_flag - which dimension to normalize by
%
% Output:
% normalized_data - data after subtracting mean and dividing by st.d.
%
function normalized_data = NormalizeData(data, orientation_flag, varargin)
ROW=1;
if(~exist('orientation_flag', 'var') || (orientation_flag == ROW))
    ncols = size(data,2);
    row_std_mat = repmat(nanstd(data,[],2),1,ncols);
    row_mean_mat = repmat(nanmean(data,2),1,ncols);
    normalized_data = (data - row_mean_mat)./row_std_mat;
else
    nrows = size(data,1);
    col_std_mat = repmat(nanstd(data,[],1),nrows,1);
    col_mean_mat = repmat(nanmean(data,1),nrows,1);
    normalized_data = (data - col_mean_mat)./col_std_mat;
end

