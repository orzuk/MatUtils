% Interface calling Kangwon Lee's 2d histogram functions
% We call hist2d to compute the histogram and Plot2dHist to draw it
%
% Input: 
% mX - 2d data (2-by-n matrix)
% vYEdge - Y bins
% vXEdge - X bins
% strLabelX - label for X axis
% strLabelY - label for Y axis
% strTitle - title string
% trunc_val - a maximal value of the histogram we can trancate by 
%
% Output: 
% mHist - 2d histogram matrix 
%
function mHist = hist2d_draw(mX, vYEdge, vXEdge, strLabelX, strLabelY, strTitle, trunc_val, colorbar_flag, varargin) 

if(length(mX) > 10000000) % set a maximal size to plot - sample above this value 
    P = randperm(length(mX)); 
    mX = mX(P(1:10000000), :);
end
if(length(vXEdge) == 1) % input is just num of bins
    vXEdge = range(mX(:,2)) / vXEdge;
    vXEdge = min(mX(:,2)):vXEdge:max(mX(:,2));
end
if(length(vYEdge) == 1) % input is just num of bins
    vYEdge = range(mX(:,1)) / vYEdge;
    vYEdge = min(mX(:,1)):vYEdge:max(mX(:,1));
end

mHist = hist2d (mX, vYEdge, vXEdge);
if(~exist('strLabelX', 'var'))
    strLabelX = '';
end
if(~exist('strLabelY', 'var'))
    strLabelY = '';
end
if(~exist('strTitle', 'var'))
    strTitle = '';
end
if((exist('trunc_val', 'var')) & (~isempty(trunc_val)))
    mHist = min(mHist, trunc_val); 
end
if(~exist('colorbar_flag', 'var'))
    colorbar_flag = 1;
end
Plot2dHist(mHist, vXEdge, vYEdge, strLabelX, strLabelY, strTitle, colorbar_flag);

