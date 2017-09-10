% Simple 1d Gaussian smoothing of a histogram. 
% Note: the function changes the length of the vector (according to sigma)
% 
% Input: 
% v - the probability values (assumed to be equaly spaced!)
% sigma - width of Gaussian smoothing
% chop_flag - chops the output to match the input's length
% 
% Output: 
% v_smoothed - the smoothed version of v 
% 
function v_smoothed = gauss_smooth(v, sigma, chop_flag, varargin)
%v_smoothed = v; 
%return; % TMP BAD!

if(~exist('chop_flag', 'var'))
    chop_flag = 0;
end
k = 50*sigma; % width should be enough (misses < 10^(-6) of histogram)
v_gauss = normpdf(  (-k:k)  ./ sigma ); % Compute the Gaussian vector 
v_gauss = v_gauss ./ sum(v_gauss); % normalize 

v_smoothed = conv(v_gauss, v); % Perform convolution with Gaussian density
if(chop_flag)
    k = length(v_smoothed) - length(v);
    v_smoothed = v_smoothed(floor(k/2):floor(k/2)+length(v)-1);
%     figure; hold on;
%     plot(v, '.'); plot(v_smoothed, 'r.'); legend('orig. data', 'smoothed');
end


