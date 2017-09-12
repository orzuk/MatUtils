function Pmean = PolishMean(data, t, varargin)
% Calculate mean after without lower/upper t%. (from Tal Shay)
% Default t is 2%. 
if nargin == 1
    t = 2;
end

low_thr = prctile(data,t);
high_thr = prctile(data,100-t);

n = size(data, 2);
Pmean = zeros(n, 1);
for i = 1:n
    A = data(:, i);
    truncA = A(find(A > low_thr(i) & A < high_thr(i)));
    Pmean(i) = mean(truncA);
end