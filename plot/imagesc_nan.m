% Display an array with nans as white squares
function h = imagesc_nan(data, t)

epsilon = 0.000000000001; % used when all data is the same 
if (nargin == 2)
    h = imagesc(data, t);
else
    data_min = min(min(data)); data_max = max(max(data)); 
    data_gap = max(epsilon, data_max-data_min);
    h = imagesc(data, [data_min-data_gap/32 data_max+data_gap/32]);
end

cmap = colormap; cmap(1,:) = [1 1 1]; colormap(cmap);

% [a b] = find(isnan(data));
% n = length(a);
% for i = 1:n
%     rectangle('position', [b(i) - 0.5, a(i) - 0.5, 1, 1], 'facecolor', 'w', 'edgecolor', 'w');
% end
% return
