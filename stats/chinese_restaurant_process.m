% Gerenate an instance of chinese restaurant process
% 
% Input: 
% n - number of people 
% num_points - how many instances to generate
% alpha - process param (default is zero)
% theta - process param (default is one)
% 
% Output: 
% s - matrix of generated partitions of n 
% 
function s = chinese_restaurant_process(n, num_points, alpha, theta)

if(~exist('alpha', 'var') || isempty(alpha)) % set default parameters 
    alpha = 0; 
end
if(~exist('theta', 'var') || isempty(theta))
    theta = 0; 
end

% either ? < 0 and ? = L? for some L ? {1, 2, ...}; or that 0 ? ? ? 1 and ? > ??.
if(alpha < 0)     % test legal parameters
    L = theta / alpha;
    if(abs(round(L) - L) > 0.000000001)
        error('Error! Illegal \alpha nad \theta parameters');
    end
else
    if( (alpha > 1) || (theta < -alpha) )
        error('Error! Illegal \alpha nad \theta parameters');
    end
end
% s = 1;  % s holds the state. s(i) is the $people in the i'th table
s = zeros(num_points,n); s(:,1) = 1; % Change: s(i) is the number of sets of size i. i=1:n
tables = ones(num_points,1);

for i=2:n % loop and add person
    empty_weights = (theta + tables.*alpha) ./ (theta + i - 1);
%    full_weights =  (s - alpha) ./ (theta + i - 1);
    full_weights = s .* (repmat(1:n, num_points, 1) - alpha) ./ (theta + i - 1); 
    
    next_table = weighted_rand([full_weights empty_weights], 1); % randomize the next table

    new_table_inds = find(next_table > tables);     old_table_inds = find(next_table <= tables);  
    tables(new_table_inds) = tables(new_table_inds)+1; 
    s(new_table_inds,1) = s(new_table_inds,1) + 1; 
    

    for j=1:length(old_table_inds) % loop on old inds 
        s(old_table_inds(j), next_table(old_table_inds(j))) = ...
            s(old_table_inds(j), next_table(old_table_inds(j)))-1;
        s(old_table_inds(j), next_table(old_table_inds(j))+1) = ...
            s(old_table_inds(j), next_table(old_table_inds(j))+1)+1;
    end
    
% %     if(next_table <= tables) % add person to a table 
% %         % s(next_table) = s(next_table)+1;
% %         s(next_table) = s(next_table)-1; 
% %         s(next_table+1) = s(next_table+1)+1; 
% %     else
% %         % s(end+1) = 1; tables = tables+1; 
% %         s(1) = s(1) + 1; tables = tables+1; 
% %     end
end

% s = hist(s, 1:n); % collect tables to a partition

