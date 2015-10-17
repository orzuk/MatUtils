% Test if a distribution is bimodal 
% Input: 
% x - values
% p - their probabilities
% err - relative error in bi-modality (we allow small fluctuations)
% 
% Output: 
% mu - mean value 
% 
function isbimodal = isbimodal_hist(x, p, relative_error)

p = normalize_hist(x, p); % normalize to sum to one 

p_diff = diff(p); % compute diff

s = sign(p_diff); % see when derivative changes sign
sign_changes = [1 1+find(s(1:end-1) ~= s(2:end))]; 

for i=1:length(sign_changes) % these are the local maxima
    maxima_height(i) = max(p(sign_changes(i)) - p(sign_changes(i)+1:end)); 
end

if(length(maxima_height) > 1)
    maxima_height = sort(maxima_height);
    % isbimodal = any(p_diff < -relative_error); %isbimodal = any(p_diff ./ p(1:end-1) < -relative_error);
    isbimodal = maxima_height(end-1) > relative_error;
    
else % just one maximum
    isbimodal = 0;
end

% figure; hold on; plot(x(1:end-1), p(1:end-1), 'r'); plot(x(1:end-1), p_diff); 
% plot(x(sign_changes), p(sign_changes), 'kX', 'linewidth', 3); 
% title('Test for Bimodality'); 
% xlabel('x'); ylabel('freq.'); 
% legend('prob. density', 'derivative'); 

