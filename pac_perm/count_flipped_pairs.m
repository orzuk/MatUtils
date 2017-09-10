% Count how many pairs have switched 
% Do it matrix way ...
function  flipped_pairs_vec = count_flipped_pairs(sample_corr_vec)

%  repmat(sample_corr_vec',  1, length(sample_corr_vec)) 
% tril(repmat(sample_corr_vec, length(sample_corr_vec), 1))
%flipped_pairs  = sum(  sum(  repmat(sample_corr_vec', 1,  length(sample_corr_vec)) <  tril(repmat(sample_corr_vec,  length(sample_corr_vec), 1) ) ) );

% Alternative : forloop
flipped_pairs_vec  = zeros(size(sample_corr_vec, 1), 1); 
% loop over indices 
% size(sample_corr_vec, 1)
% size(sample_corr_vec, 2)
% sample_corr_vec(1,:)
for i=1:size(sample_corr_vec, 2)
%     size(repmat(sample_corr_vec(:,i), 1, size(sample_corr_vec,2)-i))
%     size(sample_corr_vec(:,i+1:end))
%     size( flipped_pairs_vec)
%     size(sum(  repmat(sample_corr_vec(:,i),  1, size(sample_corr_vec,2)-i) > sample_corr_vec(:,i+1:end), 2))
    flipped_pairs_vec  = flipped_pairs_vec  + sum(  repmat(sample_corr_vec(:,i),  1, size(sample_corr_vec,2)-i) > sample_corr_vec(:,i+1:end), 2);
end