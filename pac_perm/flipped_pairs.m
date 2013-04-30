% Count how many pairs have switched 
function  flipped_pairs = count_flipped_pairs(corr_vec, sample_corr_vec)

flipped_pairs  = sum(  sum(  repmat(sample_corr_vec',  length(sample_corr_vec), 1) >  tril(repmat(sample_corr_vec, 1, length(sample_corr_vec)) ) ) );

