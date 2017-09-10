% Compute the auto-correlation at many distances 
% stratified according to the two nucleotides. The method is by calling
% 16 times to the standard autocorrelation function and normalizing
% accordingly. 
function [auto_corr_mat auto_corr_counts] = xcorr_by_nuc(seq, weights, max_lag)

% epsilon = 0.0000000001; % used to avoid division by zero 

if(~isnumeric(seq))
    seq = nt2int(seq);
end
n = length(seq); nuc_inds = zeros(4,n); 
for i=1:4 % first nuc
    nuc_inds(i,:) = double(seq == i);
end

max_lag = min(max_lag, n); auto_corr_mat = zeros(4, 4, 2*max_lag+1); auto_corr_counts = auto_corr_mat;
for i=1:4
    for j=1:4 % second nuc
        auto_corr_counts(i,j,:) = xcorr(nuc_inds(i,:), nuc_inds(j,:), max_lag); % how many times we've got nuc i and nuc j at given distances 
        auto_corr_mat(i,j,:) = xcorr(weights .* nuc_inds(i,:), weights .* nuc_inds(j,:), max_lag);
    end
end
        
        


