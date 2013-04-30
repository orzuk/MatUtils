% Now function sample_correlation
% We deal with a vector of different lengthes simultaniousely 
% The function takes the original (true) corrleations, and adds to them
% random gaussian/binomial noise.
function samp_corr_vec = sample_noisy_correlations(corr_vec, nsamples_vec, gauss_flag, sig_sig, coeffs_vec, corr_mat_sqrt)

TOL = 0.000000000000001;

if(gauss_flag == 0) % use binomial
    % generate binomial distributions from correlation. We assume that the
    % nsample_vec is sorted !!!!
    rand_mat = rand( nsamples_vec(end), length(corr_vec)); 
    
    samp_corr_vec = zeros(length(nsamples_vec), length(corr_vec)); 
    
    j=1;
    for nsamples = nsamples_vec
        samp_corr_vec(j,:) = mean(   rand_mat(1:nsamples,:) < repmat(corr_vec, nsamples, 1) );
        j=j+1;
    end
    
else % use gaussians

    % A different approach : Generate gaussian variables (approx. of binomials)
    % - Faster and breaks ties !!!!
    rand_gauss_mat = randn(length(nsamples_vec), length(corr_vec)); % Gaussian r.v.s ( much fewer)
    
    
    

    % Sort the genes such that the correlation matrix will be in the
    % correct order
    [sorted_vec sort_inds ] = sort(corr_vec);

  %%%%%%%%%%%%%%%%%  rand_gauss_mat = (corr_mat_sqrt(sort_inds,sort_inds) * rand_gauss_mat')'; % Generate correlations in the noise !!!  
    
         % Check moments after transformation
    
    % normalize to have mean C_i and std sqrt(C_i(1-C_i) / N)
%%%    samp_corr_vec = (rand_gauss_mat .* (repmat(sqrt(corr_vec .* (1-corr_vec)),length(nsamples_vec), 1) ./ ...
%%%        repmat(sqrt(nsamples_vec'), 1, length(corr_vec)))    ) + repmat(corr_vec,length(nsamples_vec), 1);
    
    
    % Constant variance approx. Why the half here ???? Currently Removed
    % !!!! 4.5.05 !!!
%    nsamples_vec
%    sissssissgisgsigisgisgig = sqrt(nsamples_vec)
%%%% ????    samp_corr_vec = (rand_gauss_mat .* 1 ./ repmat(sqrt(nsamples_vec'), 1, length(corr_vec))) + repmat(corr_vec,length(nsamples_vec), 1);
    
    
    
    % New : We let the st.d. of the noise itself to be gaussian !!! 
%     sig_mat = randn(length(nsamples_vec), length(corr_vec)) * sig_sig + 1;     
%     sig_mat = max(sig_mat, TOL);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % %     % New : Let the st.d. be  a quadratic function of the mean
% % % % % % % % % %     sig_mat = coeffs_vec(1) .* corr_vec.^2 + coeffs_vec(2) .* corr_vec + coeffs_vec(3); 
% % % % % % % % % %     % Now sample the noise  
% % % % % % % % % %     rand_gauss_mat = rand_gauss_mat .* sig_mat; 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    
    samp_corr_vec = (rand_gauss_mat .* 1 ./ repmat(sqrt(nsamples_vec'), 1, length(corr_vec))) + repmat(corr_vec,length(nsamples_vec), 1);

%%%    samp_corr_vec = rand_gauss_mat + repmat(corr_vec,length(nsamples_vec), 1);
    
end

