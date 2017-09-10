% Take the original (true) corrleations, and add to them random gaussian/binomial noise.
% We deal with a vector of different lengthes simultaniousely 
function samp_corr_vec = sample_correlation(corr_vec, nsamples_vec, gauss_flag)


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

    % normalize to have mean C_i and std sqrt(C_i(1-C_i) / N)
%%%    samp_corr_vec = (rand_gauss_mat .* (repmat(sqrt(corr_vec .* (1-corr_vec)),length(nsamples_vec), 1) ./ ...
%%%        repmat(sqrt(nsamples_vec'), 1, length(corr_vec)))    ) + repmat(corr_vec,length(nsamples_vec), 1);
    
    
    % Constant variance approx.
        samp_corr_vec = (rand_gauss_mat .* 0.5 ./ repmat(sqrt(nsamples_vec'), 1, length(corr_vec))) + repmat(corr_vec,length(nsamples_vec), 1);
    
end