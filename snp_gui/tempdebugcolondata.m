num_bins = 500; sam = 35;
COLON_DATA = (COLON.NormalizedSNPsCopyMatA(:,sam)+ COLON.NormalizedSNPsCopyMatB(:,sam)) .* ...
    COLON.NormalMeanIntensities(sam) ./ 2;
[SS,MM,PP, DimDim, LogLike] = MixtureOfGaussiansFindModelDimension(COLON_DATA, 4, 50, 1);

legends_vec = [];

MixtureOfGaussiansDraw1dGaussians( COLON_DATA, COLON.mix_gauss_params.P{sam}, COLON.mix_gauss_params.M{sam}, ...
    (COLON.mix_gauss_params.S{sam}), {['Sample # ' num2str(sam)], ''}, legends_vec, [], [], num_bins)

MixtureOfGaussiansDraw1dGaussians( COLON_DATA, PP, MM, SS, {['NEW MOG Sample # ' num2str(sam)], ''}, legends_vec, [], [], num_bins)


%% [OutputFilesNames NormalMeanIntensities r_mat] = PreLearnMoGChromFromData('E:\Research\HMM_Chromosome\zuk_hmm\hmm_chrom\data\colon', 'xba');

