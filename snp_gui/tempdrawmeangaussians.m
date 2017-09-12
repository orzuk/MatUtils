%x = 20;

num_good_inds = 13;

[mean_mat std_mat std_inv_mat] = ExpandSNPsMoments(RLMM, CHROM_MATS.snp_ids{chr}, r_mat);

good_inds = find(sum(mean_mat) > 10);
figure; hold on; x_ctr = 1;
for x=x_vec

    subplot(3,3,x_ctr); hold on; x_ctr = x_ctr+1;
    MuVecs = [mean_mat(x,good_inds(1:num_good_inds))' mean_mat(x,good_inds(num_good_inds+1:2*num_good_inds))'];
    s =  [std_mat(x,good_inds(1:num_good_inds))' std_mat(x,good_inds(num_good_inds+1:2*num_good_inds))' std_mat(x,good_inds(1:num_good_inds)+50)' ];

    for j=1:num_good_inds
        SigmaMats{j} = [s(j,1) s(j,3); s(j,3) s(j,2)];
    end

    labels_vec = {['SNP # ' num2str(x) ' chr ' num2str(chr) ' loc ' num2str(CHROM_MATS.loc{chr}(x)) ' ID ' CHROM_MATS.snp_ids{chr}{x} ], 'Y'};
    for j=1:num_good_inds
        legends_vec{j} = num2str(j);
    end
%%    legends_vec = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};

    %figure;
    MixtureOfGaussiansDraw2dGaussians(MuVecs, SigmaMats, labels_vec, []); % , legends_vec);

    hold on; plot(CHROM_MATS.data_A{chr}(x), CHROM_MATS.data_B{chr}(x), '+y');
    plot(CHROM_MATS.data_A{chr}(x), CHROM_MATS.data_B{chr}(x), 'oy');
    plot(CHROM_MATS.data_A{chr}(x), CHROM_MATS.data_B{chr}(x), '*y');
    plot(CHROM_MATS.data_A{chr}(x), CHROM_MATS.data_B{chr}(x), 'sy');
    plot(CHROM_MATS.data_A{chr}(x), CHROM_MATS.data_B{chr}(x), 'dy');

end
