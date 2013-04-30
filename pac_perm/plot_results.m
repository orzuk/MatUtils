% A script for plotting the overlap results of the data 
N_genes=size(R.dat,1);
alpha=0.012;
N_TOP=alpha*N_genes;
res=4;
min_i=10;
max_i=48;
group_sizes=[min_i:res:max_i];
figure,errorbar(group_sizes,f,var_f,'k');
hold on,errorbar(nsamples_vec, samp_frac_mean(1,:), samp_frac_std(1,:), 'r');
hold on,errorbar(nsamples_vec, inf_limit_frac(1,:), ...
    sqrt(inf_limit_frac(1,:).*(1-inf_limit_frac(1,:))/N_TOP),'g')

figure,errorbar(group_sizes,f,var_f,'k');
hold on,plot(nsamples_vec, samp_frac_mean(1,:),'r')
hold on, plot(nsamples_vec, inf_limit_frac(1,1:10));
hold on,plot(nsamples_vec, samp_frac_mean_from_data(1,:),'r')

figure;  num_bins = max(250, floor(N_genes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(SAVE_FISHER_Z{data_flag}, num_bins);
    hist(Fisher_Zs{data_flag}, num_bins); xlabel('Correlation'); ylabel('Freq.'); 
    hold on,plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag}.* N_genes .* (bin_pos(2)-bin_pos(1)), 'r');

figure, num_bins = max(250, floor(Ngenes/250)); % % num. of bins is arbitrary. Should be determined more cleverly
    [H bin_pos] = hist(rand_gauss_samp, num_bins);
    hist(rand_gauss_samp, num_bins); xlabel('Correlation'); ylabel('Freq.'); title([data_str ' Data Hist. of rand. samp.']);
    hold on, plot(Fisher_Zs_possible_vec, SAVE_FISHER_NORMAL_FIT_VEC{data_flag} .* Ngenes .* (bin_pos(2)-bin_pos(1)), 'r');
