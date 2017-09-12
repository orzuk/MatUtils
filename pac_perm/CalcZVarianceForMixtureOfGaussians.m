% Note : Here X is binary and Y is gaussian !!!

% Here we load a data, compute the correlations and get the desired
% fraction
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');


TOL = 0.00000000001;

% flags 
blocks_flag = 1; gaussian_flag = 1; % 0 - MOG 1 - Two Gaussians

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
rand_flag = GAUSSIAN;

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TRUE_AND_SAMPLED;

max_nsamples = 200; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)

res = 100; % the resulution of N that we use
alpha = 0.014; alpha_vec =  [0.012  0.12]; % corresponds to ~70 and ~700 genes

%nsamples=3;
iters = 1000;

nsamples_vec = res:res:max_nsamples;

Ncorrs = 50;
Nsamples =48;
res = 1/Ncorrs;

corrs_vec = [-1+res:res:1-res];
Fisher_corrs_vec =  0.5 * (log(1+corrs_vec) - log(1-corrs_vec));

% Compute the sigma for the MOGs
sigma_vec = corrs_vec;
bad_ind = find(sigma_vec == 0)
Ncorrs

sigma_vec(Ncorrs) = TOL; sigma_vec = sqrt(1-sigma_vec.^2) ./ sigma_vec;
Ncorrs = length(corrs_vec);


Rho_sig_mean = zeros(1, Nsamples); Rho_sig_std = zeros(1, Nsamples);
Rho_bias_mean = zeros(1, Nsamples); Rho_bias_std = zeros(1, Nsamples);
Z_sig_mean = zeros(1, Nsamples); Z_sig_std = zeros(1, Nsamples);
Z_bias_mean = zeros(1, Nsamples); Z_bias_std = zeros(1, Nsamples);




% Prepare the A matrix transforming to correlated Gaussians
AAA = zeros(2,2*Ncorrs);
plus_sqrt = sqrt(1+corrs_vec); minus_sqrt = sqrt(1-corrs_vec);
AAA(1,1:2:end-1) = plus_sqrt+minus_sqrt;
AAA(1,2:2:end) = plus_sqrt-minus_sqrt;
AAA(2,1:2:end-1) = plus_sqrt-minus_sqrt;
AAA(2,2:2:end) = plus_sqrt+minus_sqrt;
AAA = AAA ./ 2; AAA=AAA';




N_VEC = [5:1:Nsamples];
Ncorrs = length(corrs_vec);
for N=N_VEC

    if(blocks_flag == 0)
        num_blocks = 1;
    else
        num_blocks = floor(2*Nsamples/N);
    end

    N_IS = N
    Pearson_Corrs =zeros(Ncorrs,1);

    Pearson_Rhos_mean=zeros(Ncorrs,1);  Pearson_Rhos_std=zeros(Ncorrs,1);
    Fisher_Zs_mean=zeros(Ncorrs,1);   Fisher_Zs_std=zeros(Ncorrs,1);

    for i=1:iters

        % If we estimate using the blocks method, we have here an extra
        % loop on the blocks :

        chunk_Pearson_corrs = zeros(Ncorrs,num_blocks);
        chunk_Fisher_Zs = zeros(Ncorrs,num_blocks);
        for k=1:num_blocks


            if(gaussian_flag)
                %   randomize data
                Z = randn(2,N);
                X = AAA*Z;
                % Now calculate correlations.

                %Normalized data
                X = X - repmat(mean(X,2), 1, N);
                X = X ./ repmat(sqrt(sum(X.^2,2)), 1, N);

                chunk_Pearson_corrs(:,k) = sum(X(1:2:end-1,:) .* X(2:2:end,:), 2);

            else

                %   randomize data
                X = 2.*(rand(1,N) > 0.5)-1;
                while(min(X) == max(X)) % avoid equals ...
                    X = 2.*(rand(1,N) > 0.5)-1;
                end

                Z = randn(1,N);  % Get the Z's

                % Multiply by the sigma vec

                Z = Z' * sigma_vec;
                Z = Z + repmat(X', 1, Ncorrs);

                % Now calculate correlations.

                %Normalized data
                X = X - mean(X);
                X = X ./ (sqrt(sum(X.^2)));

                Z = Z - repmat(mean(Z,1),N,1);
                Z = Z ./ repmat(sqrt(sum(Z.^2,1)),N,1);

                % Adjust X by flipping it whe
%                 size(Z)
%                 size(repmat(-X', 1, 1/res-1))
%                 size(repmat(-X', 1, 1/res))
%                 size([(repmat(-X', 1, 1/res-1)) (repmat(X', 1, 1/res))] )
%                 size(chunk_Pearson_corrs)
%                 size(chunk_Pearson_corrs(:,k))
%                 size(sum(Z .* [(repmat(-X', 1, 1/res-1)) (repmat(X', 1, 1/res))]  ,1))
                chunk_Pearson_corrs(:,k) = sum(Z .* [(repmat(-X', 1, 1/res-1)) (repmat(X', 1, 1/res))]  ,1)';
            end

            chunk_Fisher_Zs(:,k) =  0.5 * (log(1+chunk_Pearson_corrs(:,k)) - log(1-chunk_Pearson_corrs(:,k)));

            Pearson_corrs = chunk_Pearson_corrs(:,k); %%sum(Z .* [(repmat(X', 1, 1/res-1)) (repmat(-X', 1, 1/res))]  ,1); Pearson_corrs=Pearson_corrs';
            Fisher_Zs = 0.5 * (log(1+Pearson_corrs) - log(1-Pearson_corrs));
        end

        if(blocks_flag == 0)
            Pearson_Rhos_mean = Pearson_Rhos_mean +  Pearson_corrs;
            Pearson_Rhos_std = Pearson_Rhos_std +  Pearson_corrs.^2;
            Fisher_Zs_mean = Fisher_Zs_mean + Fisher_Zs;
            Fisher_Zs_std = Fisher_Zs_std + Fisher_Zs.^2;

        else
%             size(Pearson_Rhos_mean)
%             size(chunk_Pearson_corrs)
%             size( mean(chunk_Pearson_corrs,1))            
%             size( mean(chunk_Pearson_corrs,2))
            Pearson_Rhos_mean = Pearson_Rhos_mean +  mean(chunk_Pearson_corrs,2);
       %     size(std(chunk_Pearson_corrs,0,2))
            Pearson_Rhos_std = Pearson_Rhos_std +  std(chunk_Pearson_corrs,0,2).^2;
            Fisher_Zs_mean = Fisher_Zs_mean + mean(chunk_Fisher_Zs,2);
            Fisher_Zs_std = Fisher_Zs_std + std(chunk_Fisher_Zs,0,2).^2;
        %    size(std(chunk_Fisher_Zs,0,2))
        end

    end   % iterations loop


    Pearson_Rhos_mean = Pearson_Rhos_mean ./ iters;
    Pearson_Rhos_std = Pearson_Rhos_std ./ iters; %iters

    if(blocks_flag == 0)
        Pearson_Rhos_std  = Pearson_Rhos_std  - Pearson_Rhos_mean.^2;
    end
    Pearson_Rhos_std = sqrt(Pearson_Rhos_std);  % transfer var to std
    
    Rho_sig_mean(N) = mean(Pearson_Rhos_std);   Z_sig_std(N) = std(Pearson_Rhos_std);
    Rho_bias_mean(N) = mean(Pearson_Rhos_mean-corrs_vec'); Rho_bias_std(N) = std(Pearson_Rhos_mean-corrs_vec');


    Fisher_Zs_mean = Fisher_Zs_mean ./ iters;
    Fisher_Zs_std = Fisher_Zs_std ./ iters;
    if(blocks_flag == 0)
        Fisher_Zs_std  = Fisher_Zs_std  - Fisher_Zs_mean.^2;
    end
    Fisher_Zs_std = sqrt(Fisher_Zs_std); % transfer var to std

    Z_sig_mean(N) = mean(Fisher_Zs_std);   Z_sig_std(N) = std(Fisher_Zs_std);
    Z_bias_mean(N) = mean(Fisher_Zs_mean-Fisher_corrs_vec'); Z_bias_std(N) = std(Fisher_Zs_mean-Fisher_corrs_vec');

end

plot_str = [];
plot_str{1} = 'MOG ';
plot_str{2} = 'MOG w. blocks';
plot_str{3} = 'Gaussians ';
plot_str{4} = 'Gaussians w. blocks';


plot_ind = 1+blocks_flag + 2*gaussian_flag;

% Now plot results
figure; subplot(2,2,1);  hold on; title(['Pearson P std, for ' plot_str{plot_ind}]);
errorbar(N_VEC, Rho_sig_mean(N_VEC),Rho_sig_std(N_VEC));
xlabel('N'); ylabel('P std'); %legend('from data', 'from Fisher');


subplot(2,2,2);  hold on; title(['Fisher Z std, for '  plot_str{plot_ind} ' and analytic']);
plot(N_VEC, 1./sqrt(N_VEC-3), 'r');
errorbar(N_VEC, Z_sig_mean(N_VEC),Z_sig_std(N_VEC));

xlabel('N'); ylabel('Z std'); legend( '(N-3)^{-1/2}','from data');

subplot(2,2,3);  hold on; title(['Pearson P bias, for '  plot_str{plot_ind}]);
errorbar(N_VEC, Rho_bias_mean(N_VEC),Rho_bias_std(N_VEC));
xlabel('N'); ylabel('P bias'); %legend('from data', 'from Fisher');


subplot(2,2,4);  hold on; title(['Fisher Z bias, for '  plot_str{plot_ind}]);
errorbar(N_VEC, Z_bias_mean(N_VEC),Z_bias_std(N_VEC));
xlabel('N'); ylabel('Z bias'); %%legend('from data', 'from Fisher');



% Plot the same thing without the error bars
figure; subplot(2,2,1);  hold on; title(['Pearson P std, for '  plot_str{plot_ind}]);
plot(N_VEC, Rho_sig_mean(N_VEC));
xlabel('N'); ylabel('P std'); %legend('from data', 'from Fisher');


subplot(2,2,2);  hold on; title(['Fisher Z std, for ' plot_str{plot_ind}]);
plot(N_VEC, 1./sqrt(N_VEC-3), 'r');
plot(N_VEC, Z_sig_mean(N_VEC));

xlabel('N'); ylabel('Z std'); legend( '(N-3)^{-1/2}','from data');

subplot(2,2,3);  hold on; title(['Pearson P bias, for '  plot_str{plot_ind}]);
plot(N_VEC, Rho_bias_mean(N_VEC));
xlabel('N'); ylabel('P bias'); %legend('from data', 'from Fisher');


subplot(2,2,4);  hold on; title(['Fisher Z bias, for '  plot_str{plot_ind}]);
plot(N_VEC, Z_bias_mean(N_VEC));
xlabel('N'); ylabel('Z bias'); %%legend('from data', 'from Fisher');





% Plot to see correlation between the 'true' value and the variance
figure; subplot(2,2,1); hold on; plot(corrs_vec, Pearson_Rhos_std, '.'); title(['True Pearson P and estimator for '  plot_str{plot_ind} ' std. N=' num2str(N)]);
xlabel('true P'); ylabel('std');
subplot(2,2,2); hold on; plot(Fisher_corrs_vec, Fisher_Zs_std, '.'); title(['True Fisher Z and estimator  for '  plot_str{plot_ind} ' std. N='  num2str(N)]);
xlabel('true Z'); ylabel('std');


subplot(2,2,3); hold on; plot(corrs_vec, Pearson_Rhos_mean-corrs_vec', '.'); title(['True Pearson P and estimator  for ' plot_str{plot_ind} ' bias.  N=' num2str(N)]);
xlabel('true P'); ylabel('bias');
subplot(2,2,4); hold on; plot(Fisher_corrs_vec, Fisher_Zs_mean-Fisher_corrs_vec', '.'); title(['True Fisher Z and estimator  for '  plot_str{plot_ind} ' bias. N='  num2str(N)]);
xlabel('true Z'); ylabel('bias');




