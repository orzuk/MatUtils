% A script for simulating throwing random points in the unit interval 
% and finding what is the maximal gap between two consecutive points

max_n = 10000; n_win = 100;

n_vec = n_win:n_win:max_n; n_num = length(n_vec);

iters = 100;

mean_interval_mat = zeros(n_num, iters);
std_interval_mat = zeros(n_num, iters);
max_interval_mat = zeros(n_num, iters);


for i=1:iters
    rand_vec = rand(1,max_n); % draw points uniformly at random
    
%    rand_vec = sort(rand_vec); % sort them 
    ctr=1;
    for j=n_vec
        diff_vec = diff(sort(rand_vec(1:j)));
        mean_interval_mat(ctr, i) =  mean(diff_vec);
        std_interval_mat(ctr, i) =  std(diff_vec);
        max_interval_mat(ctr, i) = max(diff_vec); % sort them and take maximal interval
        ctr=ctr+1;
    end
end

figure; imagesc(max_interval_mat); colorbar;    

max_interval_mean = mean(max_interval_mat,2);
std_interval_mean = mean(std_interval_mat,2);

figure; plot(n_vec, max_interval_mean); title('max interval'); xlabel('n'); ylabel('max interval length');

figure; plot(n_vec,1./(n_vec'.* max_interval_mean)); title('inverse max interval'); xlabel('n'); ylabel('1 / max interval length');

figure; plot(n_vec, std_interval_mean); title('std interval'); xlabel('n'); ylabel('std interval length');






