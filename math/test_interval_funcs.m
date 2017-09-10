% Test all functions with intervals stuff

a1 = [ 11 2]
b1 = [  20 17]
a2 = [8 15];
b2 = [ 13 16];



a1 = [4 4 4 4 4 4 4 5 25  23 6 9 9 9 2 5 7 1 1 -13 9 9 9 9 4 -1 32 43 9 9 ];
a1 = a1 + 0.001 * randn(1, length(a1));
b1 = a1 + 1; 

a2 = [1 6 11 -13 44 4 4 11 59 39 9];
a2 = a2 + 0.001 * randn(1, length(a2));
b2 = a2 + 1;

% a1 = round(80000* rand(3, 1)); a1 = a1 +  0.01 * randn(length(a1),1); b1 = a1 + 5;
% a2 = round(60000* rand(1590000,1)) - 10; a2 = a2 + 0.01 * randn(length(a2),1); b2 = a2 + 10;


% [a1 b1 merge_inds] = intervals_merge(a1, b1); 
% [a2 b2 merge_inds] = intervals_merge(a2, b2); 

% load('bad_c_example.mat'); 
inter_mat_time = cputime; [int_start int_end int_ind1 int_ind2] = intervals_intersect(a1, b1, a2, b2,0,0); inter_mat_time = cputime - inter_mat_time
%inter_mat_time2 = cputime; [int_start int_end int_ind1 int_ind2] = intervals_intersect(a2, b2, a1, b1,0); inter_mat_time2 = cputime - inter_mat_time2
intersect_length = length(int_start);
inter_c_time = cputime; [int_start_c int_end_c int_ind1_c int_ind2_c] = intervals_intersect(a1, b1, a2, b2,1,0); inter_c_time = cputime - inter_c_time
%inter_c_time2 = cputime; [int_start_c int_end_c int_ind1_c int_ind2_c] = intervals_intersect(a2, b2, a1, b1,1); inter_c_time2 = cputime - inter_c_time2
figure; plot(int_start, int_start_c, '.')
error1_is = [max(abs(int_start - int_start_c)) max(abs(int_end - int_end_c)) ...
    max(abs(vec2column(int_ind1) - vec2column(int_ind1_c))) max(abs(vec2column(int_ind2)  - vec2column(int_ind2_c)))]
time_factor = inter_mat_time / inter_c_time


%return;

% Test merge function 
[merge_start merge_end merge_inds] = intervals_merge([a1, a2], [b1 b2])
[only1_start only1_end only1_inds] = intervals_diff(a1, b1, a2, b2);
[only2_start only2_end only2_inds] = intervals_diff(a2, b2, a1, b1);
[comp1_start comp1_end ] = intervals_complement(a1, b1, min(a2), max(b2) );
[comp2_start comp2_end ] = intervals_complement(a2, b2, min(a1), max(b1) );
figure; hold on;
Dummy = intervals_plot(a1, b1, 0, 'g', 3);
Dummy = intervals_plot(a2, b2, 1, 'r', 3);
Dummy = intervals_plot(int_start, int_end, 2, 'b', 3);
Dummy = intervals_plot(int_start_c, int_end_c, 3, 'm', 3);
Dummy = intervals_plot(merge_start, merge_end, 4, 'c', 3);
Dummy = intervals_plot(only1_start, only1_end, 5, 'k', 3);
Dummy = intervals_plot(only2_start, only2_end, 6, 'y', 3);
Dummy = intervals_plot(comp1_start, comp1_end, 7, 'g', 3);
Dummy = intervals_plot(comp2_start, comp2_end, 8, 'r', 3);

legend('first', 'second', 'intersect', 'intersect-c',...
     'merge', 'only 1', 'only 2', 'comp-1', 'comp-2'); 


a1 = [4 4 4 4 4 4 4 5 25  23 6 9 9 9 2 5 7 1 1 -13 9 9 9 9 4 -1 32 43 9 9 ];
a1 = a1 + 0.001 * randn(1, length(a1));
b1 = a1 + 1; 

a2 = [1 6 11 -13 44 4 4 11 59 39 9];
a2 = a2 + 0.001 * randn(1, length(a2));
b2 = a2 + 1;

[int_start int_end int_ind1 int_ind2] = intervals_intersect(a1, b1, a2, b2)
[int_start_c int_end_c int_ind1_c int_ind2_c] = intervals_intersect(a1, b1, a2, b2,1)
figure; plot(int_start, int_start_c, '.')




figure; 
Dummy = intervals_plot(a1, b1, 0, 'g');
Dummy = intervals_plot(a2, b2, 1, 'r');
Dummy = intervals_plot(int_start, int_end, 0.5, 'b');


figure; 
Dummy = intervals_plot(a1, b1, 0, 'g');
Dummy = intervals_plot(a2, b2, 1, 'r');
Dummy = intervals_plot(int_start, int_end, 0.5, 'b');

a1 = [0 10 30];
b1 = [5 14 40]; 

a2 = [-5 12 45];
b2 = [-3 13 46];

[closest_ind closest_dist] = match_closest_vals(a1, b1, 0)
[closest_ind_c closest_dist_c] = match_closest_vals(a1, b1)

figure; 
Dummy = intervals_plot(a1, b1, 0, 'g');
Dummy = intervals_plot(a2, b2, 1, 'r');


[closest_ind closest_dist] = match_closest_intervals(a1, b1, a2, b2)


[inter I J] = intersect_all(a1, a2) 


a = round(50*rand(100, 1)); b = round(40 * rand(50,1)) - 200; 

ttt = cputime; [closest_ind closest_dist] = match_closest_vals(a, b, 0); ttt = cputime - ttt
ttt_c = cputime; [closest_ind_c closest_dist_c] = match_closest_vals(a, b); ttt_c = cputime - ttt_c

max_ind_diff = max(abs(closest_ind - closest_ind_c))
max_dist_diff = max(abs(vec2column(closest_dist) - closest_dist_c))

figure; plot(closest_ind, closest_ind_c, '.'); title('Inds diff.');
figure; plot(closest_dist, closest_dist_c, '.'); title('Dists diff.');

% Big file: Single!! load('BAD_KMERS.mat'); 
% ttt_c = cputime; [closest_ind_c closest_dist_c] = match_closest_vals(double(X), double(Y)); ttt_c = cputime - ttt_c

