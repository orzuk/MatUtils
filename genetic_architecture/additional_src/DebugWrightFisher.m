% Debug simulation of wright-fisher model
%function DebugWrightFisher()

load('debug_SFS.mat');  x_vec = freq_struct.x_vec{end-1}; p_vec = freq_struct.p_vec{end-1}; q = simulation_struct.q;  % load data without(?) conditioning on polymorphic
load('debug_SFS2.mat'); x_vec2 = freq_struct2.x_vec{end-1}; p_vec2 = freq_struct2.p_vec{end-1}; q2 = simulation_struct2.q; % load data with(?) conditioning on polymorphic
age_vec = allele_age_from_sim_struct(simulation_struct);
age_vec2 = allele_age_from_sim_struct(simulation_struct2);


figure; % Compare DAF distribution at different generations
subplot(2,2,1); tt=100; semilogy(age_vec, simulation_struct.q(:,tt), '*'); hold on;
semilogy(age_vec2, simulation_struct2.q(:,tt), 'r.'); xlabel('age'); ylabel('DAF'); legend({'All', 'Poly'}); title(['t=' num2str(tt)]);
I=find(q(:,1)); I2 = find(q2(:,1)); med_vec = zeros(size(q,2), 1); med_vec2 = zeros(size(q,2), 1); 
for i=1:size(q,2) % compute median DAF
    med_vec(i) = median(q(q(:,i)>0,i));
    med_vec2(i) = median(q2(q2(:,i)>0,i));
end
subplot(2,2,2);  semilogy(med_vec, 'b'); hold on; semilogy(med_vec2, 'r'); xlabel('t'); ylabel('DAF'); title('Median DAF of polymorphic alleles'); legend({'All', 'Poly'});
subplot(2,2,3);  semilogy(max(q), 'b'); hold on; semilogy(max(q2), 'r');  xlabel('t'); ylabel('DAF'); title('MAX DAF of polymorphic alleles'); legend({'All', 'Poly'});

% Plot what??
subplot(2,2,4);  hold on;
for i=11:30
    plot(q(I(i),:), 'b'); hold on;
    plot(q2(I(i),:), 'r');
end
legend({'All', 'Poly'}); xlabel('t'); ylabel('DAF'); title('Trajectories of random alleles');


%%%%%%%%%%%%% NEW!! Calculate everything from scratch !!! %%%%%%%%%%%%%
Demographic_model{i_pop}.iters = 10000; % number of alleles to simulate !!
Demographic_model{i_pop}.s_grid = [0 -logspace(-6, -2, 101)]; % s vector for interpolation
compute_flag = []; compute_flag.method = 'simulation'; compute_flag.smooth = 1;
Demographic_model{i_pop}.save_flag = 1; Demographic_model{i_pop}.cond_on_polymorphic_flag=0
[Demographic_model{i_pop}.SFS.x_vec, Demographic_model{i_pop}.SFS.p_vec, ...
    Demographic_model{i_pop}.SFS.L, SFS_compute_time] = ...
    compute_allele_freq_spectrum_from_demographic_model( ...
    Demographic_model{i_pop}, s_vec, compute_flag); %  s_vec([1 5 10])
save(demography_file, 'Demographic_model', 'max_LL_demographic_model');


