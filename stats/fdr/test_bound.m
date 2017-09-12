% Test various FDR bounds 
clear all
close all
clc
%solid red line p value vector
%solid color line exact FDR
%dash color line my FDR
%dash dot color line Benjamini FDR
tic
color='bgrmkcy';
figure;
Q=zeros(1,10000);
F=Q;
for w=1:1
    x=[normrnd(0,1,2000,1);normrnd(2,1,3000,1);normrnd(-1,1,5000,1)];
    P=2*normcdf(-abs(x),0,1);
    [Psort,XI]=sort(P);
    plot([1:length(P)]/length(P),Psort,'r');grid on;hold on;
    out=zeros(10000,1);
    out(find(XI<=2000))=1; % false positive (H0)
    Q = cumsum(out)' ./ [1:length(P)];
%    for i=1:length(P)
%        Q(i)=sum(out(1:i))/i;%exact FDR  
%    end
    plot([1:length(P)]/length(P),Q,color(w-7*floor((w-1)/7)));grid on;hold on; % plot exact FDR
    for i=1:length(P)
        F(i)=(length(P)/i)*2*sum(Psort(1:i))/i;%my FDR upper bound
    end
    plot([1:length(P)]/length(P),F,['--',color(w-7*floor((w-1)/7))]);grid on;hold on; % plot Amit's FDR
    q=Psort./[1:length(Psort)]'.*length(Psort);%Benjamini FDR bound
    plot([1:length(P)]/length(P),q,['-.',color(w-7*floor((w-1)/7))]);grid on;hold on; % plot BH FDR
end
legend('pvals', 'true FDR', 'Amit FDR', 'BH FDR');
toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test using the new functions 
iters = 5000; m=1000; m0 = 500; mu_gap=2; dep_flag = 0; rho=[0.3]; fdr_str = 'BH95'; q = 0.5;
[VR_hist V_moments R_moments FDR_moments] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,fdr_str,q);

VR = cell(m,1); V_m = VR; R_m = VR; FDR_m = VR; 
[VR{1} V_m{1} R_m{1} FDR_m{1}] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'min_k',0.36)
[VR{2} V_m{2} R_m{2} FDR_m{2}] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'min_k',0.68)
[VR{3} V_m{3} R_m{3} FDR_m{3}] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'min_k',1)
% [VR{4} V_m{4} R_m{4} FDR_m{4}] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'min_k',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=100; m0 = 0.5*m; mu_gap = 2; dep_flag=0;  q=0.05; oracle_q = q*m/m0; iters = 10000;
rho_vec = 0:0.01:1; n = length(rho_vec);
V_moments = zeros(2,n); 
R_moments = zeros(2,n); 
FDR_moments = zeros(2,n); 
for i=1:n    
    i_is = i
    [VR_hist V_moments(:,i) R_moments(:,i) FDR_moments(:,i)] = ...
        simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho_vec(i),'bh95',oracle_q);
end
figure; hold on; plot(rho_vec, FDR_moments(1,:), '.'); plot(rho_vec, FDR_moments(2,:), 'r.'); 
legend('mean', 'std'); xlabel('\rho');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Amit's problematic oracle example (where FDR was more than theoretic prediction)
m = 500; q=0.05; m0 = 0.2*m; mu_gap = 3.5; rho = 0.8; dep_flag = 0; iters = 10000; oracle_q = q*m/m0; two_side = 1; % bound works for one-sided test
[oracle_VR oracle_V_m oracle_R_m oracle_FDR_m] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'bh95',oracle_q ,two_side);
% figure; imagesc(oracle_VR); colorbar;title('2d Histogram of V vs. R'); xlabel('V'); ylabel('R'); 

FDR_control_should_be_negative = oracle_FDR_m(1) - q % check if we indeed get control (negative) pr violation (positive)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Try various distributions on two p-values and see if we can get
%%% violation of strong control with q (Shlomo's distribution)
m = 2; q=0.5; m0 = m; mu_gap = 3.5; rho = q/2; dep_flag = 2; iters = 100000; oracle_q = q*m/m0; two_side = 1; % bound works for one-sided test
[oracle_VR oracle_V_m oracle_R_m oracle_FDR_m] = simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,'bh95',oracle_q ,two_side);
FDR_control_should_be_negative = oracle_FDR_m(1) - q % check if we indeed get control (negative) pr violation (positive)
% figure; imagesc(oracle_VR); colorbar;title('2d Histogram of V vs. R'); xlabel('V'); ylabel('R'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the new function for computing q from the FDR 
pvals_vec = simulate_pvals(1,m,m0,mu_gap,dep_flag,rho,[],1); % simulate p-vals (sorted)
q_vec = [0.01:0.01:1]; n = length(q_vec); q_est = zeros(n,1);
for i=1:n % loop over q
   R = FDR_mat_main(pvals_vec, q_vec(i), 'bh95'); 
   q_est(i) =  q_from_R_FDR_mat_main(pvals_vec, R, 'bh95', 1);
end    

figure; hold on; plot(q_vec, q_est, '.'); 
plot([0:0.1:1], [0:0.1:1], 'r');
title('true q vs. inverse q'); xlabel('q'); ylabel('q estimated');

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New functions: for site: 
  P = rand(1,1000); % generate a vector of p-values, of which a 100 are taken from the alternative hypothesis
  P(1:500) = P(1:500) ./ 3; 
  q = FDR_R_to_q(sort(P), 150)  % Compute the FDR of the lowest 100 p-values








