% Compute convergance radius for binary HMP when p->1/2.
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
% Here we expand around p=1/2 and want the radius of convergance their

 N = 100; iters = 520; 
use_points_for_fit = [3:82];  % Which points to use for the rational fitting
hmm_use_points_for_fit = [3:12];  % Which points to use for the hmm rational fitting
res = 5; eps0 = 0.000000001; eps_jump = 0.001; 

ttt = cputime;
eps_vec = zeros(1,res);
first_order = zeros(1,res);

max_order_taken = 5;  % this says how many orders we use in the approximation

%[p_start:0.001:p_end];
eps_vec = [0.01:0.01:0.49]; % vectors of all the p's we want to calculate
TrueRadius = p_vec ./ (1-2.*p_vec); 
FittedRadius = zeros(1, length(p_vec));
i=1;
StartPoint = [1,1,1];
HMMStartPoint = [1,1,1];

% Loop over p with different resulutions
for eps = eps_vec
    eps
    
    % First compute the orders for the i.i.d. case - Not markovian !
    %HpConv = zeros(110,1); 
    HpConv = zeros(length(p_vec)+10,1); 
    HpConv_zero = -(p*log(p)+(1-p)*log(1-p));  % different since there is no zero in a matlab array
    HpConv(1) = (1-2*p) * log ( (1-p)/p );  % Half from the markovian case
    
    for n=2:100 % 110
        HpConv(n) = (-1) * ( (2*p-1)^n / p^(n-1) +  (1-2*p)^n / (1-p)^(n-1) ) / (n*(n-1)); 
    end
    
    % First plot the orders in increasing value : 
    Hp = zeros(12, 1);  % Here's one too much !!! 
    
    % Note : We have only even powers in the expansion !!! 
    % Now represent anything as a function of lambda = 1-2p
    mu = 1-2*eps; 

    Hp_zero = -log(2.0)/mu^4;  % different since there is no zero in a matlab array
    Hp(2) = 2;       
    Hp(4) = (4/3) * (7*mu^4-12*mu^2+6);        
    Hp(6) = (32/15)*(46*mu^8-120*mu^6+120*mu^4-60*mu^2+15);        
    Hp(8) = (32/21)*(1137*mu^12-4088*mu^10+5964*mu^8-4536*mu^6+1946*mu^4-504*mu^2+84);        
    Hp(10) = (512/45)* (3346*mu^16-15120*mu^14+28800*mu^12-30120*mu^10+18990*mu^8-7560*mu^6+1980*mu^4-360*mu^2+45);
    Hp = -mu^4 .* Hp;
    
    % Generate vectors of the ratio 
    %size(Hp)
    %size([Hp_zero Hp(1:end-1)])
    HpRatio = (Hp' ./ [Hp_zero' Hp(1:end-1)'])';
    HpInverseRatio = ([Hp_zero' Hp(1:end-1)'] ./ Hp')';
    HpConvRatio = (HpConv' ./ [HpConv_zero' HpConv(1:end-1)'])';
    
    % Now try to predict from some of the points the radius of convergence 
    VecToPred = HpConvRatio(use_points_for_fit);
    VecXForPred = [0,1,2, use_points_for_fit(1:end-3)]';
    
    HMMVecToPred = HpRatio(hmm_use_points_for_fit);
    HMMInverseVecToPred = HpInverseRatio(hmm_use_points_for_fit);
    HMMVecXForPred = [0,1,2, hmm_use_points_for_fit(1:end-3)]';
    
    
    % Now try to fit a line a*x+b using least-square, and here a and b are the parameters
    [RationalFittedCurve, goodness] = fit(VecXForPred, VecToPred, 'rat11');
    [HMMRationalFittedCurve, HMMgoodness] = fit(HMMVecXForPred, HMMVecToPred, 'rat11');
 %   [HMMInverseRationalFittedCurve, HMMInversegoodness] = fit(HMMVecXForPred, HMMInverseVecToPred, 'rat11');
    
    FittedRadius(i) = abs(1 / RationalFittedCurve.p1);
    HMMFittedRadius(i) = abs(1 / HMMRationalFittedCurve.p1);
%    HMMInverseFittedRadius(i) = abs(HMMInverseRationalFittedCurve.p1);
    i = i+1;
    
    % First plot the orders for a given p
    % figure; hold on; plot([0:130], log(abs([HpConv_zero HpConv(1:130)'])), '+'); plot([0:129], HpConv(1:130)' ./ [HpConv_zero HpConv(1:129)'], 'r*');
    % title(['Ind. Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('Log', 'Cosecutives Ratio');
    % 
    % figure; hold on; plot([0:11], log(abs([Hp_zero Hp(1:11)'])), '+'); plot([0:10], Hp(1:11)' ./ [Hp_zero Hp(1:10)'], 'r*');
    % title(['HMM Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('Log', 'Cosecutives Ratio');
    
    % figure; hold on; plot([0:101], log(abs([HpConv_zero HpConv(1:101)'])), '*r'); %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
    % title(['IID epsilon Log Abs Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); %legend('I.I.D.', 'HMM');
    
    
%     % Plot the fit for the HMM : 
%     HMMVecXForPredCurve = [-0.5:0.01:9.5]; 
%     HMMVecToPredCurve = (HMMRationalFittedCurve.p1 .* HMMVecXForPredCurve + HMMRationalFittedCurve.p2) ./ (HMMVecXForPredCurve + HMMRationalFittedCurve.q1);
%     figure; hold on; plot(HMMVecXForPred, HMMVecToPred, '.');  plot(HMMVecXForPredCurve, HMMVecToPredCurve, 'r'); legend('Points', 'Fitted Curve'); 
%     title(['Orders ratio and Fitted Curve for P = ' num2str(p)]); xlabel('K-2'); ylabel('Ratio'); 
    
    
end % for loop on the p's
figure; hold on; plot(p_vec, TrueRadius); plot(p_vec, FittedRadius, 'r'); plot(p_vec, HMMFittedRadius, 'm'); 
legend('True Radius', 'Predicted Radius', 'HMM Predicted Radius'); 
title(['Radius of Convergance for i.i.d. and HMM, predicting using orders ' num2str(use_points_for_fit(1)) ' - ' num2str(use_points_for_fit(end))] ); xlabel('P'); ylabel('Radius');

figure; hold on; plot(p_vec, (FittedRadius-TrueRadius) ./ TrueRadius);   
title(['Relative Radius prediction error for i.i.d., predicting using orders ' num2str(use_points_for_fit(1)) ' - ' num2str(use_points_for_fit(end))]); xlabel('P'); ylabel('(Predicted Radius - Radius)/Radius');


figure; hold on;  plot(eps_vec, HMMFittedRadius, 'b');legend('HMM Predicted Radius and Inverse '); 
title('Predicted Radius of Convergance for HMM'); xlabel('\epsilon'); ylabel('HMM Pred. Radius'); 

%%%legend('Straight', 'Inverse'); plot(p_vec, HMMInverseFittedRadius, 'r'); 

max_plot = 30; 

figure; hold on; plot([0:max_plot], HpConvRatio(1:max_plot+1), 'r*'); %%%%%   HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'], 'r*');  
plot([0:max_plot], ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] ) ./ ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )  , 'om'); 

% % % fraction = ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] ) ./ ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )
% % % num = ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] )
% % % denum = ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )

plot([0:max_plot], RationalFittedCurve([0:max_plot]), 'k');
% plot([0:100], (-[3:103]./[1:101]) .* ((1-2*p)/p), 'm'); 
title(['epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

% figure; hold on; plot(log([1:101]), log(abs(HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'] - ((1-2*p)/p) )), 'r*'); 
% title(['epsilon Log Ratio of Orders for p = ' num2str(p)]); xlabel('Log Order'); ylabel('Log Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');
% 
% figure; hold on; plot([0:100], log(abs(HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'] - ((1-2*p)/p) )), 'r*'); 
% title(['epsilon Log Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

% figure; hold on; plot([0:11], log(abs([HpConv_zero HpConv(1:11)'])), '*r'); plot([0:11], log(abs([Hp_zero Hp(1:11)'])), '+');  
% title(['epsilon Log Abs Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('I.I.D.', 'HMM');
% 
% figure; hold on; plot([0:10], HpConv(1:11)' ./ [HpConv_zero HpConv(1:10)'], 'r*'); plot([0:10], Hp(1:11)' ./ [Hp_zero Hp(1:10)'], '+'); 
% title(['epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); legend('I.I.D.', 'HMM');
% 
% figure; hold on; plot([0:10], log(HpConv(1:11)' ./ [HpConv_zero HpConv(1:10)']), 'r*'); plot([0:10], log(Hp(1:11)' ./ [Hp_zero Hp(1:10)']), '+'); 
% title(['Log epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); legend('I.I.D.', 'HMM');

total_time = cputime - ttt
