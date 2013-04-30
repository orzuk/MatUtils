% Predict convergance radius using fitting
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary

 N = 100; iters = 520; 
use_points_for_fit = [3:11];  % Which points to use for the rational fitting
hmm_use_points_for_fit = [3:11];  % Which points to use for the hmm rational fitting
res = 5; eps0 = 0.000000001; eps_jump = 0.001; 

ttt = cputime;
eps_vec = zeros(1,res);
first_order = zeros(1,res);

max_order_taken = 5;  % this says how many orders we use in the approximation


p_vec = [0.05:0.0025:0.35]; % vectors of all the p's we want to calculate
TrueRadius = p_vec ./ (1-2.*p_vec); 
FittedRadius = zeros(1, length(p_vec));
i=1;
StartPoint = [1,1,1];
HMMStartPoint = [1,1,1];


p1_vec = zeros(1, length(p_vec)); 
p2_vec = zeros(1, length(p_vec)); 
q1_vec = zeros(1, length(p_vec)); 


HMMp1_vec = zeros(1, length(p_vec)); 
HMMp2_vec = zeros(1, length(p_vec)); 
HMMq1_vec = zeros(1, length(p_vec)); 


% Loop over p with different resulutions
for p = p_vec
    
    p
    % First compute the orders for the i.i.d. case - Not markovian !
    %HpConv = zeros(110,1); 
    HpConv = zeros(length(p_vec)+10,1); 
    HpConv_zero = -(p*log(p)+(1-p)*log(1-p));  % different since there is no zero in a matlab array
    HpConv(1) = (1-2*p) * log ( (1-p)/p );  % Half from the markovian case
    
    for n=2:100 % 110
        HpConv(n) = (-1) * ( (2*p-1)^n / p^(n-1) +  (1-2*p)^n / (1-p)^(n-1) ) / (n*(n-1)); 
    end
    
    % First plot the orders in increasing value : 
    Hp = zeros(11, 1);  % Here's one too much !!! 
    
    Hp_zero = -(p*log(p)+(1-p)*log(1-p));  % different since there is no zero in a matlab array
    Hp(1) = 2 *(1-2*p) * log ( (1-p)/p );
    Hp(2) = -2 *(1-2*p) * log ( (1-p)/p ) -(1-2*p)^2/(2*p^2*(1-p)^2);
    
    % Now represent anything as a function of lambda = 1-2p
    lam = 1-2*p; 
    
    Hp(3) = (-16*(5*lam^4-10*lam^2-3) * lam^2) / ( 3*(1-lam^2)^4);
    
    Hp(4) = (8*(109*lam^8+20*lam^6-114*lam^4-140*lam^2-3)*lam^2) / (3*(1-lam^2)^6);
    
    Hp(5) = (-128*(95*lam^10+336*lam^8+762*lam^6-708*lam^4-769*lam^2-100)*lam^4) / (15*(1-lam^2)^8);
    
    Hp(6) = 128*(125*lam^14-321*lam^12+9525*lam^10+16511*lam^8-7825*lam^6- ...
        17995*lam^4-4001*lam^2-115)*lam^4/   (15*(1-lam^2)^10);
    
    Hp(7) = -256*(280*lam^18-45941*lam^16-110888*lam^14+666580*lam^12+1628568*lam^10- ...
        270014*lam^8-1470296*lam^6-524588*lam^4-37296*lam^2-245)*lam^4 / (105*(1-lam^2)^12);
    
    Hp(8) = 64*(56*lam^22-169169*lam^20-2072958*lam^18-5222301*lam^16+12116328*lam^14+ ...
        35666574*lam^12+3658284*lam^10-29072946*lam^8-14556080*lam^6- ...
        1872317*lam^4-48286*lam^2-49)*lam^4 / (21*(1-lam^2)^14);
    
    Hp(9) = 2048*(37527*lam^22+968829*lam^20+8819501*lam^18+20135431*lam^16-23482698*lam^14- ...
        97554574*lam^12-30319318*lam^10+67137630*lam^8+46641379*lam^6+8950625*lam^4+ ...
        495993*lam^2+4683)*lam^6 / (63*(1-lam^2)^16);
    
    Hp(10) = -2048*(38757*lam^26+1394199*lam^24+31894966*lam^22+243826482*lam^20+ ...
        571835031*lam^18-326987427*lam^16-2068579420*lam^14-1054659252*lam^12+1173787011*lam^10+ ...
        1120170657*lam^8+296483526*lam^6+26886370*lam^4+ 684129*lam^2+2187)*lam^6/  (45*(1-lam^2)^18);
    
    Hp(11) = 8192*(98142*lam^30-1899975*lam^28+92425520*lam^26+3095961215*lam^24+ ...
        25070557898*lam^22+59810870313*lam^20-11635283900*lam^18-173686662185*lam^16- ...
        120533821070*lam^14+74948247123*lam^12+102982107048*lam^10+35567469125*lam^8+ ...
        4673872550*lam^6+217466315*lam^4+2569380*lam^2+2277)*lam^6/ (495*(1-lam^2)^20);
    
    
    % Generate vectors of the ratio 
    %size(Hp)
    %size([Hp_zero Hp(1:end-1)])
    HpRatio = (Hp' ./ [Hp_zero' Hp(1:end-1)'])';
    HpInverseRatio = ([Hp_zero' Hp(1:end-1)'] ./ Hp')';
    HpConvRatio = (HpConv' ./ [HpConv_zero' HpConv(1:end-1)'])';
    
    % Now try to predict from some of the points the radius of convergence 
    VecToPred = HpConvRatio(use_points_for_fit);
    VecXForPred = [2, use_points_for_fit(1:end-1)]';
    
    HMMVecToPred = HpRatio(hmm_use_points_for_fit);
    HMMInverseVecToPred = HpInverseRatio(hmm_use_points_for_fit);
    HMMVecXForPred = [2, hmm_use_points_for_fit(1:end-1)]';
    
    
    % Try New Fit , for IID !!!! 
%    ggg = fittype('(p1*x+p1*p2)/(x+p2-2)');
%    [RationalFittedCurve, goodness] = fit(VecXForPred, VecToPred, ggg);
%    [HMMRationalFittedCurve, HMMgoodness] = fit(HMMVecXForPred, HMMVecToPred, ggg);
       
       
    % Now try to fit a line a*x+b using least-square, and here a and b are the parameters
    [RationalFittedCurve, goodness] = fit(VecXForPred, VecToPred, 'rat11');
    [HMMRationalFittedCurve, HMMgoodness] = fit(HMMVecXForPred, HMMVecToPred, 'rat11');
 %   [HMMInverseRationalFittedCurve, HMMInversegoodness] = fit(HMMVecXForPred, HMMInverseVecToPred, 'rat11');
    
 
    % Save the fits : 
    p1_vec(i) = RationalFittedCurve.p1;
    p2_vec(i) = RationalFittedCurve.p2;
  %  q1_vec(i) = RationalFittedCurve.q1;
    
    HMMp1_vec(i) = HMMRationalFittedCurve.p1;
    HMMp2_vec(i) = HMMRationalFittedCurve.p2;
%    HMMq1_vec(i) = HMMRationalFittedCurve.q1;

    
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
figure; hold on; plot(p_vec, TrueRadius, 'k'); plot(p_vec, FittedRadius, '-.k'); plot(p_vec, HMMFittedRadius, '--k'); 
legend('True Radius', 'Predicted Radius', 'HMM Predicted Radius'); 
title(['Radius of Convergance for i.i.d. and HMM, predicting using orders ' num2str(use_points_for_fit(1)) ' - ' num2str(use_points_for_fit(end))] ); xlabel('P'); ylabel('Radius of Convergence');

%%%figure; hold on; plot(p_vec, (FittedRadius-TrueRadius) ./ TrueRadius);   
%%%title(['Relative Radius prediction error for i.i.d., predicting using orders ' num2str(use_points_for_fit(1)) ' - ' num2str(use_points_for_fit(end))]); xlabel('P'); ylabel('(Predicted Radius - Radius)/Radius');


figure; hold on;  plot(p_vec, HMMFittedRadius, 'b');legend('HMM Predicted Radius and Inverse '); 
title('Predicted Radius of Convergance for HMM'); xlabel('P'); ylabel('HMM Pred. Radius'); 

%%%legend('Straight', 'Inverse'); plot(p_vec, HMMInverseFittedRadius, 'r'); 

% % % % % max_plot = 30; 
% % % % % 
% % % % % figure; hold on; plot([0:max_plot], HpConvRatio(1:max_plot+1), 'r*'); %%%%%   HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'], 'r*');  
% % % % % plot([0:max_plot], ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] ) ./ ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )  , 'om'); 
% % % % % 
% % % % % % % % fraction = ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] ) ./ ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )
% % % % % % % % num = ([-1:max_plot-1]./  [1:max_plot+1]) .* ((2*p-1)/p) .*    ( 1 - (p/(p-1)) .^[0:max_plot] )
% % % % % % % % denum = ( 1 - (p/(p-1)) .^[-1:max_plot-1]  )
% % % % % 
% % % % % plot([0:max_plot], RationalFittedCurve([0:max_plot]), 'k');
% % % % % % plot([0:100], (-[3:103]./[1:101]) .* ((1-2*p)/p), 'm'); 
% % % % % title(['epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

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
