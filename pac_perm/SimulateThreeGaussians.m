% A script for computing overlaps for three gaussians 
mega_iters = 100;
Fisher_Z_dev_corr = zeros(1,mega_iters);
Expression_corr_g1_s = zeros(1,mega_iters);
Expression_corr_g2_s = zeros(1,mega_iters);
Expression_corr_g1_g2 = zeros(1,mega_iters);

for mega_iter=1:mega_iters
% The correlation matrix - Only positive entries, so hopefully it is
% positive definite
C = rand(3); C = 0.5*(C+C'); C(1,1) = 1; C(2,2) = 1; C(3,3) = 1; A  = sqrtm(C);

nsamples = 100;

% Now sample
Data  = randn(3,nsamples);

% Transform to get the corrleation structure 
Data = A*Data; 

% Original data correlations
Datacorrs = zeros(3);
for i=1:3
    for j=1:3
        Datacorrs(i,j) = corr(Data(i,:)',Data(j,:)');
    end
end


% Now do the Fisher staff
block_size = 20;
iters = 1000;



 block_corr_one  =zeros(3);  block_corr_two  =zeros(3);
block_corr_one = zeros(1,iters); block_corr_two = zeros(1,iters);
 
 for i=1:iters
    randp = randperm(nsamples);
    
    block_one = Data(:,randp(1:block_size));
    %block_two = Data(:,randp(1:block_size));
    
    block_corr_one(i) = corr(block_one(1,:)', block_one(3,:)');
    block_corr_two(i) = corr(block_one(2,:)', block_one(3,:)');
    
end % iters loop 


% Do fisher transform and calc corr
% Fisher_Z_dev_corr(mega_iter) = corr(atanh(block_corr_one)'-atanh(C(1,3)), atanh(block_corr_two)'-atanh(C(2,3)));
Fisher_Z_dev_corr(mega_iter) = corr((block_corr_one)'-(C(1,3)), (block_corr_two)'-(C(2,3)));
Expression_corr_g1_s(mega_iter) = C(1,3);
Expression_corr_g2_s(mega_iter) = C(2,3);
Expression_corr_g1_g2(mega_iter) = C(1,2);




end




% Plot results
figure; imagesc(Datacorrs); colorbar; title('Corrs in data');

figure; plot(C(:), Datacorrs(:), '.'); title('true corrs vs. corrs in data'); xlabel('true'); ylabel('data');



figure; hold on; plot(Expression_corr_g1_g2, Fisher_Z_dev_corr, '.'); title('Corr in expresion vs. corr in Fisher Z'); xlabel('exp. corr.'); ylabel('Z corr.');


figure; hold on; plot(Expression_corr_g1_s, Fisher_Z_dev_corr, '.'); title('Corr in expresion g1 vs. corr in Fisher Z'); xlabel('exp. corr.'); ylabel('Z corr.');

figure; hold on; plot(Expression_corr_g2_s, Fisher_Z_dev_corr, '.'); title('Corr in expresion g2 vs. corr in Fisher Z'); xlabel('exp. corr.'); ylabel('Z corr.');


figure; hold on; plot3(Expression_corr_g1_s, Expression_corr_g2_s, Fisher_Z_dev_corr); 
title('Corr in expresion both genes  vs. corr in Fisher Z'); xlabel('exp. corr.'); ylabel('Z corr.');




