% Run the max-min program and count extremal sequences 

res = 0.01; % resulution of grid
N=3; % chain length
Y = []; % currently
p_vec = res/2-0.5:res:0.5-res/2; eps_vec = p_vec; % we take only from -1/2 to 1/2
seqs_flag = 0;

[min_seqs, min_vals, max_seqs, max_vals] = FindAllMaximalYSeqs(eps_vec, p_vec, N, seqs_flag, Y);

figure; subplot(2,2,1); imagesc(min_seqs); colorbar; title('min sequences'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,2); imagesc(sign(min_vals)); colorbar; title('min vals'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,3); imagesc(max_seqs); colorbar; title('max sequences'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,4); imagesc(sign(max_vals)); colorbar; title('max vals'); xlabel('\epsilon'); ylabel('P');

p = 0.2;  eps_real_vec = eps_vec; eps_imag_vec = eps_vec;
[min_real_seqs min_real_vals min_imag_seqs min_imag_vals  ...
    max_real_seqs max_real_vals max_imag_seqs max_imag_vals ] =  PlotComplexProbY(eps_real_vec, eps_imag_vec, p, N, seqs_flag, Y);



figure; subplot(2,2,1); imagesc(min_real_seqs); colorbar; title('min (real) sequences'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,2); imagesc(sign(min_real_vals)); colorbar; title('min real vals'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,3); imagesc(min_imag_seqs); colorbar; title('min (imag) sequences'); xlabel('\epsilon'); ylabel('P');
subplot(2,2,4); imagesc(sign(min_imag_vals)); colorbar; title('min imag vals'); xlabel('\epsilon'); ylabel('P');



% Find for each Y the best explanation X
N=7;
seqs_flag = 1; % 1 - simulate stuff (MLE) 3 - Newton Polytopes %  P(X|Y)
MLE = 0; MARGINAL = 1; inf_flag = MARGINAL;
res = 0.005;
p_vec = 0.5*res:res:1-0.5*res; 
eps_res = 0.005; eps_vec = 0.5*eps_res:eps_res:1-0.5*eps_res; ; % we take only from -1/2 to 1/2
n = length(p_vec); m = length(eps_vec);
min_seqs = zeros(2^N, m, n); max_seqs = zeros(2^N, m, n);
figure;hold on;
NF = cell(2^N,1); 
for Y=0:2^N-1
    [min_seqs(Y+1,:,:) , min_vals, max_seqs(Y+1,:,:), max_vals, NF{Y+1}] = HMPInferenceFunction(eps_vec, p_vec, N, seqs_flag, Y, inf_flag);
end
figure; hold on; % Plot empirical inf. funcs. in [p,epsilon] plane
for Y=0:2^N-1
    subplot(2^floor(N/2),2^ceil(N/2),Y+1); imagesc(reshape(max_seqs(Y+1,:,:),m,n)); colorbar; title('Maximal X sequence');
end
[InfFuncs num_funcs] = Refinement(max_seqs);

All_NF = Refinement(NF);
figure; hold on; % figure;
for i=1:size(All_NF,1)
    line(  [0 All_NF(i,1)], [0,All_NF(i,2)], 'linewidth', 3); % Plot the Newton Polytope
end
circle([0,0], 1, 500, '.');
%            axis( [min(NP(:,1))-1 max(NP(:,1))+1 min(NP(:,2))-1 max(NP(:,2))+1]);
title(['Normal Fan Refinement for all Ys has' num2str(length(NF)) ' different cones']);


figure; hold on; % Plot empirical inf. funcs. after log-transform in [p,epsilon] plane
x_vec = repmat(p_vec, 1, m); y_vec = reshape(repmat(eps_vec, n, 1), n*m,1)';
for Y=0:2^N-1
    U = reshape(max_seqs(Y+1,:,:)',  m,n);
    subplot(2^floor(N/2),2^ceil(N/2),Y+1); scatter(x_vec, y_vec, U); %imagesc(reshape(max_seqs(Y+1,:,:),m,n)); colorbar; title('Maximal X sequence');
end




for i=1:m
    InfFuncs(i,i) = 1; 
end
unique(InfFuncs)





