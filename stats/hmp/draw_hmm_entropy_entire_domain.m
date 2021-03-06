% Plot the entropy of HMP for various parameters in domain

p = 0.3;
eps = 0.1;

segs = 1000;

tolerance = 0.000000000001;
mul_ent = zeros(segs-1,segs-1);
upper_ent = zeros(segs-1,segs-1);



% New ! We here do not let p and epsilon be 0, 0.5
p_vec = [0:1/(segs):1]';
eps_vec = [-1:3/(segs):2]';
% p_vec = [0:0.005/(segs):0.005]';
% eps_vec = [0:0.005/(segs):0.005]';
p_vec = p_vec(2:end-1); 
eps_vec = eps_vec(2:end-1); 


% Now we don't need the tolerance 
% eps_vec = max(eps_vec, tolerance); eps_vec = min(eps_vec, 0.5-tolerance); 
% p_vec = max(p_vec, tolerance); p_vec = min(p_vec, 0.5-tolerance); 


% Try the new bounds and the 2nd order approximation
Ido_upper = zeros(segs-1,segs-1);
Ido_lower = zeros(segs-1,segs-1);

% The new general upper/lower bounds
C_N = zeros(segs-1,segs-1);
cx_N = zeros(segs-1,segs-1);
C_N_Plus1 = zeros(segs-1,segs-1);
cx_N_Plus1 = zeros(segs-1,segs-1);

% Do small diriclet correction to avoid zero probabilities.
dirich = 0.000000000000001;
p_vec = (p_vec + dirich) / (1+2*dirich);

M=7;

N=5; % Set the number of bits to use in the upper/lower bounds
% Now loop over eps


for i=1:segs-1  % loop over p
    if(mod(i, 100) == 0)
        i
    end
            
    T_N = HMP_entropy_finite(eps_vec(i), p_vec', N);
    tx_N = HMP_entropy_finite_X1(eps_vec(i), p_vec', N);
    C_N(:,i) = T_N-HMP_entropy_finite(eps_vec(i), p_vec', N-1); % The upper bound
    cx_N(:,i) = tx_N-HMP_entropy_finite_X1(eps_vec(i), p_vec', N-1); % The lower bound
%     C_N_Plus1(:,i) = HMP_entropy_finite(eps_vec(i), p_vec', N+1)-T_N; % The upper bound N+1
%     cx_N_Plus1(:,i) = HMP_entropy_finite_X1(eps_vec(i), p_vec', N+1)-tx_N; % The lower bound N+1
end

S_N = real(C_N);
sx_N = real(cx_N);
C_N = real(C_N .* (C_N == abs(C_N)) );
cx_N = real(cx_N .* (abs(imag(cx_N)) < tolerance));


% Distinguish real zero from close to zero
C_N ( C_N == 0) = -0; 
cx_N ( cx_N == 0) = -0; 


figure; subplot(2,2,1); hold on; 
subplot(2,2,1); hold on; imagesc(eps_vec, p_vec, (C_N+cx_N)./2); colorbar('vert'); title(['Upperbound Plus Lowerbound Average C_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(eps_vec, p_vec, max(C_N-cx_N, 0)); colorbar('vert'); title(['Upperbound Minus Lowerbound C_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(eps_vec, p_vec, max(2.*(C_N-cx_N)./(C_N+cx_N), 0)); colorbar('vert'); title('Relative Error Upperbound Minus Lowerbound / Average C_N' ); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(eps_vec, p_vec, min(max((C_N-cx_N)./(1 - (C_N+cx_N)./2), 0), 3)); colorbar('vert'); title('Relative Error Upperbound Minus Lowerbound / (1-Average) C_N'); ylabel('p'); xlabel('eps');


figure; subplot(2,2,1); hold on; imagesc(eps_vec, p_vec, C_N); colorbar('vert'); title(['Upperbound C_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(eps_vec, p_vec, cx_N); colorbar('vert'); title(['Lowerbound cx_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(eps_vec, p_vec, S_N); colorbar('vert'); title(['Upperbound S_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(eps_vec, p_vec, sx_N); colorbar('vert'); title(['Lowerbound sx_N N=' num2str(N)]); ylabel('p'); xlabel('eps');

% % Now try to plot how the upper and lower bounds approach the entropy 
% figure; subplot(2,2,1); hold on;
% 
% subplot(2,2,1); hold on; imagesc(eps_vec, p_vec, (C_N-C_N_Plus1)); colorbar('vert'); title(['Upperbound N+1 Minus Upper C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% subplot(2,2,2); hold on; imagesc(eps_vec, p_vec, (cx_N-cx_N_Plus1)); colorbar('vert'); title(['Lowerbound N+1 Minus Lower C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% subplot(2,2,3); hold on; imagesc(eps_vec, p_vec, 2.*(C_N-C_N_Plus1)./(C_N+C_N_Plus1)); colorbar('vert'); title(['Relative Error (Upperbound N+1 Minus Upper)/Average C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% subplot(2,2,4); hold on; imagesc(eps_vec, p_vec, 2.*(cx_N-cx_N_Plus1)./(cx_N+cx_N_Plus1)); colorbar('vert'); title(['Relative Error (Lowerbound N+1 Minus Lower)/Average C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% 
