% Compute finite size relative entropy for hmm;


segs = 500;

tolerance = 0.000000000001;
mul_ent = zeros(segs-1,segs-1);
upper_ent = zeros(segs-1,segs-1);




% New ! We here do not let p and epsilon be 0, 0.5
q_vec = [0:0.5/(segs):0.5]';
delta_vec = [0:0.5/(segs):0.5]';
% q_vec = [0:0.005/(segs):0.005]';
% delta_vec = [0:0.005/(segs):0.005]';
q_vec = q_vec(2:end-1); 
delta_vec = delta_vec(2:end-1); 


% New : Here we choose the model we want to compare to : 
eps = 0.1; p = 0.4; p_vec = zeros(1,length(q_vec))+p;
conv_c = eps*p + (1-eps)*(1-p);



% Now we don't need the tolerance 
% delta_vec = max(delta_vec, tolerance); delta_vec = min(delta_vec, 0.5-tolerance); 
% q_vec = max(q_vec, tolerance); q_vec = min(q_vec, 0.5-tolerance); 


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
q_vec = (q_vec + dirich) / (1+2*dirich);

M=7;




N=4; % Set the number of bits to use in the upper/lower bounds
% Now loop over eps




% Here we do KL !!!!!
for i=1:segs-1  % loop over p
    if(mod(i, 100) == 0)
        i
    end
            
    T_N = HMP_KL_finite(eps,p_vec, delta_vec(i), q_vec', N);
    tx_N = HMP_KL_finite_X1(eps,p_vec, delta_vec(i), q_vec', N);
    C_N(:,i) = T_N-HMP_KL_finite(eps,p_vec, delta_vec(i), q_vec', N-1); % The upper bound
    cx_N(:,i) = tx_N-HMP_KL_finite_X1(eps,p_vec, delta_vec(i), q_vec', N-1); % The lower bound
    C_N_Plus1(:,i) = HMP_KL_finite(eps,p_vec, delta_vec(i), q_vec', N+1)-T_N; % The upper bound N+1
    cx_N_Plus1(:,i) = HMP_KL_finite_X1(eps,p_vec, delta_vec(i), q_vec', N+1)-tx_N; % The lower bound N+1
end


% Do transpose to change p/eps
% C_N = C_N';  cx_N = cx_N'; 

% Transfer everything to binary entropy !!! 
% Ido_upper = Ido_upper ./ log(2.0); 
% Ido_lower = Ido_lower ./ log(2.0); 
% upper_ent = upper_ent  ./ log(2.0);
% mul_ent = mul_ent  ./ log(2.0);

% figure; subplot(3,3,1); hold on; 
% subplot(3,3,2); hold on; imagesc((Ido_upper)); colorbar('vert'); title('HMM Ido 2-step Upperbound'); ylabel('p'); xlabel('eps');
% subplot(3,3,3); hold on; imagesc((Ido_lower)); colorbar('vert'); title('HMM Ido 2-step Lowerbound'); ylabel('p'); xlabel('eps');
% subplot(3,3,5); hold on; imagesc((Ido_upper-Ido_lower)); colorbar('vert'); title('Upperbound Minus Lowerbound'); ylabel('p'); xlabel('eps');
% subplot(3,3,8); hold on; imagesc(2*(Ido_upper-Ido_lower)./(Ido_upper+Ido_lower)); colorbar('vert'); title('Upperbound Minus Lowerbound Relative Error'); ylabel('p'); xlabel('eps');

conv_eps_vec = (q_vec + conv_c-1) ./ (2.*q_vec - 1);
conv_eps_vec = conv_eps_vec(find(conv_eps_vec > 0));


figure; subplot(2,2,1); hold on; 
subplot(2,2,1); hold on; imagesc(delta_vec, q_vec, max(log(max(C_N,tolerance)),-15)); colorbar('vert'); title(['KL log  Upperbound C_N N=' num2str(N) ' (p,\epsilon) = (' num2str(p) ' , ' num2str(eps), ')']); ylabel('q'); xlabel('delta');
plot(q_vec(1:length(conv_eps_vec)), conv_eps_vec, 'm');


subplot(2,2,3); hold on; imagesc(delta_vec, q_vec, max(log(max(cx_N,tolerance)),-15)); colorbar('vert'); title(['KL log Lowerbound C_N N=' num2str(N)]); ylabel('q'); xlabel('delta');

plot(q_vec(1:length(conv_eps_vec)), conv_eps_vec, 'm');


figure; subplot(2,2,1); hold on; 
subplot(2,2,1); hold on; imagesc(delta_vec, q_vec, (C_N+cx_N)./2); colorbar('vert'); title(['KL Upperbound Plus Lowerbound Average C_N N=' num2str(N)  ' (p,\epsilon) = (' num2str(p) ' , ' num2str(eps), ')']); ylabel('q'); xlabel('delta');
subplot(2,2,2); hold on; imagesc(delta_vec, q_vec, C_N-cx_N); colorbar('vert'); title(['KL Upperbound Minus Lowerbound C_N N=' num2str(N)]); ylabel('q'); xlabel('delta');
subplot(2,2,3); hold on; imagesc(delta_vec, q_vec, 2.*(C_N-cx_N)./(C_N+cx_N)); colorbar('vert'); title('KL Relative Error Upperbound Minus Lowerbound / Average C_N' ); ylabel('q'); xlabel('delta');
subplot(2,2,4); hold on; imagesc(delta_vec, q_vec, min(max((C_N-cx_N)./(1 - (C_N+cx_N)./2), 0), 3)); colorbar('vert'); title('KL Relative Error Upperbound Minus Lowerbound / (1-Average) C_N'); ylabel('q'); xlabel('delta');



% Now try to plot how the upper and lower bounds approach the entropy 
figure; subplot(2,2,1); hold on;

subplot(2,2,1); hold on; imagesc(delta_vec, q_vec, (C_N-C_N_Plus1)); colorbar('vert'); title(['KL Upperbound N+1 Minus Upper C N=' num2str(N)  ' (p,\epsilon) = (' num2str(p) ' , ' num2str(eps), ')']); ylabel('q'); xlabel('delta');
subplot(2,2,2); hold on; imagesc(delta_vec, q_vec, (cx_N-cx_N_Plus1)); colorbar('vert'); title(['KL Lowerbound N+1 Minus Lower C N=' num2str(N)]); ylabel('q'); xlabel('delta');
subplot(2,2,3); hold on; imagesc(delta_vec, q_vec, 2.*(C_N-C_N_Plus1)./(C_N+C_N_Plus1)); colorbar('vert'); title(['KL Relative Error (Upperbound N+1 Minus Upper)/Average C N=' num2str(N)]); ylabel('q'); xlabel('delta');
subplot(2,2,4); hold on; imagesc(delta_vec, q_vec, 2.*(cx_N-cx_N_Plus1)./(cx_N+cx_N_Plus1)); colorbar('vert'); title(['KL Relative Error (Lowerbound N+1 Minus Lower)/Average C N=' num2str(N)]); ylabel('q'); xlabel('delta');



% subplot(2,2,1); hold on; imagesc((C_N_Plus1-cx_N_Plus1)); colorbar('vert'); title(['Upperbound Minus Lowerbound C N=' num2str(N+1)]); ylabel('p'); xlabel('eps');
% subplot(2,2,2); hold on; imagesc((C_N_Plus1-cx_N_Plus1)); colorbar('vert'); title(['Upperbound Minus Lowerbound C N=' num2str(N+1)]); ylabel('p'); xlabel('eps');
% 
% subplot(2,2,3); hold on; imagesc((C_N-C_N_Plus1)); colorbar('vert'); title(['Upperbound N+1 Minus Upper C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% subplot(2,2,4); hold on; imagesc((cx_N-cx_N_Plus1)); colorbar('vert'); title(['Lowerbound N+1 Minus Lower C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% 
% subplot(2,2,3); hold on; imagesc((Ido_upper-upper_ent)./2); colorbar('vert'); title('Upperbound Diff '); ylabel('p'); xlabel('eps');
% subplot(2,2,4); hold on; imagesc((Ido_lower-mul_ent)./2); colorbar('vert'); title('Lowerbound Diff'); ylabel('p'); xlabel('eps');



