% Compare lower and upper bounds on HMM entropy 

p = 0.3;
eps = 0.1;

segs = 1000;

tolerance = 0.000000000001;
mul_ent = zeros(segs-1,segs-1);
upper_ent = zeros(segs-1,segs-1);



% New ! We here do not let p and epsilon be 0, 0.5
p_vec = [0:0.5/(segs):0.5]';
eps_vec = [0:0.5/(segs):0.5]';
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

for i=1:segs-1  % loop over p
    i
    mul_ent(i,:) = entropy( [p_vec(i) + eps_vec - 2*p_vec(i)*eps_vec, 1 - (p_vec(i) + eps_vec - 2*p_vec(i)*eps_vec)]' );
    temp = p_vec(i) + 2*eps_vec - 4*p_vec(i)*eps_vec + 4*p_vec(i)*(eps_vec .^ 2) - 2* (eps_vec .^ 2);
    if(~isempty(find(temp > 1)))
        eps_vec
        ppp = p_vec(i)
        temp
    end
    
    upper_ent(i,:) = entropy( [temp, 1-temp]' );
    
    p = min(max(p_vec(i), tolerance), 1-tolerance);
    
   % Do the upper bound H(Y_n | Y_{n-1},Y_{n-2})
   IDO_PA = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*(eps_vec(i)^3 + (1-eps_vec(i))^3 + eps_vec(i)*(1-eps_vec(i))) + p_vec.^2*eps_vec(i)*(1-eps_vec(i));
   IDO_PB = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + 2*p_vec.*(1-p_vec)*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*(eps_vec(i)^3 + (1-eps_vec(i))^3);
   IDO_PC = (1-p_vec).^2*(eps_vec(i)^3 + (1-eps_vec(i))^3) + 2*p_vec.*(1-p_vec)*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*eps_vec(i)*(1-eps_vec(i));
   Ido_upper(:,i) = -( IDO_PA.*log2(IDO_PA) + 0.5*IDO_PB.*log2(IDO_PB) + 0.5*IDO_PC.*log2(IDO_PC));
   
   
   % Now do the lower bound H(Y_n | Y_{n-1},Y_{n-2}, X_{n-1})
   IDO_PA = (1-p_vec).^2*(1-eps_vec(i))^2 + p_vec.^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*eps_vec(i);
   IDO_PB = (1-p_vec).^2*eps_vec(i)^2 + p_vec.^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*(1-eps_vec(i));
   IDO_PC = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*(1-eps_vec(i))^2 + p_vec.*(1-p_vec)*eps_vec(i);
   IDO_PD = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*eps_vec(i)^2 + p_vec.*(1-p_vec)*(1-eps_vec(i));
   Ido_lower(:,i) = -0.5*( IDO_PA.*log2(IDO_PA) + IDO_PB.*log2(IDO_PB) + IDO_PC.*log2(IDO_PC) + IDO_PD.*log2(IDO_PD));
   
end


N=4; % Set the number of bits to use in the upper/lower bounds
% Now loop over eps


for i=1:segs-1  % loop over p
    if(mod(i, 100) == 0)
        i
    end
            
    T_N = HMP_entropy_finite(eps_vec(i), p_vec', N);
    tx_N = HMP_entropy_finite_X1(eps_vec(i), p_vec', N);
    C_N(:,i) = T_N-HMP_entropy_finite(eps_vec(i), p_vec', N-1); % The upper bound
    cx_N(:,i) = tx_N-HMP_entropy_finite_X1(eps_vec(i), p_vec', N-1); % The lower bound
    C_N_Plus1(:,i) = HMP_entropy_finite(eps_vec(i), p_vec', N+1)-T_N; % The upper bound N+1
    cx_N_Plus1(:,i) = HMP_entropy_finite_X1(eps_vec(i), p_vec', N+1)-tx_N; % The lower bound N+1
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


% New figures : Only upper and lower bounds 
figure; subplot(2,2,1); hold on; 
subplot(2,2,1); hold on; imagesc( eps_vec, p_vec, (upper_ent+mul_ent)./2 ); colorbar('vert'); title('HMM Lower+Upper Bound Average [H(Y_n/Y_{n-1})+H(Y_n/X_{n-1})]/2'); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(eps_vec, p_vec,(Ido_upper-Ido_lower)); colorbar('vert'); title('Upperbound Minus Lowerbound'); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(eps_vec, p_vec,(C_N-cx_N)); colorbar('vert'); title(['Upperbound Minus Lowerbound C N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(eps_vec, p_vec,(C_N-Ido_upper)); colorbar('vert'); title('Upperbound C Minus Ido Upperbound '); ylabel('p'); xlabel('eps');



figure; subplot(2,2,1); hold on; 
subplot(2,2,1); hold on; imagesc(eps_vec, p_vec, (C_N+cx_N)./2); colorbar('vert'); title(['Upperbound Plus Lowerbound Average C_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(eps_vec, p_vec, max(C_N-cx_N, 0)); colorbar('vert'); title(['Upperbound Minus Lowerbound C_N N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(eps_vec, p_vec, max(2.*(C_N-cx_N)./(C_N+cx_N), 0)); colorbar('vert'); title('Relative Error Upperbound Minus Lowerbound / Average C_N' ); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(eps_vec, p_vec, min(max((C_N-cx_N)./(1 - (C_N+cx_N)./2), 0), 3)); colorbar('vert'); title('Relative Error Upperbound Minus Lowerbound / (1-Average) C_N'); ylabel('p'); xlabel('eps');



% Now try to plot how the upper and lower bounds approach the entropy 
figure; subplot(2,2,1); hold on;

subplot(2,2,1); hold on; imagesc(eps_vec, p_vec, (C_N-C_N_Plus1)); colorbar('vert'); title(['Upperbound N+1 Minus Upper C N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(eps_vec, p_vec, (cx_N-cx_N_Plus1)); colorbar('vert'); title(['Lowerbound N+1 Minus Lower C N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(eps_vec, p_vec, 2.*(C_N-C_N_Plus1)./(C_N+C_N_Plus1)); colorbar('vert'); title(['Relative Error (Upperbound N+1 Minus Upper)/Average C N=' num2str(N)]); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(eps_vec, p_vec, 2.*(cx_N-cx_N_Plus1)./(cx_N+cx_N_Plus1)); colorbar('vert'); title(['Relative Error (Lowerbound N+1 Minus Lower)/Average C N=' num2str(N)]); ylabel('p'); xlabel('eps');



% subplot(2,2,1); hold on; imagesc((C_N_Plus1-cx_N_Plus1)); colorbar('vert'); title(['Upperbound Minus Lowerbound C N=' num2str(N+1)]); ylabel('p'); xlabel('eps');
% subplot(2,2,2); hold on; imagesc((C_N_Plus1-cx_N_Plus1)); colorbar('vert'); title(['Upperbound Minus Lowerbound C N=' num2str(N+1)]); ylabel('p'); xlabel('eps');
% 
% subplot(2,2,3); hold on; imagesc((C_N-C_N_Plus1)); colorbar('vert'); title(['Upperbound N+1 Minus Upper C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% subplot(2,2,4); hold on; imagesc((cx_N-cx_N_Plus1)); colorbar('vert'); title(['Lowerbound N+1 Minus Lower C N=' num2str(N)]); ylabel('p'); xlabel('eps');
% 
% subplot(2,2,3); hold on; imagesc((Ido_upper-upper_ent)./2); colorbar('vert'); title('Upperbound Diff '); ylabel('p'); xlabel('eps');
% subplot(2,2,4); hold on; imagesc((Ido_lower-mul_ent)./2); colorbar('vert'); title('Lowerbound Diff'); ylabel('p'); xlabel('eps');



