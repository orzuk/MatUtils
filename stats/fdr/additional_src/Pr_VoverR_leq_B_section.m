% Compute and plot cross-sections for Prob. of V/R >= B 
clear all
tic
rand_seed;
color='bgrmkcy';
m_vec=[1000];  % m
res = 0.001; m0_m=[res:res:1]; % m0/m vector
q_vec=[0.05];  % q (why needed?)
rho_vec=[0,0.5];  % rho (correlations)
mu1_vec=[2.5];   % how far apart are the null and non-null hypothesis
figure;
for i_m=1:length(m_vec)
    m=m_vec(i_m);
    W=50*m; % number of simulations
    y=single(randn(W,m+1)); % draw Gaussians
%     x=zeros(W,m);
    for i_rho=1:length(rho_vec)  % loop on rho
        rho=rho_vec(i_rho);
        % prepare enough p-values for both true and false hypothesis
        x_true =sqrt(rho)*repmat(y(:,m+1),1,m)+sqrt(1-rho)*y(:,1:m);
        P_true = single(2*normcdf(-abs(x_true)));


        for i_q=1:length(q_vec)  % loop on q
            q=q_vec(i_q);
            Pr_VR_leq_m0q_m=zeros(length(mu1_vec),length(m0_m));
            Pr_VR_leq_m0hatq_m=zeros(length(mu1_vec),length(m0_m));
            V = {}; V_over_R = {}; m0_hat = {}; 
            for i_mu1=1:length(mu1_vec)    % loop on mu1
                mu1=mu1_vec(i_mu1);
                P_false = single(2*normcdf(-abs(x_true+mu1)));
                for i_m0_m=1:length(m0_m)   % loop on m0/m
                    num_tru=floor(m*m0_m(i_m0_m));
                    num_false=m-num_tru;
%                      x(:,1:num_tru)=sqrt(rho)*repmat(y(:,m+1),1,num_tru)+sqrt(1-rho)*y(:,1:num_tru);
%                      x(:,num_tru+1:end)=sqrt(rho)*repmat(y(:,m+1),1,num_false)+sqrt(1-rho)*y(:,num_tru+1:end-1)+mu1;
%                    x = x_true; x(:,num_tru+1:end) = x(:,num_tru+1:end)+mu1;                    
                   %   P = 2*normcdf(-abs(x));
                    [Psort,XI]=sort( [P_true(:,1:num_tru) P_false(:,num_tru+1:end)], 2 ); % Generated p-values 
%                     out=zeros(W,m);
%                     out(find(XI<=num_tru))=1;
%                     out = (XI
                   
                    m0_hat{i_m0_m}=2*sum(Psort,2);
                    %                     F_BH95=fdr_proc_mat(Psort,q);
                    %                     F_BKY=fdr_BKY_mat(Psort,q);
                    %                     F_IBH=fdr_IBH_mat_step_up(Psort,q);
                    F_BH95=fdr_proc_mat(Psort,q);
                    V{i_m0_m} = zeros(W,1);
%                     Q=cumsum(XI<=num_tru,2)./repmat([1:m],W,1);
%                     temp(find(F_BH95>0))=Q(W*(F_BH95(find(F_BH95>0))-1)+find(F_BH95>0));
                    for j=1:W                     
                        V{i_m0_m}(j) = sum(XI(j,1:F_BH95(j)) <= num_tru);
                    end
                    V_over_R{i_m0_m} = V{i_m0_m} ./ max(F_BH95, 1);
                    Pr_VR_leq_m0q_m(i_mu1,i_m0_m)=sum((V_over_R{i_m0_m}  <=m0_m(i_m0_m)*q)); % how many times V/R <= q
                    Pr_VR_leq_m0hatq_m(i_mu1,i_m0_m)=sum(((V_over_R{i_m0_m} -(m0_hat{i_m0_m}/m)*q)<=0));
                    %                     EVR_orc(i_mu1,i_m0_m)=sum(temp)/W;
                    %                     temp=zeros(W,1);
                    %                     temp(find(F_IBH>0))=Q(W*(F_IBH(find(F_IBH>0))-1)+find(F_IBH>0));
                    %                     EVR_IBH(i_mu1,i_m0_m)=sum(temp)/W*q/EVR_orc(i_mu1,i_m0_m);
                end
            end
            Pr_VR_leq_m0q_m = Pr_VR_leq_m0q_m ./ W;
            Pr_VR_leq_m0hatq_m = Pr_VR_leq_m0hatq_m ./ W;
            
            if i_rho==1
                set(gcf,'color','w');
                subplot(1,2,i_q+2*(i_m-1));
                plot(m0_m,Pr_VR_leq_m0q_m,'linewidth',2);grid on;hold on;
                plot(m0_m,Pr_VR_leq_m0hatq_m,'m--','linewidth',2); grid on;
                set(gca,'fontname','times','fontweight','demi','fontsize',12');
                xlabel('m_0/m'); ylabel('Pr(V/R<=B)');
                title(['m=',num2str(m_vec(i_m)),', \mu_1=',num2str(mu1),', q=',num2str(q),', \rho=',num2str(rho),', BH95']);
                legend('B=q','B=m_{0,hat}q/m');
            else
                subplot(1,2,i_q+1+2*(i_m-1));
                plot(m0_m,Pr_VR_leq_m0q_m,'linewidth',2);grid on;hold on;
                plot(m0_m,Pr_VR_leq_m0hatq_m,'m--','linewidth',2); grid on;
               %% plot(m0_m,m0_m *q ./ m ,'r:','linewidth',2); grid on;
                
                set(gca,'fontname','times','fontweight','demi','fontsize',12');
                xlabel('m_0/m'); ylabel('Pr(V/R<=B)');
                title(['m=',num2str(m_vec(i_m)),', \mu_1=',num2str(mu1),', q=',num2str(q),', \rho=',num2str(rho),', BH95']);
            end
        end
    end
end

saveas(gcf, ['Pr_V_R_leq_B_m_' num2str(m) '_res_' num2str(res) '_iters_' num2str(W) '.fig']);   % save figure as .fig
saveas(gcf, ['Pr_V_R_leq_B_m_' num2str(m) '_res_' num2str(res) '_iters_' num2str(W) '.jpg']);   % save figure as .jpg
save(['V_R_m_' num2str(m) '_res_' num2str(res) '_iters_' num2str(W) '.mat'], 'V', 'V_over_R', 'm0_hat');

toc
