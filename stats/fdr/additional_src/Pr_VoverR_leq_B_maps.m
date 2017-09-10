% Compute and plot 2d maps for Prob. of V/R >= B 
close all
clear all
clc

tic
rand_seed;
color='bgrmkcy';
m_vec=[100,500];
m0_m=[0.05:0.05:1];
q_vec=[0.05];
rho_vec=[0];
mu1_vec=[1:0.2:4];
for i_m=1:length(m_vec)
    m=m_vec(i_m);
    W=100*m;
    y=normrnd(0,1,W,m+1);
    x=zeros(W,m);
    for i_rho=1:length(rho_vec)
        rho=rho_vec(i_rho);
        for i_q=1:length(q_vec)
            q=q_vec(i_q);
            Pr_VR_leq_m0q_m=zeros(length(mu1_vec),length(m0_m));
            Pr_VR_leq_m0hatq_m=zeros(length(mu1_vec),length(m0_m));
            for i_mu1=1:length(mu1_vec)
                mu1=mu1_vec(i_mu1);
                for i_m0_m=1:length(m0_m)
                    num_tru=floor(m*m0_m(i_m0_m));
                    num_false=m-num_tru;
                    x(:,1:num_tru)=sqrt(rho)*repmat(y(:,m+1),1,num_tru)+sqrt(1-rho)*y(:,1:num_tru);
                    x(:,num_tru+1:end)=sqrt(rho)*repmat(y(:,m+1),1,num_false)+sqrt(1-rho)*y(:,num_tru+1:end-1)+mu1;
                    P=2*normcdf(-abs(x),0,1);
                    [Psort,XI]=sort(P,2);
                    out=zeros(W,m);
                    out(find(XI<=num_tru))=1;
                    Q=cumsum(out,2)./repmat([1:m],W,1);
                    m0_hat=2*sum(P,2);
                    %                     F_BH95=fdr_proc_mat(Psort,q);
                    %                     F_BKY=fdr_BKY_mat(Psort,q);
%                     F_IBH=fdr_IBH_mat_step_up(Psort,q); 
                    F_BH95=fdr_proc_mat(Psort,q);
                    temp=zeros(W,1);
                    temp(find(F_BH95>0))=Q(W*(F_BH95(find(F_BH95>0))-1)+find(F_BH95>0));
                    Pr_VR_leq_m0q_m(i_mu1,i_m0_m)=length(find(temp<=m0_m(i_m0_m)*q))/W;
                    Pr_VR_leq_m0hatq_m(i_mu1,i_m0_m)=length(find((temp-m0_hat/m*q)<=0))/W;
%                     EVR_orc(i_mu1,i_m0_m)=sum(temp)/W;                                       
%                     temp=zeros(W,1);
%                     temp(find(F_IBH>0))=Q(W*(F_IBH(find(F_IBH>0))-1)+find(F_IBH>0));
%                     EVR_IBH(i_mu1,i_m0_m)=sum(temp)/W*q/EVR_orc(i_mu1,i_m0_m);
                end
            end
            figure(i_rho);
            set(gcf,'color','w');
            subplot(2,2,i_q+2*(i_m-1));
            [C,h] = contour(m0_m,mu1_vec,Pr_VR_leq_m0q_m,[0.1:0.1:1],'linewidth',2);
            clabel(C,h);grid on;
            set(gca,'fontname','times','fontweight','demi','fontsize',12');
            xlabel('m_0/m'); ylabel('\mu_1');
            title(['Pr(V/R<=m_0q/m), m=',num2str(m_vec(i_m)),', q=',num2str(q),', \rho=',num2str(rho),' BH95']);
            subplot(2,2,i_q+1+2*(i_m-1));
            [C,h] = contour(m0_m,mu1_vec,Pr_VR_leq_m0hatq_m,[0.1:0.1:1],'linewidth',2);
            clabel(C,h);grid on;
            set(gca,'fontname','times','fontweight','demi','fontsize',12');
            xlabel('m_0/m'); ylabel('\mu_1');
            title(['Pr(V/R<=m_0^{hat}q/m), m=',num2str(m_vec(i_m)),', q=',num2str(q),', \rho=',num2str(rho),' step up']);
%             figure(i_rho+1);
%             set(gcf,'color','w');
%             subplot(2,2,i_q);
%             [C,h] = contour(m0_m,mu1_vec,EVR_orc,'linewidth',2);
%             clabel(C,h);grid on;
%             set(gca,'fontname','times','fontweight','demi','fontsize',12');
%             xlabel('m_0/m'); ylabel('\mu_1');
%             title(['E(V/R)_{oracle}, m=',num2str(m_vec(i_m)),', q=',num2str(q),', \rho=',num2str(rho)]);
%             figure(i_rho+2);
%             set(gcf,'color','w');
%             subplot(2,2,i_q);
%             [C,h] = contour(m0_m,mu1_vec,EVR_orc-EVR_IBH,'linewidth',2);
%             clabel(C,h);grid on;
%             set(gca,'fontname','times','fontweight','demi','fontsize',12');
%             xlabel('m_0/m'); ylabel('\mu_1');
%             title(['E(V/R)_{oracle-IBH}, m=',num2str(m_vec(i_m)),', q=',num2str(q),', \rho=',num2str(rho)]);
        end
    end
end

toc
