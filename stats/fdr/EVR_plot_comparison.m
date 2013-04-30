% Script for plotting E[V/R] 
close all
clear all
clc

tic
rand_seed;
color='bgrmkcy';
m_vec=[500];
m0_m=[0.05:0.05:1];
q_vec=[0.05];
rho_vec=[0];
mu1_vec=[2.5];
for i_m=1:length(m_vec)
    m=m_vec(i_m);
    W=100*m;
    y=normrnd(0,1,W,m+1);
    x=zeros(W,m);
    for i_rho=1:length(rho_vec)
        rho=rho_vec(i_rho);
        for i_q=1:length(q_vec)
            q=q_vec(i_q);
            EVR_IBH=zeros(length(mu1_vec),length(m0_m));
            EVR_STH=zeros(length(mu1_vec),length(m0_m));
            EVR_BKY=zeros(length(mu1_vec),length(m0_m));
            EVR_orc=zeros(length(mu1_vec),length(m0_m));
            EVR_BH95=zeros(length(mu1_vec),length(m0_m));
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
%                     m0_hat=min(m*ones(W,1),2*sum(P,2)+sqrt(m/6/pi));
                    %                     F_BH95=fdr_proc_mat(Psort,q);
                    %                     F_BKY=fdr_BKY_mat(Psort,q);
                    F_IBH=fdr_IBH_mat_step_up(Psort,q); 
                    F_orc=fdr_proc_mat(Psort,q/m0_m(i_m0_m));
                    F_BKY=fdr_BKY_mat(Psort,q);
                    F_STH=fdr_STH_mat(Psort,q);
                    temp=zeros(W,1);
                    temp(find(F_orc>0))=Q(W*(F_orc(find(F_orc>0))-1)+find(F_orc>0));
                    EVR_orc(i_mu1,i_m0_m)=sum(temp)/W;
                    EVR_BH95(i_mu1,i_m0_m)=EVR_orc(i_mu1,i_m0_m)*m0_m(i_m0_m);
                    temp=zeros(W,1);
                    temp(find(F_IBH>0))=Q(W*(F_IBH(find(F_IBH>0))-1)+find(F_IBH>0));
                    EVR_IBH(i_mu1,i_m0_m)=sum(temp)/W;%*q/EVR_orc(i_mu1,i_m0_m);
                    temp=zeros(W,1);
                    temp(find(F_BKY>0))=Q(W*(F_BKY(find(F_BKY>0))-1)+find(F_BKY>0));
                    EVR_BKY(i_mu1,i_m0_m)=sum(temp)/W;%*q/EVR_orc(i_mu1,i_m0_m);
                    temp=zeros(W,1);
                    temp(find(F_STH>0))=Q(W*(F_STH(find(F_STH>0))-1)+find(F_STH>0));
                    EVR_STH(i_mu1,i_m0_m)=sum(temp)/W;%*q/EVR_orc(i_mu1,i_m0_m);
                end
            end
            figure(i_rho);
            set(gcf,'color','w');
%             subplot(2,2,i_q);
%             [C,h] = contour(m0_m,mu1_vec,EVR_IBH,'linewidth',2);
%             clabel(C,h);grid on;
            plot(m0_m,EVR_orc,'linewidth',2);grid on; hold on;
            plot(m0_m,EVR_BH95,'--r','linewidth',2);grid on; hold on;
            plot(m0_m,EVR_BKY,'g','linewidth',2);grid on; hold on;
            plot(m0_m,EVR_STH,'c','linewidth',2);grid on; hold on;
            plot(m0_m,EVR_IBH,'m','linewidth',2);grid on; hold on;
            set(gca,'fontname','times','fontweight','demi','fontsize',12');
            xlabel('m_0/m'); ylabel('E(V/R)');
            title(['m=',num2str(m_vec(i_m)),', q=',num2str(q),', \mu_1=',num2str(mu1_vec),', \rho=',num2str(rho)]);
            legend('ORC','BH95','BKY','STH','IBH');
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
