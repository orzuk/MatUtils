% Test when two inequalities for FDR hold 
close all
clear all
clc

tic
rand_seed;
color='bgrmkcy';
m_vec=[500];
m0_m=[0.05:0.05:1];
q_vec=[0.02,0.05,0.1,0.2];
rho_vec=[0];
mu1_vec=[1:0.2:4];
for i_m=1:length(m_vec)
    m=m_vec(i_m);
    c_vec=[1.099757389219602
        1.080178948712779
        1.068806352812296
        1.061166943357300
        1.055589032386557
        1.051288991415590
        1.047844786256652
        1.045006533413207
        1.042615688768783
        1.040566207513654
        1.038784140148840
        1.037216146320969
        1.035822660423759
        1.034573635848637
        1.033445792217690
        1.032420776342130
        1.031483898991471
        1.030623245956020
        1.029829039083024];
    c_corr=interp1([100:50:1000],c_vec,m,'spline');
    W=10*m;
    y=normrnd(0,1,W,m+1);
    x=zeros(W,m);
    for i_rho=1:length(rho_vec)
        rho=rho_vec(i_rho);
        for i_q=1:length(q_vec)
            q=q_vec(i_q);
            %             EVR_IBH=zeros(length(mu1_vec),length(m0_m));
            Pr1and2=zeros(length(mu1_vec),length(m0_m));
            Pr1and2no=zeros(length(mu1_vec),length(m0_m));
            Pr1noand2=zeros(length(mu1_vec),length(m0_m));
            Pr1noand2no=zeros(length(mu1_vec),length(m0_m));
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
                    %                     F_IBH=fdr_IBH_mat_step_up(Psort,q);
                    m0_hat=c_corr*min(m*ones(W,1),2*sum(Psort,2));
                    F_BH95=fdr_proc_mat(Psort,q);
                    temp=zeros(W,1);
                    temp(find(F_BH95>0))=Q(W*(F_BH95(find(F_BH95>0))-1)+find(F_BH95>0));
                    Pr1and2(i_mu1,i_m0_m)=length(find(temp<=m0_m(i_m0_m)*q & m0_hat>=m*m0_m(i_m0_m)))/W;
                    Pr1and2no(i_mu1,i_m0_m)=length(find(temp<=m0_m(i_m0_m)*q & m0_hat<m*m0_m(i_m0_m)))/W;
                    Pr1noand2(i_mu1,i_m0_m)=length(find(temp>m0_m(i_m0_m)*q & m0_hat>=m*m0_m(i_m0_m)))/W;
                    Pr1noand2no(i_mu1,i_m0_m)=length(find(temp>m0_m(i_m0_m)*q & m0_hat<m*m0_m(i_m0_m)))/W;
                    %                     temp=zeros(W,1);
                    %                     temp(find(F_IBH>0))=Q(W*(F_IBH(find(F_IBH>0))-1)+find(F_IBH>0));
                    %                     EVR_IBH(i_mu1,i_m0_m)=sum(temp)/W*q/EVR_orc(i_mu1,i_m0_m);
                end
            end
            figure(i_rho);
            set(gcf,'color','w');
            subplot(2,2,i_q);
            [C,h] = contour(m0_m,mu1_vec,Pr1and2./(Pr1and2+Pr1and2no),'linewidth',2);
            clabel(C,h);grid on;
            set(gca,'fontname','times','fontweight','demi','fontsize',12');
            xlabel('m_0/m'); ylabel('\mu_1');
            title(['Pr(m_0<= m_{0,hat} | V/R <= qm_0/m), m='...
                ,num2str(m_vec(i_m)),', q=',num2str(q),', \rho=',num2str(rho),' step up']);
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
