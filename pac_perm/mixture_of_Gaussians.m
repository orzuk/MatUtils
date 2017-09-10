% Temp MoG function (not needed anymore since we have a package) 
function [S,M,P, LogLike]=mixture_of_Gaussians(x,num_of_Gaussians,num_of_iterations,N_samples)

min_x=min(x);
max_x=max(x);
gap_x=max_x-min_x;
step_x=gap_x/(num_of_Gaussians-1);
miu=min_x:step_x:max_x;
% sigma=ones(1,num_of_Gaussians)/sqrt(N_samples-3);
prior=ones(1,num_of_Gaussians)/num_of_Gaussians;
N=length(x);
p=ones(N,num_of_Gaussians)/num_of_Gaussians;

maxKL=-1000000000;
for replica=1:20
    sigma=ones(1,num_of_Gaussians).*rand(1,num_of_Gaussians)*step_x;
    for itt=1:num_of_iterations
        % % % % %         rep_sigma=repmat(sigma,N,1);
        % % % % %         p=1./(rep_sigma*sqrt(2*pi)).*...
        % % % % %             exp(repmat(x,1,num_of_Gaussians)-repmat(miu,N,1)).^2./(2*rep_sigma.^2);
        % % % % %         z=p.*repmat(prior,N,1)./repmat(p*prior',1,num_of_Gaussians);


        % % % %                 for i=1:N
        % % % %                     for m=1:num_of_Gaussians
        % % % %                         if(sigma(m))
        % % % %                             p(i,m)=1/(sqrt(2*pi)*sigma(m))*exp(-(x(i)-miu(m))^2/(2*sigma(m)^2));
        % % % %                         else
        % % % %                             p(i,m)=1;
        % % % %                         end
        % % % %                         z(i,m)=p(i,m)*prior(m)/(p(i,:)*prior');
        % % % %                     end
        % % % %                 end

        for m=1:num_of_Gaussians
            if(sigma(m))
                p(:,m)=1./(sqrt(2.*pi).*sigma(m)).*exp(-(x-miu(m)).^2./(2.*sigma(m).^2));
            else
                p(:,m)=1;
            end
            z(:,m)=(p(:,m).*prior(m))./((p*prior'));
        end



        %     sum_z=sum(z);
        %     for m=1:num_of_Gaussians
        %         sigma(m)=sqrt(z(:,m)'*(x-miu(m)).^2)/sum_z(m);           %x is a column vector
        %         miu(m)=z(:,m)'*x/sum_z(m);
        %         prior(m)=sum_z(m)/N;
        %     end
        sum_z=sum(z);
        sigma_old=sigma;
        miu_old=miu;
        prior_old=prior;
        sum_z(sum_z==0)=1;
        sigma=sqrt(sum(z.*(repmat(x,1,num_of_Gaussians)-repmat(miu,N,1)).^2)./sum_z);           %x is a column vector
        miu=(z'*x)'./sum_z;
        prior=sum_z/N;
        %         norm(sigma_old-sigma)+norm(miu_old-miu)+norm(prior_old-prior)

%         clear y
%         for m=1:num_of_Gaussians
%             y(m,:)=prior(m)*1/(sqrt(2*pi)*sigma(m))*exp(-(x-miu(m)).^2/(2*sigma(m)^2));
%         end
%         KL=sum(log(sum(y,1)))


    end

    clear y
    for m=1:num_of_Gaussians
        y(m,:)=prior(m)*1/(sqrt(2*pi)*sigma(m))*exp(-(x-miu(m)).^2/(2*sigma(m)^2));
    end
    KL=sum(log(sum(y,1)))
    if(KL>maxKL)
        maxKL=KL;
        P=prior;
        S=sigma;
        M=miu;
    end
end

% return also the best log-likelihood
LogLike = maxKL;
