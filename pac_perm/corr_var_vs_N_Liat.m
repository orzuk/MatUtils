% A short script testing variance vs. correlation 
%load Rosetta_data.mat
data=R.dat;
N_TOP=267;
N_genes=size(data,1);
N_samples=size(data,2);
idx_poor=find(samples==1);
idx_good=find(samples==0);
N_poor=length(idx_poor);
N_good=N_samples-N_poor;
ratio=N_poor/N_samples;
N_iterations=100;
min_i=4;
max_i=floor(N_samples/2);
group_sizes=[min_i:5:max_i];
num_of_sizes=length(group_sizes);
mean_corr=zeros(N_genes,num_of_values);
var_corr=zeros(N_genes,num_of_values);
f=zeros(1,num_of_values);
for i=1:num_of_sizes
    group_size=group_size_values(i);
    clear corr curr_samples;
    N_groups=floor(N_samples/group_size);%number of groups
    N_couples=N_groups*(N_groups-1)/2;
    curr_N_poor(1:N_groups)=floor(N_poor/N_groups);%number of poor samples in each group
    res_poor=mod(N_poor,N_groups);
    res_good=mod(N_good,N_groups);
    if(res_poor+res_good>=N_groups)
        curr_N_poor(1:res_poor)=curr_N_poor(1:res_poor)+1;%the last group 
        curr_N_good=group_size-curr_N_poor;
    else
       curr_N_good(1:N_groups)=floor(N_good/N_groups); 
    end
    clear S
    for j=1:N_iterations
        j
        x=randperm(N_poor);
        y=randperm(N_good);
        start_poor=1;
        start_good=1;
        S(1).idx=[idx_poor(x(start_poor:curr_N_poor(1)));...
                idx_good(y(start_good:curr_N_good(1)))];
        for k=2:N_groups
            start_poor=start_poor+curr_N_poor(k-1);
            start_good=start_good+curr_N_good(k-1);
            S(k).idx=[idx_poor(x(start_poor:start_poor+curr_N_poor(k)-1));...
                idx_good(y(start_good:start_good+curr_N_good(k)-1))];
        end
        ff=0;
        for k=1:N_groups
            [norm_data,mean_genes,norm_genes]=dna_N_new(data(:,S(k).idx));
            curr_samples=samples(S(k).idx);
            norm_samples=curr_samples-mean(curr_samples);
            norm_samples=norm_samples/sqrt(sum(norm_samples.^2));
            corr(1:N_genes,k)=norm_data*norm_samples;
            [val,ind]=sort(-abs(corr(:,k)));
            list(1:N_TOP,k)=ind(1:N_TOP);
            for kk=1:k-1
                ff=ff+length(intersect(list(:,k),list(:,kk)));
            end
        end
        f(i)=f(i)+ff/(N_couples);
        mean_corr(:,i)=mean_corr(:,i)+mean(corr,2);
        var_corr(:,i)=var_corr(:,i)+mean(corr.^2,2)-mean(corr,2).^2;
    end
end
mean_corr=mean_corr/N_iterations;
var_corr=var_corr/N_iterations;
f=f/(N_iterations*N_TOP);    
