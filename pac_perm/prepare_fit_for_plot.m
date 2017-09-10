% A temp script fitting correlations to plot 
[hieght,bin_loc]=hist(Fisher_Zs,250);
figure,hist(Fisher_Zs,250)
clear y

for i=1:num_of_Gaussians
    y(i,:)=prior(i)*1/(sqrt(2*pi)*Fisher_Zs_std(i))*exp(-(bin_loc-miu(i)).^2/(2*Fisher_Zs_std(i)^2));
end
hold on,
plot(bin_loc,sum(y,1)*(bin_loc(2)-bin_loc(1))*5000,'r-')
