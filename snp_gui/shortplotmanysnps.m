% Plot several SNPs to see how they behave ... 
I = floor(rand(1,7).* 57244);

color_vec = 'bgrkymc';

figure; hold on;
for i=1:length(I)
    plot(NormalizedSNPsCopyMatA(I(i),:), NormalizedSNPsCopyMatB(I(i),:), ['.' color_vec(i)]);
end
xlabel('A Intensity'); ylabel('B Intensity');






