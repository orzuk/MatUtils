% Compare RLMM intensities for many populations
function Dummy = HapmapComparePopulationsRLMMParams(hapmap_population1, hapmap_population2, chip_type)

AssignAllGlobalConstants;

POP1 = load(['..\database\RLMM_' pop_str_vec{hapmap_population1} '_' chip_type '.mat']);
POP2 = load(['..\database\RLMM_' pop_str_vec{hapmap_population2} '_' chip_type '.mat']);
%YRI = load('..\database\RLMM_YRI_xba.mat');
%JPT_CHB = load('..\database\RLMM_JPT_CHB_xba.mat');

figure; title_str = [pop_str_vec{hapmap_population1} ' vs. ' pop_str_vec{hapmap_population2} ' on ']; 
subplot(2,3,1);
plot(POP1.RLMM.MuMats.AA(:,1), POP2.RLMM.MuMats.AA(:,1), '.'); title([title_str 'AA - A Intensity ' chip_type]); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});
subplot(2,3,4);
plot(POP1.RLMM.MuMats.AA(:,2), POP2.RLMM.MuMats.AA(:,2), '.'); title([title_str 'AA - B Intensity']); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});
subplot(2,3,2);
plot(POP1.RLMM.MuMats.AB(:,1), POP2.RLMM.MuMats.AB(:,1), '.'); title([title_str 'AB - A Intensity']); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});
subplot(2,3,5);
plot(POP1.RLMM.MuMats.AB(:,2), POP2.RLMM.MuMats.AB(:,2), '.'); title([title_str 'AB - B Intensity']); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});
subplot(2,3,3);
plot(POP1.RLMM.MuMats.BB(:,1), POP2.RLMM.MuMats.BB(:,1), '.'); title([title_str 'BB - A Intensity']); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});
subplot(2,3,6);
plot(POP1.RLMM.MuMats.BB(:,2), POP2.RLMM.MuMats.BB(:,2), '.'); title([title_str 'BB - B Intensity']); 
xlabel(pop_str_vec{hapmap_population1}); ylabel(pop_str_vec{hapmap_population2});


TTT = 98;

Dummy = 0;


