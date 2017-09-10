figs_dir = 'C:\Users\oorzu\Google Drive\HUJI\Teaching\StatisticalEvolutionaryGenomics_52887\Docs\figs';
kin = 0.2; dom = 0.3;  % set kinship and fraternity values;

p = 0.34; % allele frequency

p_fam = [1-4*kin+dom, 4*kin-2*dom,  dom];

y = [1 2 5.5]; y_AA = y(1); y_Aa = y(2); y_aa = y(3);
y_bar = p^2 * y_AA + 2*p*(1-p)* y_Aa + (1-p)^2 * y_aa;

sigma_sqr = p^2 * (y_AA-y_bar)^2 + 2*p*(1-p)* (y_Aa-y_bar)^2 + (1-p)^2 * (y_aa-y_bar)^2;
sigma_a_sqr = 2*p*(1-p)*(p*y_AA+(1-2*p)*y_Aa - (1-p)*y_aa)^2;
sigma_d_sqr = p^2*(1-p)^2 * (2*y_Aa - y_AA - y_aa)^2;

sigma_a_sqr + sigma_d_sqr - sigma_sqr


cov_fam = 2*kin*sigma_a_sqr + dom*sigma_d_sqr

cov_fam_long = p^2 * (p*p_fam(2)+p_fam(3))*y_AA^2 + 2*p*(1-p) * (0.5*p_fam(2)+p_fam(3)) * y_Aa^2 + (1-p)^2 * ((1-p)*p_fam(2)+p_fam(3)) * y_aa^2 + ...
    2*p^2*(1-p)*p_fam(2)*y_AA*y_Aa + 2*(1-p)^2*p*p_fam(2)*y_Aa*y_aa + ...
    (p_fam(1)-1) * y_bar^2

cov_fam - cov_fam_long


% Test HW:

init_p = {[0.1 0.2 0.5], [0.25 0.25 0.25]};

for similar_flag = 0:1
    if(similar_flag)
        similar_inbreeding=1; dissimilar_inbreeding=0; sim_str = 'similar';
    else
        similar_inbreeding=0; dissimilar_inbreeding=1; sim_str = 'dissimilar';
    end
    figure;
    for i_init = 1:2
        gen=20; p_AA = zeros(gen,1); p_Aa=p_AA; p_aa=p_AA; p=p_AA;
        p_AA(1) = init_p{i_init}(1);  p_Aa(1) = init_p{i_init}(2); p_aa(1) = init_p{i_init}(3);
        %    p_AA(1) = (sqrt(5)-1)/4+0.19; p_aa(1) = p_AA(1)+0.0000; p_Aa(1) = (1-p_AA(1)-p_aa(1))/2;
        for i=1:gen
            p(i) = p_AA(i) + p_Aa(i);
            p_AA(i+1) = p(i)^2;
            p_Aa(i+1) = p(i)*(1-p(i));
            p_aa(i+1) = (1-p(i))^2;
            
            if(similar_inbreeding)
                p_AA(i+1) = p_AA(i+1) / (1-2*p_AA(i)*p_aa(i));
                p_Aa(i+1) = (p_Aa(i+1)-p_AA(i)*p_aa(i)) / (1-2*p_AA(i)*p_aa(i));
                p_aa(i+1) = p_aa(i+1) / (1-2*p_AA(i)*p_aa(i));
            end
            
            if(dissimilar_inbreeding)
                p_AA(i+1) = (p_AA(i+1)-p_AA(i)^2) / (1-p_AA(i)^2-p_aa(i)^2);
                p_Aa(i+1) = p_Aa(i+1) / (1-p_AA(i)^2-p_aa(i)^2);
                p_aa(i+1) = (p_aa(i+1)-p_aa(i)^2) / (1-p_AA(i)^2-p_aa(i)^2);
            end
        end
        
        subplot(2,1,i_init); hold on;
        plot(p_AA, 'bo'); plot(p_Aa, 'r+'); plot(p_aa, 'g*'); 
        if(i_init==1)
            legend('AA', 'Aa', 'aa'); legend('boxoff');
        end
        title(['(' num2str(init_p{i_init}(1)) ', ' num2str(init_p{i_init}(2)) ', ' num2str(init_p{i_init}(3)) ')']); 
%        plot(p_AA + 2*p_Aa + p_aa, 'k--');
        
        
    end % loop on init conditions
    my_saveas(gcf, fullfile(figs_dir, [sim_str '_dynamics']), {'jpg', 'epsc'});
end % similar/dissimilar inbreeding

%x = 0:0.001:0.8;
%f_x = 1+2*x+2*x.^2 - (1-5*x+8*x.^2+8*x.^3);
%figure; plot(x, f_x);



syms pAA pAa paa
eqn = [pAA*(1-pAA^2-paa^2) == pAa*(2*pAA+pAa), pAa*(1-pAA^2-paa^2) == (pAA+pAa)*(paa+pAa), paa*(1-pAA^2-paa^2) == pAa*(2*paa+pAa), pAA+2*pAa+paa==1];
S = solve(eqn, [pAA, pAa, paa]);

%fsolve(z^3 - z^2/2 - z/2 + 1/8, z)
RR = roots([1, -1/2, -1/2, 1/8])


