% Remind what is every variable
genes_mu_vec = zeros(2,Ngenes);
all_pos_ind = find(R.Labels==1);
all_neg_ind = setdiff([1:Nsamples], all_pos_ind);

genes_mu_vec(1,:) = mean(R.dat(:,all_pos_ind), 2);
genes_mu_vec(2,:) = mean(R.dat(:,all_neg_ind), 2);
genes_sigma_vec(1,:) = std(R.dat(:,all_pos_ind), [], 2);
genes_sigma_vec(2,:) = std(R.dat(:,all_neg_ind), [], 2);


p = length(all_pos_ind) / ( length(all_pos_ind) + length(all_neg_ind) );

E_x_vec = p; E_y_vec = (1-p)*genes_mu_vec(2,:) + p*genes_mu_vec(1,:); E_xy_vec = p.*genes_mu_vec(1,:);
sigma_x_vec = sqrt(p*(1-p));

sigma_y_vec = sqrt(p .* (genes_sigma_vec(1,:).^2+genes_mu_vec(1,:).^2) + ...
               (1-p) .* (genes_sigma_vec(2,:).^2+genes_mu_vec(2,:).^2) - ...
               (1-p)^2 .*  genes_mu_vec(2,:).^2 - p^2 .* genes_mu_vec(1,:).^2 - ...
               2*p*(1-p) .* genes_mu_vec(2,:).*genes_mu_vec(1,:));


rho_vec = (E_xy_vec - E_x_vec .* E_y_vec) ./ (sigma_x_vec .*  sigma_y_vec);

Z_vec = atanh(rho_vec);
