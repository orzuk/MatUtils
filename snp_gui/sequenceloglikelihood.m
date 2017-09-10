% The function computes the log-likelihood of a given (X,Y) sequence, and a
% given HMP model
function LL = SequenceLogLikelihood(X, Y, Y2, HMP)

DISCRETE = 0; GAUSSIAN = 1;

n = length(Y); % the vector length
HMP.Y_TYPE
% This is the loglikelihood of the first symbol
if(HMP.SPECIAL_MODELS_FLAG == 1)  % Here deal with the SNPs special model ...
    % Transfer X to its 'seperate' values
    X_alpha_genotype = bitget(X,1);
    X_beta_genotype = bitget(X,5);
    X_alpha_copynumber = bitget(X,2) + 2*bitget(X,3) + 4*bitget(X,4);
    X_beta_copynumber = bitget(X,6) + 2*bitget(X,7) + 4*bitget(X,8);
    X_copy = X_alpha_copynumber+X_beta_copynumber;
    X_A = X_alpha_genotype.*X_alpha_copynumber + ...
        X_beta_genotype.*X_beta_copynumber;
    X_B = (1-X_alpha_genotype).*X_alpha_copynumber + ...
        (1-X_beta_genotype).*X_beta_copynumber;

%    figure; hold on; plot(X_A-0.05, '.'); plot(X_B+0.05, 'r.'); plot(X_copy, 'g.');
        
    % The copy number contribution
    LL = log(HMP.PI(X_alpha_copynumber(1)+1)) + log(HMP.PI(X_beta_copynumber(1)+1));
    for i=2:n
        LL = LL + log(HMP.M(X_alpha_copynumber(i-1)+1, X_alpha_copynumber(i)+1)) + ...
                  log(HMP.M(X_beta_copynumber(i-1)+1, X_beta_copynumber(i)+1));
    end
    LL_COPY = LL;
    % The SNP genotype contribution (First we have SNP contribution, negligible ...)
    LL_SNP = 0;
    for i=2:n
        LL_SNP = LL_SNP + log(HMP.PLACE_M(X_alpha_genotype(i-1)+2*X_alpha_genotype(i)+1, i-1)) + ...
            log(HMP.PLACE_M(X_beta_genotype(i-1)+2*X_beta_genotype(i)+1, i-1));
    end
    LL_SNP;
    LL = LL + LL_SNP;
else
    LL = log(HMP.PI(X(1)));
    for i=2:n
        LL = LL + log(HMP.M(X(i-1),X(i)));
    end
end


% Separate Discrete and Gaussian cases
if(HMP.Y_TYPE == DISCRETE)
    for i=1:n
        LL = LL + log(HMP.N(X(i),Y(i)));
    end
else % Gaussian
    if(HMP.SPECIAL_MODELS_FLAG == 1)  % Here deal with the SNPs special model ...
        LL_Y = 0; LL_Y_end = 0;
        for i=1:n
            local_prob = 0; local_prob2 = 0;
            for j=1:HMP.y_dim
                local_prob = local_prob + HMP.N(1,j) * exp( -(Y2(i) - HMP.MU(X_A(i)+1,j))^2 / (2*HMP.SIGMA(j))^2 )  / sqrt(2*pi*HMP.SIGMA(j));
                local_prob2 = local_prob2 + HMP.N(1,j) * exp( -(Y(i) - HMP.MU(X_B(i)+1,j))^2 / (2*HMP.SIGMA(j))^2 )  / sqrt(2*pi*HMP.SIGMA(j));
            end
            if(  (i < 2500) && (i > 2000) ) 
                LL_Y = LL_Y + log(local_prob) + log(local_prob2);
            else
                LL_Y_end = LL_Y_end + log(local_prob) + log(local_prob2);
            end
        end
%         figure; subplot(2,1,1); hold on; plot(Y(2000:2500), '.'); plot(X_B(2000:2500), 'r.'); title('B');
%         subplot(2,1,2); hold on; plot(Y2(2000:2500), '.'); plot(X_A(2000:2500), 'r.'); title('A');        
        
        LL_Y;
        LL_Y_end;
        LL = LL + LL_Y + LL_Y_end;
    else % 'Normal' HMP
        for i=1:n
            local_prob = 0;
            for j=1:HMP.y_dim
                local_prob = local_prob + HMP.N(X(i),j) * exp( -(Y(i) - HMP.MU(j))^2 / (2*HMP.SIGMA(j))^2 )  / sqrt(2*pi*HMP.SIGMA(j));
            end
            LL = LL + log(local_prob);
        end
    end % if special
end  % if discrete
