% Joint density of effect size and selection coefficient
% Input: 
% s - vector of selection coefficients
% beta - vector of effect sizes
% N - effective population size 
% coupling_str - string representing the coupling density
% coupling_param - vector with parameters representing the density 
% 
function g = selection_effect_size_coupling_density(s, beta, N, coupling_str, coupling_param)

if(isvector(s))
    num_s = length(s); num_beta = length(beta);
    [s_mesh, beta_mesh] = meshgrid(s, beta); % create matrices
else
    num_s = size(s,2); num_beta = size(s,1);
end
coupling_density = zeros(num_s, num_beta);

S = 4.*N.*s; % compute big S


switch lower(coupling_str)
    case 'linear' % here no noise at all. Direct strong coupling
        slope = coupling_param;
        if(isvector(s))
            for i=1:num_s
                [~, j] = min(abs(beta - s(i)));
                coupling_density(i,j) = 1;
            end
        else
            s_of_beta = [beta -beta] .*slope;
%            s_of_beta2 = -beta.*slope;
            %            diag_inds = sub2ind([num_s num_beta], 1:num_s, I);
        end
    case 'sigmoid'
        slope = coupling_param;
        offset = 0.5; % slope = 30; 
        s_of_beta = -0.2+0.2.*[ 1 - 1 ./ (1+exp(-(((beta)-offset).*slope)))  1 - 1 ./ (1+exp(-(((-beta)-offset).*slope)))]; %  1 ./ (1+exp(abs((beta-0.4).*slope)))-0.5];
        figure; plot(beta(:,1), s_of_beta(:,1), '.');   xlabel('\beta'); ylabel('s'); title('Sigmoid Curve');
        
    case 'quadratic'
        
    case 'uniform' % no coupling
        coupling_density(:) = 1;
        
        
    case 'eyre-walker' % use Eyre-Walker coupling density 
        g = EyreWalker_dist_pdf(coupling_param.tau, coupling_param.sigma_epsilon, coupling_param.delta, ... % parameters relating s and beta
               coupling_param.theta, coupling_param.k, s, beta, link_function);         
        
        
end

if(max(coupling_density(:)) == 0)
    epsilon = 0.000001;
    [~, I] = min(abs(s_of_beta(:,1:num_s)-s));
    if(size(s_of_beta, 2) > num_s)
        [~, J] = min(abs(s_of_beta(:,num_s+1:end)-s)); % (end:-1:1,:)- s)); % here flip beta and minus-beta (assume symmetry with respect to beta)
        diag_inds = union( sub2ind([num_s num_beta], 1:num_s, I), sub2ind([num_s num_beta], 1:num_s, J) );
    else
        diag_inds = sub2ind([num_s num_beta], 1:num_s, I);
    end
    more_diag_inds = find(abs(s_of_beta(:,1:num_s)-s) < epsilon);
    diag_inds = union(diag_inds, more_diag_inds);
    
    coupling_density(diag_inds) = 1;
end

g = coupling_density'; % 2-d return density matrix 

% coupling_density = normalize_hist2d(s, beta, coupling_density); % Normalize density (causes some problems ...)
