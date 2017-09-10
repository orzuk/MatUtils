% Convert heritability to a differnet scale
% Use Falconer approximate formula This is for narrow-sense heritability.
% For broad sense we develop our own approximate formula
%
% Input:
% h_input - heritability on input scale
% h_output_scale - what scale to output heritability ('liability' or 'binary')
% mu - prevalence
%
% Output:
% h_output - heritability on output scale
% h_broad_output - broad sense heritability on output scale (binary only)
% h_output_exact - heritability on output scale, computed exactly from Gaussian integral (no approximation)
%
function [h_output h_broad_output h_output_exact] = ...
    heritability_scale_change(h_input, h_output_scale, mu)

threshold = norminv(1-mu);
z_score = normpdf(threshold); % Get Gaussian height at the incidence threshold
TOL = 0.00000000000000001;
num_inputs = max(length(threshold), length(h_input));
if(length(threshold) == 1)
    threshold = repmat(threshold, 1, num_inputs);
end
if(length(mu) == 1)
    mu = repmat(mu, 1, num_inputs);
end
if(length(h_input) == 1)
    h_input = repmat(h_input, 1, num_inputs);
end
if(iscolumn(h_input))
    mu = vec2column(mu);
end
switch h_output_scale
    case 'liability' % convert binary to liability
        h_output =  vec2row(h_input) .*  vec2row(mu .* (1-mu)) ./ vec2row(z_score.^2); % why always row vector?        
        if(nargout > 1) % heavy part. compute only if neccessary
            %h_broad_output = h_output; % this is h NARROW on liability scale assuming that the INPUT was broad-sense. This is computed by binary search
            h_min = 0; h_max = 1; h_mid = 0.5;
            while (h_max - h_min > TOL)
                [h_binary_output h_binary_broad_output] = ...
                    heritability_scale_change(h_mid, 'binary', mu);
                if(h_binary_broad_output > h_input)
                    h_max = h_mid;
                else
                    h_min = h_mid;
                end
                h_mid = (h_max + h_min)/2;
            end
            h_broad_output = h_mid;
            h_output_exact = h_output; % temp. (approximation)
        end
        if(iscolumn(h_input))
            h_output = vec2column(h_output);
            if(exist('h_broad_output', 'var'))
                h_broad_output = vec2column(h_broad_output);
                h_output_exact = vec2column(h_output_exact);
            end
        end
        
    case 'binary' % convert liability to binary
        h_output =  z_score.^2 .* h_input ./  (mu .* (1-mu)); % narrow-sense computed using Falconer's formula
        
        if(nargout > 1) % heavy part. compute only if neccessary
            h_output_exact = zeros(size(h_output));
            %        h_broad_output = zeros(size(h_output));
            triple_integral = zeros(size(h_output));
            
            for i=1:length(h_input) % Slow part: integral doesn't vectorize
                if(nargout > 2)
                    h_output_exact(i) = quadl(@(x)narrow_liability_conversion_fun(x, threshold(i), sqrt(h_input(i))), ...
                        -10*max(h_input(i), 1-h_input(i)), 10*max(h_input(i), 1-h_input(i)), TOL) .^ 2 .* ... % double integral. Take +/- 10 st.d.
                        1 ./  (mu(i) .* (1-mu(i))); % compute the exact integral
                    %            h_input(i) ./  (mu(i) .* (1-mu(i))); % compute the exact integral
                end
                triple_integral(i) = quadl(@(x)broad_liability_conversion_fun(x, threshold(i), sqrt(h_input(i))), ...
                    -10*max(h_input(i), 1-h_input(i)), 10*max(h_input(i), 1-h_input(i)), TOL); % triple integral. Take +/- 10 st.d.
                %            h_broad_output(i) = 1 - (mu(i) - triple_integral) / (mu(i)*(1-mu(i))); % get fraction of binary variance explained
            end
            h_broad_output = 1 - (mu - triple_integral) ./ (mu.*(1-mu));
        end        
    case 'plot-gaussian' % Here just plot some gaussian integrals
        mu = 0.0000001; x_mu = norminv(1-mu); c = 0.5;
        res = x_mu / 1000;
        x_vec = (-5*x_mu):res:(5*x_mu);
        y_vec = normpdf(x_vec);
        figure; hold on; plot(x_vec, y_vec);
        z_vec = normcdf((c*x_vec - x_mu) ./ sqrt(1-c^2));
        plot(x_vec, z_vec, 'r');
        w_vec = y_vec .* z_vec; w_vec = 0.1*w_vec ./ max(w_vec);
        plot(x_vec, w_vec, 'g', 'linewidth', 3);
        legend('\phi', '\Phi', 'mult.');
end



