% Internal function for deciding when to stop trying (if number of failures
% is too large). Currently be conservative: stop if at least alpha*N
% failures (we can do better by getting confidence intervals on the
% observed alpha)
function break_flag = decide_success_or_failure_online( ...
    current_successes, current_failures, alpha, iters)

break_method = 'binomial_confidence_interval'; % if we're confidence number is out of range 
switch break_method
    case 'reach_alpha_failures'
        if(current_failures >= alpha*iters)
            break_flag = 1;
        else
            break_flag = 0;
        end
        
    case 'binomial_confidence_interval' % can be improved by lowering variance on the extremes
        N = current_successes + current_failures;
        p = current_failures / N; % estimator of frequency of errors
        beta = 0.001; % probability we're outside the confidence interval 
        confidence_sigma = norminv(1-beta/2) * sqrt(1/(4*N));
        confidence_interval = [p - confidence_sigma p + confidence_sigma]; 
        
        break_flag = 0; 
        if(alpha < confidence_interval(1)) % here we've failed 
            break_flag = 1;
        end
        if(alpha > confidence_interval(2)) % here we've succeeded
           break_flag = 2;  
        end
end
