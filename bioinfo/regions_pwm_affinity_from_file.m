% Compute affinity of a DNA sequence to a pwm, with input from files.
% Currently we are not given the pwm itself but rather file names with the input 
%
% Input:
% regions_file - name of file containing the regions and the pwms
% pwms_file - name of file containing the pwms 
% pwm_ind - which pwm to use
% do_log - flag saying if to do log to the pwm
% strand - which strand to compute (default should be both)
% beta - the temperature. -1 means threshold some affinity, -2 means take best site
% in_matlab_flag - run everything directly through matlab or send to the farm
% shuffle_pwm_flag - optional flag saying if we want to randomly shuffle the pwm
% iters - optional flag saying if we want to run many iters (in case shuffle flag is on)
% compute_pval_flag - do we want to give also p-value
% good_inds - ??? (currently not used) 
% ready_affinities_file - load the affinities from a file (rather than computing them) - DOESN'T WORK YET !!!
% (alpha - top-fraction for conservation cutoff CANCELLED FOR NOW !!! (should be in regions file))
%
% Output:
% regions_affinity - a number representing the total affinity of the regions to the pwm
% regions_affinity_pval_vec - a pvalue for each pwm for a set of regions 
%
function [regions_affinity_mat regions_affinity_pval_vec] = ...
    regions_pwm_affinity_from_file(regions_file, pwms_file, pwm_ind, do_log, strand, beta, ...
    in_matlab_flag, shuffle_pwm_flag, iters, compute_pval_flag, good_inds, ready_affinities_file, varargin)

load(regions_file); 
Assign24MammalsGlobalConstants;
if(machine == PC) % on pc we can't use jobs
    in_matlab_flag = 1;
end

tic;
if(~exist('regions', 'var'))
    regions = load(regions_file); % load everything to struct 
end
if(~exist('pwms', 'var'))
    load(pwms_file);
end
if(~isfield(regions, 'packed_seqs'))
   if(~isfield(regions, 'seqs'))
%       chr_inds = find(regions.chr_vec == 13)
       regions = ExtractSeqsByPositions(regions.chr_vec, regions.pos_start_vec, regions.pos_end_vec, ...
           seqs_dir, genome_version, 1); % extract only sequences 
   end
    [regions.packed_seqs regions.seqs_lens] = pack_seqs(regions.seqs);
    save(regions_file, 'regions', 'genome_version'); % save in seqs format : overwrite !!!! 
end

num_pwms = size(pwms,1); num_regions = length(regions.packed_seqs);
if(~exist('pwm_ind', 'var'))
    pwm_ind = 1:num_pwms;
end
if(isempty(pwm_ind))
    pwm_ind = 1:num_pwms;
end
pwms = pwms(pwm_ind,:); num_pwms = size(pwms,1); % change pwms to include onlt the ones we called to
if(~exist('shuffle_pwm_flag', 'var'))  % check if the shuffling of pwms flag is on
    shuffle_pwm_flag = 0;
end
if(~exist('iters', 'var')) % default is one iteration
    iters = 1;
end
if(~exist('compute_pval_flag', 'var'))
    compute_pval_flag = 0;
end
if(~exist('good_inds', 'var'))
    good_inds = [];
end
if(~isfield(regions, 'conservation_cutoff'))
    regions.conservation_cutoff = 0.05;
end

if(compute_pval_flag) % here just compute the p-value and throw away the affinities
    regions_affinity_pval_vec = zeros(num_pwms,1);
end
label_str = remove_dir_from_file_name(regions_file(1:end-4)); % get just the file name without the dir and the .mat ending
if(in_matlab_flag) % run everything within matlab - on one pwm!!
    in_matlab = cputime
    temp_pwm_ind = 1:num_pwms;
    if(exist('ready_affinities_file', 'var'))
        load(ready_affinities_file);
        regions_affinity = regions_affinity(:,good_inds); 
    else
        regions_affinity = zeros(iters*num_pwms, num_regions, 'single'); % prepare space for regions affinity matrix (single to save memory)
    end
    for p=1:num_pwms
        doing_pwm = p
        for iter = 1:iters % loop on number of iterations
            if(mod(iter, 100) == 0)
                sprintf('doing iter %d out of %d\n', iter, iters)
            end
            if(shuffle_pwm_flag  && (iter > 1)) % new! no shuffling on first iteration!
                L = size(pwms{temp_pwm_ind(p),2},2); P = randperm(L); pwms{temp_pwm_ind(p),2} = pwms{temp_pwm_ind(p),2}(:,P);  % shuffle pwm
            end
            if(~exist('ready_affinities_file', 'var'))
                if(~isfield(regions, 'loglike'))   % no conservation
                    regions_affinity(iter+(p-1)*iters,:) = regions_pwm_affinity(pwms{temp_pwm_ind(p),2}, regions.packed_seqs, ...
                        regions.seqs_lens, do_log, strand, beta);
                else % use conservation
                    % prepare the loglikelihood of the regions (need to sum over the motif length)
                    loglike_mat = cell(length(regions.loglike),1);
                    for i=1:length(regions.loglike)
                        if(length(regions.loglike{i}) >= size(pwms{temp_pwm_ind(p),2}, 2))
                            loglike_mat{i} = my_smooth(regions.loglike{i}, size(pwms{temp_pwm_ind(p),2}, 2));
                        end
                    end
                    regions_affinity(iter+(p-1)*iters,:) = regions_pwm_affinity(pwms{temp_pwm_ind(p),2}, regions.packed_seqs, ...
                        regions.seqs_lens, do_log, strand, beta, loglike_mat, regions.conservation_cutoff);
                end
            end
        end % loop on number of iterations
        if(compute_pval_flag) % here just compute the p-value and throw away the affinities
            true_affinity_vec = regions_affinity(1+(p-1)*iters,:); % the first row is without shuffling
            good_inds = [];
            regions_affinity_pval = regions_affinity_mat_to_pval(true_affinity_vec, ...
                good_inds, regions_affinity(1+(p-1)*iters:p*iters,:));
            regions_affinity_pval_vec(p) = regions_affinity_pval;
        end
    end   % loop on pwms
    save(fullfile('temp', ['affinity_pwm_' num2str(pwm_ind(1)) '_' label_str '.mat']), 'regions_affinity'); % here the matrix is saved !!!
    regions_affinity_mat = regions_affinity;
    if(compute_pval_flag)
        save(fullfile('temp', ['affinity_pwm_' num2str(pwm_ind(1)) '_' label_str '.mat']), ...
            'regions_affinity_pval', '-append'); % here the matrix is saved !!!
        regions_affinity_mat = []; % we don't need the affinities ... 
    end
    if(in_matlab_flag == 2) % a code for removing the big matrix
        regions_affinity_mat = [];
    end
else % send a job to the farm - enumerate over pwms (ignore pwm_ind)
    for i=1:num_pwms % first submit all jobs. We loop on the pwms since we can't store all scores at one time
        if(shuffle_pwm_flag)
            L = size(pwms{i,2},2); P = randperm(L); pwms{i,2} = pwms{i,2}(:,P);  % shuffle pwm
        end
        if(exist(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat']), 'file'))
            delete(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat'])); % remove any old files with the same name
        end
        % Important! no spaces in job_str
        %%% regions_affinity = 
        job_str = ['regions_affinity=regions_pwm_affinity_from_file(''' regions_file ''',''' ...
            pwms_file ''',' num2str(i) ',' ...
            num2str(do_log) ','  num2str(strand) ',' num2str(beta) ',2,' ...
            num2str(shuffle_pwm_flag) ',' num2str(iters) ',' ...
            num2str(compute_pval_flag) ');']; % run it through the farm, lightening the within matlab flag
        SubmitMatlabJobToFarm(job_str, fullfile('temp', ['job_output_affinity_pwm_' num2str(i) '_' label_str '.out']), ...
            'short', []); % for now we do not use the bioinformatics toolbox
    end
    if(~compute_pval_flag)
        regions_affinity_mat = zeros(num_pwms*iters, num_regions, 'single'); % new! we take into account also many iterations
    else
        regions_affinity_mat = zeros(num_pwms, 1, 'single');
    end
    sec_ctr = 0;
    finished_jobs_vec = zeros(1,num_pwms); % indicator variables saying which jobs had finished already
    while(sum(finished_jobs_vec) < num_pwms) % loop until all jobs are successfully finished
        unfinished_inds = find(finished_jobs_vec == 0); % all unfinished jobs. Must be a row vector
        for i=unfinished_inds % now collect all results. Note: we need to account for the fact that sometimes the files a job fails to run due to license problems
            temp_file_name = fullfile('temp', ['job_output_affinity_pwm_' num2str(i) '_' label_str '.out']);
            if(exist(temp_file_name, 'file')) %check if the job output log file is created (don't wait!)
                text_vec = textread(temp_file_name, '%s', 'delimiter', '\n');
                if( (~isempty(strmatch('License Manager Error', text_vec))) || ...
                        (~isempty(strmatch('Error in', text_vec))) || ...
                        (~isempty(strfind_cell(text_vec, 'syntax error'))) || ...
                        (~isempty(strfind_cell(text_vec, 'job killed'))) ) % this means that we need to run the job again
                    %%% regions_affinity = 
                    job_str = ['regions_affinity=regions_pwm_affinity_from_file(''' regions_file ...
                        ''',''' pwms_file ''',' num2str(i) ',' ...
                        num2str(do_log) ','  num2str(strand) ',' num2str(beta) ',2,' ...
                        num2str(shuffle_pwm_flag) ',' num2str(iters) ...
                        ',' num2str(compute_pval_flag) ');']; % run it through the farm, lighting the within matlab flag
                    SubmitMatlabJobToFarm(job_str, fullfile('temp', ['job_output_affinity_pwm_' num2str(i) '_' label_str '.out']), ...
                        'short', []); % for now we do not use the bioinformatics toolbox
                    sprintf('License Problem! Need to Submit Job for PWM %d again!', i)
                    text_vec_is = text_vec 
                    license_error = strmatch('License Manager Error', text_vec)
                    other_errors = strmatch('Error in', text_vec)
                else % this means that the job was probably successful
                    while(~exist(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat']), 'file')) % wait until the file is created ...
                        if(mod(sec_ctr, 100) == 0)
                            looking_for_file = fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat'])
                            secs_passed = sec_ctr
                        end
                        sec_ctr = sec_ctr + 0.1;
                        pause(0.1); % check every one second (don't want to check to often)
                    end
                    clear regions_affinity;
                    while(~exist('regions_affinity', 'var'))
                        load(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat']));
                        pause(0.1);
                    end
                    if(~compute_pval_flag)
                        regions_affinity_mat((i-1)*iters+1:i*iters,:) = regions_affinity;
                    else
                        %                        regions_affinity_pval_vec(i) = regions_affinity_pval; % this should be in the same file as regions_affinity
                        regions_affinity_mat(i) = regions_affinity_pval; % this should be in the same file as regions_affinity
                    end
                    finished_jobs_vec(i) = 1; % this job was finished successfully
                    if(exist(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat']), 'file'))
                        delete(fullfile('temp', ['affinity_pwm_' num2str(i) '_' label_str '.mat']));  % delete the temporary file
                    end
                    sprintf('Done job for pwm #%d. Finished %d out of %d pwms\n', i, sum(finished_jobs_vec), num_pwms)
                end
            end
        end
        pause(0.1); % avoid too rapid checkings of job status
    end
end
toc
