% Just take some binding sites and save them nicely as a .txt file which
% can be read in excel. Note that we have to deal here with two cases:
% one in which we only have one pwm and the second in which we have
% multiple pwms and then BS_seqs is a cell array
%
% Input:
% SEQS - structure holding the sequence information
% BS_regions - regions indexes
% BS_positions - positions within regions
% BS_scores - binding sites scores
% BS_seqs - sequences at binding sites
% BS_strand - binding sites strand
% BS_loglike - conservation log-likelihood measure
% pwms - pwms names
% bs_outfile - where to save the stuff
% data - any additional data to be saved (can be of any number of columns)
% field_names - names of fields in variable data
%
function Dummy = save_bs_mat_as_text(SEQS, BS_regions, BS_positions, BS_scores, ...
    BS_seqs, BS_strand, BS_loglike, pwms, bs_outfile, data, field_names, varargin)

if(~exist('BS_loglike', 'var'))
    BS_loglike = [];
end
if(iscell(BS_scores))  % get number of TFs
    TFs = length(BS_scores);
    n_vec = length_cell(BS_regions); % allow different lengths for different TFs
    n = max(n_vec);
else
    TFs = size(BS_scores, 2);
    n = size(BS_regions,1);
    n_vec = repmat(n, TFs, 1); % all lengths are the same
end
m = 10; % take maximum m (later we'll reduce it)
R = cell(sum(n_vec)+1,m);
R{1,1} = 'chr'; R{1,2} = 'start'; R{1,3} = 'end'; R{1,4} = 'BS-position';
R{1,5} = 'strand'; R{1,6} = 'seq'; R{1,7} = 'BS-score';
m = 7;
if(isfield(SEQS, 'loglike'))
    if(~isempty(SEQS.loglike))
        m=m+1;
        R{1,m} = 'pi-LODs-score';
    end
end
if(isfield(SEQS, 'score_vec'))
    m=m+1;
    if(isfield(SEQS, 'score_str'))
        R{1,m} = SEQS.score_str;
    else
        R{1,m} = 'p53-induced'; % default! non-generic!!!
    end
end
if(TFs > 1)
    m=m+1;
    R{1,m} = 'pwm';
end
R = R(:,1:m); % reduce R to include only relevant data slots

for t=1:TFs
    L = size(pwms{t,2}, 2);
    if(iscell(BS_seqs))
        BS_seqs_nt = int2nt(unpack_seqs(BS_seqs{t}, L));
    else
        BS_seqs_nt = int2nt(unpack_seqs(BS_seqs, L));
    end
    if(iscell(BS_regions))
        if(isfield(SEQS, 'pos_start_vec'))
            BS_genomic_pos = vec2column(SEQS.pos_start_vec(BS_regions{t})) + BS_positions{t};
        else
            BS_genomic_pos = BS_positions{t};
        end
        if(isfield(SEQS, 'chr_vec'))
            BS_chr_vec = SEQS.chr_vec(BS_regions{t}); % New! removed the transpose from here
        else
            BS_chr_vec = repmat('-', length(BS_regions{t}), 1);
        end
    else
        if(isfield(SEQS, 'pos_start_vec'))
            BS_genomic_pos = vec2column(SEQS.pos_start_vec(BS_regions(:,t))) + BS_positions(:,t);
        else
            BS_genomic_pos = BS_positions(:,t);
        end
        if(isfield(SEQS, 'chr'))
            BS_chr_vec = SEQS.chr(BS_regions(:,t)); % New! removed the transpose from here
        else
            BS_chr_vec = repmat('-', size(BS_regions,1), 1);
        end
    end
    if(iscell(BS_strand))
        cur_BS_strand = BS_strand{t};
    else
        cur_BS_strand = BS_strand(:,t);
    end
    for i=1:n_vec(t)
        if(cur_BS_strand(i) == 0)
            BS_seqs_nt(i,:) = seqrcomplement(BS_seqs_nt(i,:));
        end
    end
    if(t == 1)
        cur_start_ind = 1;
    else
        cur_start_ind = sum(n_vec(1:t-1)) + 1;
    end
    for i=1:n_vec(t)
        j = i+cur_start_ind; % index in R to store data
        R{j,1} = BS_chr_vec(i);
        if(iscell(BS_regions))
            if(isfield(SEQS, 'pos_start_vec'))
                R{j,2} = SEQS.pos_start_vec(BS_regions{t}(i));
                R{j,3} = SEQS.pos_end_vec(BS_regions{t}(i));
            else
                R{j,2} = '-';
                R{j,3} = '-';
            end
            R{j,5} = BS_strand{t}(i);
            R{j,7} = num2str(BS_scores{t}(i));
        else
            if(isfield(SEQS, 'pos_start_vec'))
                R{j,2} = SEQS.pos_start_vec(BS_regions(i,t));
                R{j,3} = SEQS.pos_end_vec(BS_regions(i,t));
            else
                R{j,2} = '-';
                R{j,3} = '-';
            end
            R{j,5} = BS_strand(i,t);
            R{j,7} = num2str(BS_scores(i,t));
        end
        R{j,4} = BS_genomic_pos(i);
        
        R{j,6} = BS_seqs_nt(i,:);
        m=7;
        if(isfield(SEQS, 'loglike'))
            if(~isempty(SEQS.loglike))
                m=m+1;
                if(~isempty(BS_loglike))
                    if(iscell(BS_loglike))
                        R{j,m} = BS_loglike{t}(i);
                    else
                        R{j,m} = BS_loglike(i,t);
                    end
                else
                    if(iscell(BS_regions))
                        R{j,m} = sum(SEQS.loglike{BS_regions{t}(i)}(BS_positions{t}(i):BS_positions{t}(i)+L-1));
                    else
                        R{j,m} = sum(SEQS.loglike{BS_regions(i,t)}(BS_positions(i,t):BS_positions(i,t)+L-1));
                    end
                end
            end
        end
        if(isfield(SEQS, 'score_vec'))
            m=m+1;
            if(iscell(BS_regions))
                R{j,m} = SEQS.score_vec(BS_regions{t}(i));
            else
                R{j,m} = SEQS.score_vec(BS_regions(i,t));
            end
        end
        if(TFs > 1)
            m=m+1;
            R{j,m} = pwms{t,1};
        end
    end
    if(exist('data', 'var')) % add aditional data
        K = length(field_names); % number of columns we need to add
        for i=m+1:m+K % set headers
            R{1,i} = field_names{i-m};
            if(iscell(data))
                for j=1:n_vec(t)
                    R{j+cur_start_ind,i} = data{j,i-m};
                end
            else
                data_size = size(data)
                R_size = size(R)
                n_vec_is = n_vec(t)
                start_ind = cur_start_ind
                for j=1:n_vec(t)
                    R{j+cur_start_ind,i} = data(j,i-m);
                end
            end
        end
        
    end
end % loop on TFs

savecellfile(R,  bs_outfile); % save best pwms matches

Dummy = 0;
