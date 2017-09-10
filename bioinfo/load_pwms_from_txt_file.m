% Load a set of pwms in a format convenient for humans to write in txt.
%
% Input:
% pwms_file - a name of a text file to output
% pwms_output_mat_file - (optional) output file to save the pwms in .mat format
% fasta_format - which format should we use: 1 - fasta, 0 - not-fasta, -1 - pouya's format 
%
% Output:
% pwms - a structure of the pwms
%
function pwms = load_pwms_from_txt_file(pwms_file, pwms_output_mat_file, fasta_format, varargin)

q = loadcellfile(pwms_file); % load sequences
% load('tmp_pouya_pwms.mat'); 

if(~exist('fasta_format', 'var')) % determine format based on data 
    if(isnumeric(q{2})) % set fasta format according to file type (8 lines/5 lines per pwm)
        fasta_format = 1;
    else
        fasta_format = 0;
    end
end

switch fasta_format % set number of pwms
    case 1 % fasta file
        n = floor(size(q, 1) / 5);
    case 0 % not fasta
        n = floor(size(q, 1) / 8); 
    case -1 % pouya's format: we don't know here how many pwms since they're on rows rather than columns 
        if(size(q,2) == 1) % NEW format
            pwms_lines = strfind_cell(q, '>'); % assume all names are with > (like fasta)
        else
            pwms_lines = strfind_cell(q(:,2), 'mer'); % assume all names are with mer
        end
        n = length(pwms_lines); pwms_lines = [pwms_lines size(q,1)+1]; % add a dummy line for the last index 
    case -2 % Nir's format: NEW
        pwms_lines = 1:6:length(q); 
        n = floor(size(q, 1) / 6); 
end
pwms = cell(n,4);

ctr = 1;
for i=1:n
    switch fasta_format
        case 1  % fasta
            pwms{i,1} = strdiff(q{ctr}, '>'); % name
            pwms{i,3} = ''; % description
            pwms{i,4} = ''; % long-description
        case 0 % non-fasta
            pwms{i,1} = q{ctr};
            pwms{i,3} = q{ctr+1};
            pwms{i,4} = q{ctr+2};
            ctr = ctr+2;
        case -1 % pouya's format. No need for counter here 
            if(size(q,2) == 1)
                pwms{i,1} =  q{pwms_lines(i), 1}(2:end); % get name
                pwms{i,3} = '';
            else
                pwms{i,1} = q{pwms_lines(i), 2}; % get name
                pwms{i,3} = q{pwms_lines(i)+1, 1}; % get branch-length cutoff
            end
            pwms{i,4} = '';
        case -2 % Nir's format
            pwms{i,1} = q{pwms_lines(i), 1}; % get name 
            pwms{i,3} = q{pwms_lines(i), 2}; % get description (database: transfac/bulyk/etc.)
            pwms{i,4} = q{pwms_lines(i), 3}; % longer description
    end
    switch fasta_format
        case {0,1,-2}
            for j=1:4
                if(size(q, 2) == 1) % all were put in the same cell
                    pwms{i,2}(j,:) = q{ctr+j};
                else % different numbers in different cells 
                    L = find(1-isempty_cell(q(ctr+1,:)), 1, 'last');
                    for k=1:L % size(q,2)
                        pwms{i,2}(j,k) = q{ctr+j,k};
                    end
                end
            end
        case -1 % pouya's format
            L = pwms_lines(i+1) - pwms_lines(i) - 2;
            for k=1:L
                pwms{i,2}(:,k) = str2nums(q{pwms_lines(i)+k+1,1});
            end
    end
    
    switch fasta_format
        case 1
            ctr = ctr+5;
        case {0,-2}
            ctr = ctr+6;
        case -1
            ctr = ctr+L; 
    end
end
if(exist('pwms_output_mat_file', 'var') && (~isempty(pwms_output_mat_file)))  % save pwms in output file
    save(pwms_output_mat_file, 'pwms');
end

