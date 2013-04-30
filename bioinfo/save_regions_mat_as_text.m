% Save a vector of regions as .txt file
%
% Input:
% outfile - where to save stuff
% chr_vec - chromosomes
% pos_start_vec - position starts
% pos_end_vec - position ends
% genome_version - what genome
% data_struct - additional data structure
% add_chr_flag - flag saying if to add the 'chr' flag before stuff
%
function R = save_regions_mat_as_text(outfile, chr_vec, pos_start_vec, pos_end_vec, ...
    genome_version, data_struct, add_chr_flag, varargin)

Assign24MammalsGlobalConstants();

if(~exist('add_chr_flag', 'var') || isempty(add_chr_flag))
    add_chr_flag = 0;
end

if(ischar(chr_vec)) % here accept file name as input 
    load(chr_vec);  % just load all data 
    if(exist('data', 'var'))
        data_struct = data; 
    end
end

n = length(chr_vec);
if(exist('data_struct', 'var') && (~isempty(data_struct)))
    m = 3+size(data_struct, 2);
else
    m = 3;
end
R = cell(n,m);

if(~isempty(chr_vec)) % check if we have something to save
    if(length(chr_vec) == 1)
        R{:,1} = chr_num2str(vec2column(chr_vec), genome_version);
    else
        R(:,1) = chr_num2str(vec2column(chr_vec), genome_version);
    end
    R(:,2) = num2str_cell(vec2column(pos_start_vec));
    R(:,3) = num2str_cell(vec2column(pos_end_vec));

    if(add_chr_flag) % change format to: chr start stop
        for i=1:n
            R{i,1} = ['chr' R{i,1}];
        end
    end
    if(exist('data_struct', 'var') && (~isempty(data_struct)))
        if(iscell(data_struct))
            for i=1:n
                if(mod(i, 1000) == 0)
                    sprintf('copying region %ld out of %ld', i, n)
                end
                for j=4:m
                    R{i,j} = data_struct{i,j-3};
                end
            end
        else % data is a matrix
            for i=1:n
                if(mod(i, 1000) == 0)
                    sprintf('copying region %ld out of %ld', i, n)
                end
                for j=4:m
                    R{i,j} = num2str(data_struct(i,j-3));
                end
            end
        end
    end
    if(add_chr_flag == 2) % format is: chrX:start-stop
        for i=1:n
            R{i,1} = [R{i,1} ':' R{i,2} '-' R{i,3}];
        end
        R = R(:,[1,4:m]); % remove column 2 and 3
    end
end % if empty

savecellfile(R, outfile);
