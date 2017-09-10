% Parse the genomic regions file obtained from Mitch,
% and put the information in a .mat file
%
% Input:
% data_dir - directory of regions file
% regions_file - name of regions file
% organism_str - which organism are we dealing with (relevant to the X,Y,M chroms.). Optional (default is human)
% tf_names - not used now 
%
function Dummy = GenomicRegionsFiletoMat(data_dir, regions_file, ...
    organism_str, tf_names, varargin)

H = loadcellfile(fullfile(data_dir, regions_file)); % load the .txt file  - we cannot just use 'load'


if(nargin == 2)
    organism_str = 'HUMAN';
end
if(size(H,2) > 2) % here assume the format to be:  chr5 1341234 1341235
    if( (~isnumeric(H{end,3})) || (isempty(H{end,3})) ) % make sure at the last line we're doing ok
        H = H(1:end-1,:);
    end
    if(H{end,3} < H{end,2}) % make sure end position is larger than start position 
        H = H(1:end-1,:);
    end

    chr_vec = H(:,1);
    pos_start_vec = cell2mat(H(:,2));
    pos_end_vec = cell2mat(H(:,3));
    if(size(H,2) > 3) % look for another field. We assume that the next one is the scores/affinities
        for j=4:size(H,2)
            if( min(isnumeric_cell(H(:,j))) )
                scores_vec = cell2mat(H(:,j));
                data = H(:,[4:j-1,j+1:end]); % save everything else as non-numeric
                break;
            end
        end
    end
    new_chr_vec = zeros(length(chr_vec), 1);
    for i=1:length(chr_vec)
        if(~isnumeric(chr_vec{i}))
            if(~isempty(strmatch('chr', chr_vec{i})))
                chr_vec{i} = chr_vec{i}(4:end);
            end
            if(~isempty(strfind(chr_vec{i}, '_random')))
                chr_vec{i} = chr_vec{i}(1:end-7);
            end
            new_chr_vec(i) = chr_str2num(chr_vec{i}, organism_str);
        else
            new_chr_vec(i) = chr_vec{i};
        end
    end
    chr_vec = new_chr_vec; % copy the numeric stuff
else % here assume the format to be:  chr5:1341234-1341235
    n = size(H,1)
    chr_vec = zeros(n,1); pos_start_vec = zeros(n,1); pos_end_vec = zeros(n,1);
    if(size(H,2) == 2) % look for another field. We assume that the next one is the scores/affinities
        if( min(isnumeric_cell(H(:,2))) )
            scores_vec = cell2mat(H(:,2));
        end
    end

    for i=1:n
        if(mod(i,100) == 0)
            i_is = i
        end
        first_stop = strfind(H{i,1}, ':');
        second_stop = strfind(H{i,1}, '-');
        if(~isempty(strfind(H{i,1}, '_random'))) % get rid of random stuff
            chr_vec(i) = chr_str2num(H{i,1}(4:first_stop-8), organism_str );
        else
            chr_vec(i) = chr_str2num(H{i,1}(4:first_stop-1), organism_str );
        end
        pos_start_vec(i) = str2num(H{i,1}(first_stop+1:second_stop-1));
        pos_end_vec(i) = str2num(H{i,1}(second_stop+1:end));
    end
end
output_file = fullfile(data_dir, [regions_file(1:end-4) '.mat']);
save(output_file, 'chr_vec', 'pos_start_vec', 'pos_end_vec'); % save overall and chromosome-specific structure ...
if(exist('tf_names', 'var')) % save also output file
    save(output_file, 'tf_names', '-append');
end
if(exist('scores_vec', 'var')) % save also output file
    save(output_file, 'scores_vec', '-append');
end
if(exist('data', 'var'))
    save(output_file, 'data', '-append');
end


Dummy = 0;


