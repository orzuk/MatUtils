% Read a wig formatted file into mat file.
% This should also work for bed files
%
% Input:
% wig_file - name of wig file
% genome_version (optional) - what genome version
% mat_outfile - file where to write output (.mat)
% numeric_flag - data is already numeric, faster (default - NO)
% field_names - NEW! decompose the variable 'data' into many fields
%
% Output:
% chr_vec - chromosomes
% pos_start_vec - start positions
% pos_end_vec - end positions
% data - additional data matrix
%
function [chr_vec pos_start_vec pos_end_vec data] = ...
    ReadWigFileToMat(wig_file, genome_version, mat_outfile, ...
    numeric_flag, field_names, varargin)

addutils;
Assign24MammalsGlobalConstants;

wig_suffix = suffix_from_file_name(wig_file);
if( (~strcmp(wig_suffix, 'txt')) && (~strcmp(wig_suffix, 'bed'))  )    % txt file. bed also works
    if(~exist([wig_file '.txt'], 'file'))
        %        system(['cp ' wig_file ' ' wig_file '.txt']);
        copyfile(wig_file, [wig_file '.txt']);
    end
    wig_file = [wig_file '.txt'];
end

if(exist('numeric_flag', 'var') & numeric_flag)
    data = load(wig_file);
    chr_vec = data(:,1);
    pos_start_vec = data(:,2);
    pos_end_vec = data(:,3);
    data = data(:,4:end); % take forth column and more if there is some data
else
    data = textread(wig_file, '%s', 'delimiter', '\n');
    
    n = length(data);
    chr_vec = zeros(n,1, 'single');
    pos_start_vec = zeros(n,1);
    pos_end_vec = zeros(n,1);
    
    for i=1:n % first find the starting index
        if(strmatch('chr', data{i}))
            start_ind = i;
            break;
        end
    end
    tab = sprintf('\t');
    if(isempty(data)) % empty file
        chr_vec = []; pos_start_vec = []; pos_end_vec = [];
    else % non-empty file
        if(~isempty(strfind(data{start_ind}, tab)))
            delim = tab;
        else
            delim = ' ';
        end
        tmp = strsplit(data{start_ind}, delim); m = length(tmp) - 3;
        for i=start_ind:n
            if(mod(i,1000) == 0)
                sprintf('reading region %d out of %d', i, n)
            end
            tmp = strsplit(data{i}, delim);
            chr_vec(i) = chr_str2num(tmp{1}(4:end), genome_version);
            pos_start_vec(i) = str2double(tmp{2});
            pos_end_vec(i) = str2double(tmp{3});
            for j=1:m
                data{i,j} = tmp{3+j};
            end
        end
        data = data(start_ind:end,:);
        if(min(isnumeric_cell(data(1,:))))
            data = cell2mat(num2str_cell(data));
        end
        
        chr_vec = chr_vec(start_ind:end);
        pos_start_vec = pos_start_vec(start_ind:end);
        pos_end_vec = pos_end_vec(start_ind:end);
    end % if empty file
end % if numeric_flag

if(~exist('mat_outfile', 'var'))
    mat_outfile = file_name_to_mat(wig_file);
end
if(~isempty(mat_outfile)) % giving an empty out-file-name means you don't want to save
    save(mat_outfile, 'chr_vec', 'pos_start_vec', 'pos_end_vec', 'data');
    
    if(exist('field_names', 'var')) % here convert data to many fields
        save(mat_outfile, 'chr_vec', 'pos_start_vec', 'pos_end_vec'); % get rid of 'data'
        for i=1:length(field_names)
            eval([field_names{i} ' = cell2mat(str2num_cell(data(:,' num2str(i) ')));']);
            eval(['save(''' mat_outfile ''',''' field_names{i} ''', ''-append'');']);
        end
    end
    
end

