% Convert number to string with special attention to x,y and m chroms
% Input:
% chrom - number of chromosome
% organism_str - organism (optional. Default is human)
% add_chr_flag - add the string 'chr' before
%
% Output:
% chr_str - string with the chromosome name
%
function chr_str = chr_num2str(chrom, organism_str, add_chr_flag, varargin)

if(~exist('organism_str', 'var') || isempty(organism_str))  % no input organism means human
    x_ind = 23;
else
    organism_str_tmp = genome_ver_to_organism(organism_str);
    Assign24MammalsGlobalConstants; organism_str = organism_str_tmp; % saved the input organism string
    organism_ind = -1;
    for i_org_str=1:length(genome_versions_vec)
        if(strcmp(organism_str, genome_versions_vec{i_org_str}{1}))
            organism_ind = i_org_str;
            break;
        end
    end
    organism_chr = genome_versions_vec{organism_ind}{3};
    x_ind = organism_chr;
end
y_ind=x_ind+1; m_ind=x_ind+2; % Y and Mitochondria

if(~exist('add_chr_flag', 'var') || isempty(add_chr_flag))
    add_chr_flag = 0;
end
n = length(chrom);
if(n > 1)  % New implementation: improving efficiency for large vectors
    chr_str = num2str_cell(chrom); % first move everything to str
    x_inds = vec2row(find(chrom == x_ind)); % deal only with sex chromosomes seperately
    y_inds = vec2row(find(chrom == y_ind)); 
    m_inds = vec2row(find(chrom == m_ind)); 
    for i=1:length(x_inds)
        chr_str{x_inds(i)} = 'X';
    end
    for i=1:length(y_inds)
        chr_str{y_inds(i)} = 'Y';
    end
    for i=1:length(m_inds)
        chr_str{m_inds(i)} = 'M';
    end
    if(add_chr_flag)
        for i=1:n
            chr_str{i} = ['chr' chr_str{i}];
        end
    end
else % here do just one
    chr_str = num2str(chrom);
    switch chrom
        case x_ind
            chr_str = 'X';
        case y_ind
            chr_str = 'Y';
        case m_ind
            chr_str = 'M';
    end
    if(add_chr_flag)
        chr_str = ['chr' chr_str];
    end
    
end
