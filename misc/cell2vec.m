% Concatenate all cell array vectors into one big vector
%
% NEW: We also try to let a matrix be concatenated (provided all cells are
% of the same width!!)
%
% Input:
% c - a cell array
% delim - (optional) a delimiter between two consecutive cell array (for strings)
% orientation - (optional) concatenate as rows or as columns
%
% Output:
% v - one vector
% lens - the lengths of the vectors in c (so we know where to start from each element)
%
function [v, lens] = cell2vec(c, delim, orientation, varargin)

if(~iscell(c))
    v = c;
    lens = length(c);
    return; 
end    

n = length(c);
if(isa(c{1}, 'char')) % strings
    v = char(1);
    if(exist('delim', 'var') && (~isempty(delim)))
        for i=1:n-1 % don't put delimiter on the last cell
            c{i} = [c{i} delim];
        end
    end
    lens = length_cell(c'); % here lens includes delimiter!!
else % numeric. enable also matrix (not only vec!!!)
    size_cell = cell2mat(vec2column(cellfun(@size, c, 'uniformoutput', false)));
    w = unique(size_cell(:,1));
    l = unique(size_cell(:,2));
    if(~exist('orientation', 'var'))
        do_rows = (length(setdiff(w,0)) <= 1); % here we go all the way with rows (otherwise columns)
    else
        do_rows = orientation;
    end
    if(do_rows)
        lens = size_cell(:,2); width = w;
        t_inds = find(size_cell(:,1) == 0);
    else
        lens = size_cell(:,1); width = l;
        t_inds = find(size_cell(:,2) == 0);
    end
    c(t_inds) = cellfun(@transpose, c(t_inds), 'uniformoutput', false);
    if(~isempty(t_inds))
        size_cell = cell2mat(vec2column(cellfun(@size, c, 'uniformoutput', false)));
        w = unique(size_cell(:,1));
        l = unique(size_cell(:,2));
        if(do_rows)
            lens = size_cell(:,2); width = max(w);
        else
            lens = size_cell(:,1); width = max(l);
        end
    end
    if( (min(length(w),length(l)) > 1) && (min(max(w),max(l)) > 1) )
        sprintf('Error! Problem with merging: elements not of same width/length!')
        return;
    end
    
    switch class(c{1})
        case 'single'
            v = zeros(width, sum(lens), 'single');
        case 'vpi' % variable length integers
            v = vpi(zeros(width, sum(lens)));
        otherwise
            v = zeros(width, sum(lens));
    end
    if(~do_rows)
        v=v';
    end
end
ctr=0;

if( (isnumeric(c{1}) || isa(c{1}, 'vpi')) && (~do_rows) ) % everything's numeric AND columns
    for i=1:n
        if(~isempty(c{i}))
            v(ctr+1:ctr+lens(i),:) = c{i};
        end
        ctr = ctr+lens(i);
    end
else % strings OR numeric rows
    for i=1:n
        if(~isempty(c{i}))
            v(:,ctr+1:ctr+lens(i)) = c{i};
        end
        ctr = ctr+lens(i);
    end
end
lens = cumsum(lens); % take cumulative sum so we know where to start each vector
if(isrow(c))
    lens = lens';
end

%if(iscolumn(c{1})) % output a column vector
%    v=v';
%end
