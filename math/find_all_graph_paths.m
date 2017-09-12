% Find all directed paths between any couple of vertices
%
% Input:
% A - graph adjancy matrix
% output_format - whether to output indices of vertices seperated by length
% (default), or binary vectors and then they're combined by length (BINARY)
% full_path_flag - 0: take all paths (including partial, default). 1: just from source to sink
% 
% Output:
% paths - structure of all paths. Can be either binary or indices 
% num_paths_pairwise - how many paths connect each pair of nodes
% num_paths_per_length - how many paths are there of each length
%
function [paths num_paths_pairwise num_paths_per_length] = ...
    find_all_graph_paths(A, output_format, full_path_flag, varargin)

if(~exist('output_format', 'var') || isempty(output_format))
    output_format = 'indices';
end
if(~exist('full_path_flag', 'var') || isempty(full_path_flag)) % default is to take all parial paths
    full_path_flag = 0; 
end
N = length(A);
% paths = cell(N);

%source_nodes = get_graph_sources(A);
% sink_nodes = get_graph_sinks(A);
% cur_level_nodes = source_nodes;
% for i=1:levels
%     for j=1:length(cur_level_nodes)
%
%
%     end
% end

paths = cell(N,1); % paths{i} contains all paths of length i
for i=1:N         % initilize paths
    paths{i} = cell(N);
    for j=1:N
        if(A(i,j))
            paths{1}{i,j} = [i j]; % zeros(1,N); paths{1}{i,j}(i) = 1;  paths{1}{i,j}(j) = 1;
        end
    end
end
init_paths = paths{1};
for i=2:N % loop up to N
    paths{i} = path_matrix_mult(paths{i-1}, init_paths);
end

if(full_path_flag) % take only paths from source to sink
    source_nodes = get_graph_sources(A);
    sink_nodes = get_graph_sinks(A);
    for i=1:N % loop on path length
        for j=1:N
            for k=1:N
                if(~any(j == source_nodes) || ~any(k == sink_nodes))
                    paths{i}{j,k} = [];
                end
            end
        end
    end
end

num_paths_per_length = zeros(N,1);
num_paths_pairwise = zeros(N);
for i=1:N
    cur_num_paths = cellfun('size', paths{i}, 1); % 'uniformoutput', false);
    num_paths_per_length(i) = sum(cur_num_paths(:));
    num_paths_pairwise = num_paths_pairwise + cur_num_paths;
end

switch output_format
    case {'binary', 'BINARY'}
        new_paths = cell(N); 
        for i=1:N % unite paths of different lengths - do we really need this?
            %   paths{1} = union_cell(paths{1}, paths{i}, 'rows');
            for j=1:N
                for k=1:N
                    tmp_paths = zeros(size(paths{i}{j,k}, 1), N); 
                    for r = 1:size(paths{i}{j,k}, 1)
                        tmp_paths(r,paths{i}{j,k}(r,:)) = 1;
                    end
                    new_paths{j,k} = [new_paths{j,k}' tmp_paths']'; 
                end
            end
        end% unite different lengths paths
        paths = new_paths;
end


% get a list of paths from two lists of paths
function C = path_matrix_mult(A, B)

N = length(A); C = cell(N);
for i=1:N
    for j=1:N
        C{i,j} = [];
        for k=1:N
            new_path = combine_paths(A{i,k}, B{k,j});
            if(~isempty(new_path))
                if(isempty(C{i,j}))
                    C{i,j} = new_path;
                else
                    %                    C{i,j} = union(C{i,j}, new_path, 'rows'); % old version
                    C{i,j} = [C{i,j}' new_path']'; % unique_cell([C{i,j}, new_path])
                end
            end
        end
    end
end


function C = combine_paths(A, B)

if(isempty(A) || isempty(B))
    C = [];
    return;
end
% Old version (keeps binary vectors - faster but doesn't save directions)
n1 = size(A, 1); n2 = size(B, 1);
% C = max(repmat(A, n2, 1), repmat(B, n1, 1));

C = [repmat(A(:,1:end-1), n2, 1) repmat(B, n1, 1)]; % new version: no need for cell array (all paths of same size)

% New version (keeps indices - slower but saves directions)
% C = cell(length(A)*length(B),1);
% for i=1:length(A)
%     for j=1:length(B)
%         C{(i-1)*length(B)+j} = [A{i} B{j}(2:end)];
%     end
% end
%C = unique(C); % probably not needed ...




