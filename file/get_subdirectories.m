% Get the entire path under a given directory
% Currently this is done in a BFS manner (so we go over directories level by level)
% Input:
% root_dir - an input directory
% depth - (optional) maximal depth to go in the tree (e.g. 1 gives only the immediate sub-directories). Default is infinity (go down the tree all the way)
% search_mode - Breadth-First-Search (default, 0) or Depth-First-Search (1)
% permissions - (new) get only directories with certain permissions
%
% Output:
% f - a .cell file with all sub-directories
%
function f = get_subdirectories(root_dir, depth, search_mode, permissions, varargin)

BFS=0; DFS=1; 

if(~exist('root_dir', 'var'))
    root_dir = pwd;
end
if(~exist('depth', 'var') || isempty(depth)) % no real limit on depth
    depth = 999999;
end
if(~exist('search_mode', 'var') || isempty(search_mode)) % order of searched directories
    search_mode = BFS;
end
AssignGeneralConstants;
machine = get_machine_type(); 
switch machine
    case PC
        g = setdiff(GetFileNames( fullfile(root_dir, '*'), 1), GetFileNames(fullfile(root_dir, '..*'), 1));
        g = setdiff(g, fullfile(root_dir, '.'));
    case UNIX
        g = setdiff(GetFileNames( fullfile(root_dir, '*'), 1), GetFileNames(fullfile(root_dir, '.*'), 1));
end
if(exist('permissions', 'var') && (~isempty(permissions))) % Take only directories with certain permissions
    [file_names file_permissions] = GetFilePermissions(root_dir, 1); % see for which dirs do we have permission
%     file_names_are = file_names
%     file_permissions_are = file_permissions
    g = intersect(g, file_names(find(file_permissions(:,permissions))) );
else
    permissions = []; 
end


n = length(g); ctr=1;
f = cell(1,n);
for i=1:n
    if(exist(g{i}, 'dir'))
        f{ctr} = g{i};
        ctr = ctr+1;
    end
end
f = f(1:ctr-1);

g_ctr=2; % now counter is where to insert the new sub-directories
if(search_mode == DFS)
    tmp_f = f;
end
if(depth > 1) % loop on all sub-dirs
    for i=1:ctr-1
%        run_f = f{i}
        g = get_subdirectories(f{i}, depth-1, search_mode, permissions);        
        if(~isempty(g))
            if(search_mode == BFS) % default
                f = [f g];
            else % DFS
                tmp_f = insert_cell(tmp_f, g, g_ctr);
            end
        end
        g_ctr = g_ctr+1+length(g);
    end
end
if(search_mode == DFS)
    f = tmp_f;
end





