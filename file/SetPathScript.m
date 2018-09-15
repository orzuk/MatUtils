% Add variuos directories to Matlab path
if(~exist('SetPathFlag', 'var')) % set a flag to make sure this script is run only once
    SetPathFlag = 0;
else
    SetPathFlag = 1;
end

if(~SetPathFlag)     % Add all function paths to matlab path
    AssignGeneralConstants;
    %    Assign24MammalsGlobalConstants;
    
    % Add many programs to matlab path - works for PC for now
    if(isempty(strmatch('/seq/orzuk/24mammals/src', pwd)) && (~exist('/cs/cbio/orzuk/software/Code/', 'dir'))) % home PC
        if(isempty(strmatch('\\oxygen\seq_orzuk\24mammals\src', pwd)) && ...
                isempty(strmatch('T:\24mammals\src', pwd))) % home
            machine = PC;
            %                matlab_libs_root_dir = 'C:\s\Or\Research\matlab\libs'; % old desktop
            cur_dir = pwd;
            my_path(cur_dir); % 'c:\research\24mammals\src'); % add mammals src
%            matlab_libs_root_dir = strrep(cur_dir, ...
%                '24mammals\src', 'matlab\libs'); % new laptop
            %                    matlab_libs_root_dir = 'C:\research\matlab\libs'; % new laptop
            my_path('C:/research/src/cvx/'); % NEW! Add cvx !!!
        else % pc broad
            machine = PC;
            html_outdir = 'Y:\public_html\data\';
%            matlab_libs_root_dir = 'Y:\public_html\matlab\libs';
            my_path('\\oxygen\seq_orzuk\24mammals\src'); % add mammals src
        end
        suite_sparse_dir = 'c:\research\src\SuiteSparse\';
    else % unix 
        machine = UNIX;
        html_outdir = '~orzuk/public_html/data/';
        cur_dir = pwd;
        my_path(cur_dir); % 'c:\research\24mammals\src'); % add mammals src
        if(exist('/cs/cbio/orzuk/software/Code/', 'dir')) % unix huji
            matlab_libs_root_dir = '/cs/cbio/orzuk/software/Code/Matlab'; % should be assigned in AssignGeneralConstants
            suite_sparse_dir = '/cs/cbio/orzuk/software/Code/Matlab'; % Temp, empty !!! 
        else % unix broad
            matlab_libs_root_dir = '~orzuk/public_html/matlab/libs';
            %        my_path('\\oxygen\seq_orzuk\24mammals\src'); % add mammals src
            suite_sparse_dir = '/seq/orzuk/src/SuiteSparse/';            
            my_path('/seq/orzuk2/stephens_course/new_HMT/hmt/code/Matlab/'); % add WHMT project
        end
    end % if unix
    
    
    %    my_path('Y:\public_html\matlab\MoG\new'); % add MoG package
    
    
    
    %    my_path( '/seq/orzuk/zuk_weizmann/research/snp_gui/src'); % snp_gui funcs
    %    my_path('T:\zuk_weizmann\research\snp_gui\src'); % snp_gui funcs
    
    % path(path, 'MoG'); % other MoG location
    
    % matlab_libs = {'Mog', 'str', 'stats', 'math', 'bioinfo', 'utils', '};
    % matlab_libs = GetFileNames(matlab_libs_root_dir);
    
    matlab_libs = get_subdirectories(matlab_libs_root_dir);
    suite_sparse_libs = get_subdirectories(suite_sparse_dir);
    
    for i=1:length(matlab_libs)
        %        add_to_path_lib = matlab_libs{i};
        my_path(matlab_libs{i});  %     my_path(fullfile(matlab_libs_root_dir, matlab_libs{i}));
    end
    for i=1:length(suite_sparse_libs)
        %        add_to_path_lib = suite_sparse_libs{i}
        my_path(suite_sparse_libs{i});
    end
    add_to_path_lib = suite_sparse_dir
    my_path(suite_sparse_dir); % add also root suite_spares directory
    my_path(fullfile( dir_from_file_name(matlab_libs_root_dir), 'bugs'))
    
    addutils; % add also Peter J. Acklam's matlab utils to the path
else
    sprintf('Path already set.')
end

% Add from other projects (need to put it in the matlab libs part !!! just don't add it to html!!! )
my_path('../../cs_snps/src')

if(machine == PC)
    github_dir = fullfile('c:\Users\', user_str, 'Documents\GitHub');
else % on huji unix 
    github_dir = '/cs/cbio/orzuk/software/GitHub';
end    
    
github_repos = {'MatUtils', 'bioinformatics-utilities', 'graphviz4matlab', 'from_others', 'HMP', 'SVLS'}; % Do NOT Include following: COMPASS, hmt, SPARCWave, ComSeq, BCS
for j=1:length(github_repos)
    my_path(fullfile(github_dir, github_repos{j}));
    github_subdirs = get_subdirectories(fullfile(github_dir, github_repos{j}), 3);
    for i=1:length(github_subdirs)
        my_path(github_subdirs{i});
    end
end
% New: add graphviz:
my_path('C:\Program Files (x86)\graphviz-2.38\release\bin')

% New: don't add directories in dropbox (take from github path)
%matrix_completion_dir = 'C:\Users\user\Dropbox\matrix completion\src\my_thesis'; my_path(matrix_completion_dir);
%matrix_completion_subdirs = get_subdirectories(matrix_completion_dir, 3);
%for i=1:length(matrix_completion_subdirs)
%    if(~strmatch(remove_dir_from_file_name(matrix_completion_subdirs{i}), ))
%        my_path(matrix_completion_subdirs{i});
%    end
%end

% New: add intensix directories
intensix_dir = fullfile('C:\Users\', user_str, 'Google Drive\Intesnsix\src\Matlab');
my_path(intensix_dir); % add to path

% add path for courses directories
my_path(fullfile('C:\Users\', user_str, 'Google Drive\HUJI\Teaching\Regression_52320\Code')); 

% add path for wavelets+scattering directories
wavelets_libs = get_subdirectories(fullfile('C:\Users\', user_str, 'Dropbox\wavelet\src'));
for i=1:length(wavelets_libs)
    my_path(wavelets_libs{i});  %     my_path(fullfile(matlab_libs_root_dir, matlab_libs{i}));
end
%my_path('C:\Users\user\Dropbox\wavelet\src\MATLAB_src\')

dropbox_dir = fullfile('C:\Users\', user_str, 'Dropbox'); 
google_drive_dir = fullfile('C:\Users\', user_str, 'Google Drive'); 
