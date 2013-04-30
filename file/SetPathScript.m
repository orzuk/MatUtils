% Add variuos directories to Matlab path
if(~exist('SetPathFlag', 'var')) % set a flag to make sure this script is run only once
    SetPathFlag = 0;
else
    SetPathFlag = 1;
end

if(~SetPathFlag)     % Add all function paths to matlab path
    Assign24MammalsGlobalConstants;
    
    % Add many programs to matlab path - works for PC for now
    if(isempty(strmatch('/seq/orzuk/24mammals/src', pwd)))
        if(isempty(strmatch('\\oxygen\seq_orzuk\24mammals\src', pwd)) && ...
                isempty(strmatch('T:\24mammals\src', pwd))) % home
            machine = PC;
            %                matlab_libs_root_dir = 'C:\Users\Or\Research\matlab\libs'; % old desktop
            matlab_libs_root_dir = 'C:\research\matlab\libs'; % new latop
            my_path('c:\research\24mammals\src'); % add mammals src
        else % pc broad
            machine = PC;
            html_outdir = 'Y:\public_html\data\';
            matlab_libs_root_dir = 'Y:\public_html\matlab\libs';
            my_path('\\oxygen\seq_orzuk\24mammals\src'); % add mammals src
        end
    else % unix broad
        machine = UNIX;
        html_outdir = '~orzuk/public_html/data/';
        matlab_libs_root_dir = '~orzuk/public_html/matlab/libs';
        my_path('\\oxygen\seq_orzuk\24mammals\src'); % add mammals src
    end
    
    
    %    my_path('Y:\public_html\matlab\MoG\new'); % add MoG package
    
    
    
    %    my_path( '/seq/orzuk/zuk_weizmann/research/snp_gui/src'); % snp_gui funcs
    %    my_path('T:\zuk_weizmann\research\snp_gui\src'); % snp_gui funcs
    
    % path(path, 'MoG'); % other MoG location
    
    % matlab_libs = {'Mog', 'str', 'stats', 'math', 'bioinfo', 'utils', '};
    % matlab_libs = GetFileNames(matlab_libs_root_dir);
    
    matlab_libs = get_subdirectories(matlab_libs_root_dir);
    
    for i=1:length(matlab_libs)
        my_path(matlab_libs{i});  %     my_path(fullfile(matlab_libs_root_dir, matlab_libs{i}));
    end
    my_path(fullfile( dir_from_file_name(matlab_libs_root_dir), 'bugs'))
    
    addutils; % add also Peter J. Acklam's matlab utils to the path
else
    sprintf('Path already set.')
end

% Add from other projects (need to put it in the matlab libs part !!! just don't add it to html!!! )
my_path('../../cs_snps/src')
