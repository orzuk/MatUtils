% Make automatically an html file telling where exactly are all functions.
% This will still be needed to be edited manually
%
% Input:
% matlab_libs_dir - directory with all matlab functions
% matlab_html_template_file - a template file to generate html from
%
function GenerateLibHTML(matlab_libs_dir, matlab_html_template_file)
addutils;

Assign24MammalsGlobalConstants;
if(machine == PC)
    delim = '\';
else
    delim = '/';
end
if(matlab_libs_dir(end) ~= delim)
    matlab_libs_dir = [matlab_libs_dir delim];
end

comment_convertion_struct = { ... % a vector with comments for each sub-directory
    'BCS', 'Utilities for the BCS algorithm.', ...
    ['  <li> <b> Overview: </b> ', ...
    '<br> This page accompanies the paper:' ...
    '<br> <A &nbsp;  href="http://www.broadinstitute.org/~orzuk/publications/AmirZuk_BCS_JCB.pdf">Bacterial Community Reconstruction Using Compressed Sensing</A>, A. Amir and O. Zuk, Journal of Computational Biology, 18(11):1-20 (2011) <br>', ... %    'The online supplamentary information can be viewed here <A &nbsp;  href="http://www.broadinstitute.org/~orzuk/matlab/libs/BCS/BCS_Recomb_supp_info.pdf"> here </A>.<br>', ...
    'The paper presents the Bacterial Compressed-Sensing (BCS) algorithm. The source code is available as a Matlab package below. <br>' ...
    'The package contains functions for mixture chromatogram simulations, pre-processing and reconstruction.', ...
    '<br> Additional functions are used for handling the database of s16-RNA gene, and for displaying results. <br>', ...
    'Please see the file "readme.txt" for a coprehensive description of the package. <br>', ...
    'In addition, a few data files related to the project are available for download <A &nbsp;  href="http://www.broadinstitute.org/~orzuk/matlab/libs/BCS/BCS_data.zip"> here </A>.', ... 
    ' <br> <li> <b> About: </b> ' ...
    ' <br> BCS is a collaboration between the Broad Institute and the Weizmann Institute of Science.', ...
    ' <br> <li> <b> Installation: </b> ', ...
    ' <br> To install the package simply download the file BCS.tgz to a directory of your choice, extract it, and add to your Matlab path.' ...
    ' <br> The package uses the GPSR algorithnm, which is available through <A &nbsp;  href="http://www.lx.it.pt/~mtf/GPSR/"> here </A>.' ...      
    ' <br> For questions, suggestions, bug-reports etc. please contact Amnon Amir (<A href="mailto:amnon.amir@weizmann.ac.il">amnon.amir@weizmann.ac.il</A>) or ' ...
    'Or Zuk (<A href="mailto:orzuk@broadinstitute.org">orzuk@broadinstitute.org</A>). <br>', ...
    'Below we list all functions available in the package with their description: ']; ...
    'examples', 'Example utilities', ''; ...
    'str', 'String utilities', ''; ...
    'bioinfo', 'Bioinformatics utilities', ''; ...
    'cs_snps', 'Utilities used for compressed-sensing pooling design', ...    
    ['  <li> <b> Overview: </b> ', ...
    ' <br> cs_snps is a set of utilities used for simulating and analyzing compressed-sensing pooling design',  ...
    ' <br> experiments for next-generation sequencing projects whose purpose is to detect rare allele carriers in a population', ...
    ' <br> These are used in the paper: ', ...
    ' <br> <i> Rare-Allele Detection Using Compressed Se(que)nsing, N. Shental, A. Amir and O. Zuk (2009) </i>', ... %    ' <br> avaiable currently on the <A &nbsp;  href="http://arxiv.org/abs/0909.0400"> arxiv </A>.', ... 
    ' <br> The package enables an in-silico simulation of a sequencing pooling experiments, including', ...
    ' <br> modeling of sequencing error parameters. The package also enables reconstructing genotypes of', ... 
    ' <br> individuals profiled in the pools. Detailed instructions and usage examples are given in the README file', ...
    ' <br> <li> <b> About: </b> ' ...
    ' <br> cs_snps is a collaboration between the Broad Institute, the <A &nbsp;  href="http://www.openu.ac.il/home/shental/ "> Open University of Israel </A>, and the Weizmann Institute ' ...
    ' <br> of Science.', ...
    ' <br> <li> <b> Installation: </b> ', ...
    ' <br> To install the package simply download the file cs_snps.tgz to a directory of your choice, extract it, and add to your Matlab path.' ...
    ' <br> The package uses the GPSR algorithnm, which is available through <A &nbsp;  href="http://www.lx.it.pt/~mtf/GPSR/"> here </A>.' ...      
    ' <br> List of functions available with their description: ']; ...
    'misc', 'Misceleneous utilities', ''; ...
    'stats', 'Statistics utilities', ''; ...
    'math', 'Mathematical utilities', ''; ...
    'file', 'Files utilities', ''; ...
    'plot', 'Plotting and Imaging utilities', ''; ...
    'from_others', 'Utilities written by other people which I use ', ''; ...
    'pac_perm', 'Utilities used for ranking stability project.',  ...
    ['They were used for the following papers:  ' ...
    ' <br> Thousands of Samples are Needed to Generate a Robust Gene List for Predicting Outcome in Cancer, ' ...
    'L. Ein-Dor, O. Zuk and E. Domany. PNAS 103(15):5923-5928 (2006)', ...
    ' <br> Ranking Under Uncertainty, O. Zuk, L. Ein-Dor and E. Domany. UAI 2007: 466-473']; ...
    'ipod', 'Utilities used for percolation-on-networks project  ', ''; ...
    'aging', 'Utilities used for Aging gene-expression project  ', ''; ...
    'snp_gui', 'Utilities used for snp-gui application. ', ''; ...
    'mog', 'Utilities for Mixture-of-Gaussians package ', ''; ...
    'fdr', 'Utilities for False-Discovery-Rate package.',  ...
    ['They were used in the following paper: ' ...
    ' <br> <i>FDR Control with adaptive procedures and FDR monotonicity (2009), A. Zeisel, O. Zuk and E. Domany.</i>' ...
    ' Please cite this paper if you use them', '']; ...
    'hmp', 'Utilities for Hidden-Markov-Processes package ', ''; ...
    'bnt', 'Utilities for Bayesian Networks package (auxillary functions for Kevin Murphy''s BNT) ', ''; ...
    'bnt_dim', 'Utilities used for Bayesian Networks Dimension determining package ', ''; ...
    };

R = loadcellfile(matlab_html_template_file, [], 9); % load template file. Should be manyXone cell-array !!! 
start_ind = strfind_cell(lower(R), '<tr>');
R_start = R(1:start_ind(1)); R_end = R(start_ind(end):end);

title_line_inds = strfind_cell(R, 'Template');

% Prepare the master file
matlab_master_html_libs_file = fullfile(matlab_libs_dir, 'matlab_utils.html');
S = R_start;     ctr = length(S)+1;
S{ctr+1} = '<ul><li> This page contains a set of Matlab utilities which I and others have written for various projects. <br>';
S{ctr+2} = 'I have decided to put them here with the thought that they might be useful for other people. <br>';
S{ctr+3} = ['For the vast majority of functions you just need to download them and they will work fine, but a few functions <br>' ...
    'currently work only on the broad enviourment due to verious data requiements (these will be removed soon) <br> ' ...
    'The functions are supplied as is, with no warranty. You can use them freely - I would appreciate if you send me <br>' ...
    ' a note in case you use them so I can keep track of an estimated number of users. </li> </ul>'];
S{ctr+4} = ['<ul><li> I was motivated mainly by the excellent Matlab utilities page by ' ...
    '<a href="http://home.online.no/~pjacklam/matlab/software/util/index.html">Peter J. Acklam</a> which I found very useful. <br>'];
S{ctr+5} = ['Other great Matlab packages which I have used are in <a href="http://www.cs.ubc.ca/~murphyk/Software/index.html">Kevin P. Murphy''s</a>, ' ...
    ' <a href="http://wise-obs.tau.ac.il/~eran/matlab.html">Eran O. Ofek''s</a> and <br>' ...
    '<a href="http://carmelab.huji.ac.il/software.html#general"> Liran Carmel''s</a> websites, ' ...
    ' in <a href="http://stommel.tamu.edu/~baum/toolboxes.html"> S. Baum''s</a> index and of course in <a href="http://www.mathworks.com/matlabcentral/">Matlab Central</a>. <br> </li> </ul>'];
S{ctr+6} = '<ul><li> Many people have contributed to these functions. I apologize for any omissions of credits. In particular, <br>';
S{ctr+7} = '<b><i> Assif Yitzhaki, Libi Hertzberg, Liat Ein-Dor, Yuval Tabach, Tal Shay, Noam Shental, Amnon Amir </i></b>and <b><i>Amit Zeisel </i></b> <br> have written/co-written many functions. <br>';
S{ctr+8} = 'Many functions also call/are modifications of code from other people (see the ''from_others'' sub-directory). <br>';
S{ctr+9} = 'If you want your code to be added/removed/modified/credited, please let me know. <br>';
S{ctr+10} = 'Any other comments, bug reports, suggestions, etc. are mostly welcome. If you report a bug/have a question<br>';
S{ctr+11} = 'please be specific and describe the problem in detail (e.g. input, output, I expect XXX but get YYY etc.) <br> </li> </ul>';
S{ctr+12} = '<ul><li> Some functions require mex files generated by c code for increased speed. These appear <br>';
S{ctr+13} = 'in the directory of the relevant library as .mexw32 (windows) or .mexa64 (unix). If they don''t work <br>';
S{ctr+14} = ' for you, you will need to compile the sources yourself - these .cpp and .h files are also available in the same directories </li> </ul>';
S{ctr+15} = '<ul><li> Installation should be quite simple: (contact me if you have troubles installing) <br>';
S{ctr+16} = '&nbsp - To install the whole package (recommended!), download the file or_zuk_matlab_utils.tgz at the bottom of this page to <br>';
S{ctr+17} = '&nbsp a directory of your choice. Unzip it to get a tree containing all functions. Open the file ''SetPathScript.m''<br>';
S{ctr+18} = '&nbsp which is located in the ''file'' sub-directory and modify it by setting the path of the matlab libs (which is in <br>';
S{ctr+19} = '&nbsp the variable matlab_libs_root_dir) to your preffered directory. <br>';
S{ctr+20} = '&nbsp Then type ''SetPathScript'' in Matlab and this will add the package to your path. <br>';
S{ctr+21} = '&nbsp - To install a specific package, go to its link below, and download the .gz file at the bottom of the its page.<br>';
S{ctr+22} = '&nbsp Unzip to a directory of your choice in Matlab''s path. Make sure you download all other packages on which <br>';
S{ctr+23} = '&nbsp this package depends (some functions in certain packages depend on functions from other packages). <br>';
S{ctr+24} = '&nbsp - To install a specific function, just download it to a directory included by your Matlab path. Make sure to<br>';
S{ctr+25} = '&nbsp download all other functions on which this function depends.<br>';

ctr = length(S)+1;
lib_dirs = get_subdirectories(matlab_libs_dir, 1) % here get only the immediate libraries at the root
[file_names file_permissions] = GetFilePermissions(matlab_libs_dir, 1); % see for which dirs do we have permission
lib_dirs = intersect(lib_dirs, file_names(find(file_permissions(:,7))) );

root_lib_dirs = lib_dirs; lib_ctr=2;
for i=1:length(root_lib_dirs)
    lib_sub_dirs = get_subdirectories( fullfile(matlab_libs_dir, remove_dir_from_file_name(root_lib_dirs{i})), ...
	[], DFS, 7 ); % read only sub-dirs with read permission for everybody

    lib_dirs = insert_cell(lib_dirs, lib_sub_dirs, lib_ctr);
    lib_ctr = lib_ctr + 1 + length(lib_sub_dirs);
end


% good_lib_dirs = intersect(lib_dirs, file_names(find(file_permissions(:,7))) );
% good_lib_vec = zeros(length(lib_dirs), 1);
% for i=1:length(good_lib_dirs)
%     good_lib_vec(strmatch(good_lib_dirs{i}, lib_dirs)) = 1;
% end
% lib_dirs = lib_dirs(good_lib_vec);


num_dirs = length(lib_dirs)
lib_dirs = [lib_dirs 'or_zuk_matlab_utils.tgz']; % add another line for a zip file containing all package

for j=1:num_dirs+1
    matlab_html_lib_file = fullfile( strdiff(lib_dirs{j}, matlab_libs_dir), ...
        ['matlab_' remove_dir_from_file_name( lib_dirs{j} ) '_utils.html']);
    dir_depth = sum(strdiff(lib_dirs{j}, matlab_libs_dir) == delim);
    lib_dir_str = [repmat('&nbsp; ', 1, dir_depth*2) strdiff(lib_dirs{j}, matlab_libs_dir)];
    
    S{ctr} = '<tr>';
    if(j < num_dirs+1)
        S{ctr+1} = ['<td><a href="' matlab_html_lib_file '" name="' ...
            lib_dir_str ...
            '" id="' strdiff(lib_dirs{j}, matlab_libs_dir) '">' lib_dir_str '</a></td>'];
    else
        S{ctr+1} = ['<td><a href="' 'or_zuk_matlab_utils.tgz' '" name="' ...
            'or_zuk_matlab_utils' ...
            '" id="' 'or_zuk_matlab_utils' '">' 'or_zuk_matlab_utils.tgz' '</a></td>'];
    end
    S{ctr+2} = '';
    S{ctr+3} = '<td>&nbsp;-&nbsp;</td>';
    S{ctr+4} = '';
    if(j < num_dirs+1)
        f = strmatch(remove_dir_from_file_name(lib_dirs{j}), comment_convertion_struct(:,1), 'exact');
        if(isempty(f))
            f = [remove_dir_from_file_name(lib_dirs{j}) ' utilities '];
        else
            f = comment_convertion_struct{f(1),2};
        end
        S{ctr+5} = ['<td>' f  '</td>'];
        %                S{ctr+5} = ['<td>' remove_dir_from_file_name(lib_dirs{j}) ' utilities '  '</td>'];
        
    else
        S{ctr+5} = ['<td><b><i>Download the entire library. (' num2str(num_dirs) ' libs) </i></b></td>'];
    end
    S{ctr+6}  = '</tr>';
    ctr=ctr+7;
end
S = [S' R_end']';
S = strrep_cell(S, 'Matlab Template', 'Or Zuk''s Matlab ');
i = strfind_cell(S, '<UL>'); j = strfind_cell(S, '</UL>'); S = [S(1:i(1)-1)' S(j(1)+1:end)'];% Remove first automatic bullet line


% Add counter
counter_str = {'<a href="http://www.digits.com" target="_blank">', ...
    '<img src="http://counter.digits.com/?counter={6b8439e1-a985-ffa4-b56f-591db7277b21}&template=simple"', ...
    'alt="Hit Counter by Digits" border="0"  />', ...
    '</a> <br>'};

S =  insert_cell(S, counter_str, my_strmatch('Contact', S));
savecellfile(vec2column(S), matlab_master_html_libs_file); % save master html file
SS = S; % save for later 



cur_dir = pwd; % zip all files of all libraries, including docomuntation into one file
eval(['cd ' matlab_libs_dir]);

tar_dirs = remove_dir_from_file_name(lib_dirs(1:end-1)) % don't unzip data big data dirs 
data_dirs = tar_dirs(strfind_cell(tar_dirs, 'data'))
tar_dirs = setdiff(tar_dirs, data_dirs)

% return; % TEMP FOR DEBUG
zip_str = ['tar cvzf '  'or_zuk_matlab_utils.tgz --exclude "*/data/*" ' cell2vec(tar_dirs, ' ') ]; % Zip all directories %  fullfile(sub_dirs{i},'*.m')]
system(zip_str); % need to check that it works
eval(['cd ' cur_dir]);

sub_dirs = get_subdirectories(matlab_libs_dir, [], [], 7);
num_dirs = length(sub_dirs);  % Now prepare a file for each sub-directory
num_mat_files_vec = zeros(num_dirs,1);
for i=1:num_dirs % loop on libraries 
    lib_str = remove_dir_from_file_name( sub_dirs{i} ); lib_str(1) = upper(lib_str(1)); % Force upper 
%     if(lib_str(1) ~= 'G') % deal only with genetic architecture
%         continue;
%     end
    S = R_start;
    S = strrep_cell(S, 'Template', lower(lib_str));
    S = strrep_cell(S, 'Utilities', 'utilities'); 
    
    % Get some titles specific to library
    lib_ind = strfind_cell(lower(comment_convertion_struct(:,1)), lower(lib_str))
%    comment_struct_is = comment_convertion_struct
    if( ((~isempty(lib_ind)) & (length(lib_ind) == 1)) & ...
            (strmatch(lower(comment_convertion_struct{lib_ind,1}), lower(lib_str), 'exact') & ...                
        (~isempty(comment_convertion_struct{lib_ind,3}))) )
%        S{21}  = [S{21}(1:end-5) ' <br> ' comment_convertion_struct{lib_ind,3} ' </LI>']; % taking S{21} gives header 
        S{34}  = [ comment_convertion_struct{lib_ind,3} ' </LI>']; % taking S{21} gives header 

    end    
    
    mat_files = GetFileNames(fullfile(sub_dirs{i}, '*.m'));
    other_file_extensions = {'*.c', '*.cpp', '*.mex*', '*.pl', '*.txt', '*.py'}; % enable other source files 
    other_files = []; 
    for j=1:length(other_file_extensions)
        other_files = [other_files GetFileNames(fullfile(sub_dirs{i}, ...
            other_file_extensions{j}))]; 
    end    
%     other_files = union( union(GetFileNames(fullfile(sub_dirs{i}, '*.cpp')), ...
%         GetFileNames(fullfile(sub_dirs{i}, '*.c'))), ...
%         union( GetFileNames(fullfile(sub_dirs{i}, '*.h')), ...
%         GetFileNames(fullfile(sub_dirs{i}, '*.mex*')) ) ); %  a list of other files needed for library (not only .mat)
%     other_files = union(other_files, GetFileNames(fullfile(sub_dirs{i}, '*.pl')) ); 
%     other_files = union(other_files, GetFileNames(fullfile(sub_dirs{i}, '*.txt')) ); 
    num_mat_files = length(mat_files); 
    mat_files = [mat_files other_files]; % add other non-.mat files 
    num_files = length(mat_files);
    num_mat_files_vec(i) = num_files;
    mat_files = [mat_files [remove_dir_from_file_name( sub_dirs{i} ) '.tgz']]; % add another line for a zip file containing all package
    
    % Add links to upwards backward directories
    link_inds = strfind_cell(R, 'navbar');
    dir_words = strsplit(sub_dirs{i}, delim);
    start_ind = length(strsplit(matlab_libs_dir, delim))-1;
    up_dirs = cell(length(dir_words) - start_ind, 1);
    up_dirs{1} = remove_dir_from_file_name(matlab_libs_dir(1:end-1)); % remove the delimiter
    for j=2:length(up_dirs)
        up_dirs{j} = [up_dirs{j-1} '/' dir_words{j+start_ind}];
    end
    tmp = cell(length(up_dirs),1);
    for j=1:length(up_dirs)
        if(j == 1)
            tmp_file =  'matlab_utils.html';
        else
            tmp_file =  ['matlab_' remove_dir_from_file_name( up_dirs{j} ) '_utils.html'];
        end
        
        %         <a href="http://www.broadinstitute.org/~orzuk/matlab/index.html">Home</a> :
        tmp{j} = ['<a href="http://www.broadinstitute.org/~orzuk/matlab/' ...
            fullfile(up_dirs{j}, tmp_file)  '">' ...
            remove_dir_from_file_name(up_dirs{j}) '</a>  :'];
    end
    sub_sub_dirs = get_subdirectories(sub_dirs{i}, [], DFS, 7); % Add all subdirectories (DFS order)
    [sub_file_names sub_file_permissions] = GetFilePermissions(sub_dirs{i}, 1); % see for which dirs do we have permission
    sub_sub_dirs = intersect(sub_sub_dirs, sub_file_names(find(sub_file_permissions(:,7))) );

    % Make sure we've got permissions
    
    num_sub_dirs = length(sub_dirs);
    mat_files = [mat_files sub_sub_dirs]; % add possible packages in sub-directories 
    %     for j=1:length(sub_sub_dirs)
    %
    %     end
    
    
    S = insert_cell(S, tmp, link_inds(1)+1);
    ctr = length(S)+1;
    for j=1:length(mat_files) % num_files+1 % loop over matlab files and add a line per file
        %        j_is = j
        S{ctr} = '<tr>';
        if(j < num_files+1)
            if(j <= num_mat_files) % .mat file
                S{ctr+1} = ['<td><a href="' fullfile(sub_dirs{i}, mat_files{j}) ...
                    '" name="' mat_files{j}(1:end-2) ...
                    '" id="' mat_files{j}(1:end-2) '">' mat_files{j} '</a></td>'];
            else % other file
                S{ctr+1} = ['<td><a href="' fullfile(sub_dirs{i}, mat_files{j}) ...
                    '" name="' mat_files{j}(1:end-2) ...
                    '" id="' mat_files{j}(1:end-2) '"><font color="003300">' mat_files{j} '</a></td>'];
            end
        else
            if(j == num_files+1) % ??? 
                S{ctr+1} = ['<td><a href="' fullfile(sub_dirs{i}, mat_files{j}) ...
                    '" name="' mat_files{j}(1:end-4) ...
                    '" id="' mat_files{j}(1:end-4) '"><STRONG><i>' mat_files{j} '</i></STRONG></a></td>'];
            else % sub-directories
                dir_depth = sum(strdiff(mat_files{j}, sub_dirs{i}) == delim);
                lib_subdir_str = [repmat('&nbsp; ', 1, dir_depth*2) strdiff(mat_files{j}, [sub_dirs{i} delim])];
                
                matlab_html_sub_lib_file = ['matlab_' remove_dir_from_file_name( mat_files{j} ) '_utils.html'];
                S{ctr+1} = ['<td><a href="' fullfile(strdiff(mat_files{j}, [sub_dirs{i} delim]), matlab_html_sub_lib_file) ...
                    '" name="' lib_subdir_str '" id="' lib_subdir_str '"><STRONG><i>' lib_subdir_str '</i></STRONG></a></td>'];
                
                %                 S{ctr+1} = ['<td><a href="' fullfile(remove_dir_from_file_name(mat_files{j}), matlab_html_sub_lib_file) '" name="' ...
                %                     remove_dir_from_file_name(mat_files{j}(1:end-4)) ...
                %                     '" id="' fullfile(remove_dir_from_file_name(mat_files{j}(1:end-4)), matlab_html_sub_lib_file) '"><STRONG><i>' ...
                %                     remove_dir_from_file_name(mat_files{j}) '</i></STRONG></a></td>'];
                
            end
        end
        S{ctr+1} = strrep(S{ctr+1}, '~orzuk/public_html', 'http://www.broadinstitute.org/~orzuk'); 
%         replaced_S = S{ctr+1}
                
        S{ctr+2} = '';
        S{ctr+3} = '<td>&nbsp;-&nbsp;</td>';
        S{ctr+4} = '';
        if(j < num_files+1)
            fid = fopen( fullfile(sub_dirs{i}, mat_files{j}), 'r');
            tmp = fgetl(fid);
            while(isempty(tmp) || (tmp(1) ~= '%'))
                tmp = fgetl(fid);
                if(isnumeric(tmp)) % this means we reached the end of file - so read again the first line
                    fclose(fid);
                    fid = fopen( fullfile(sub_dirs{i}, mat_files{j}), 'r');
                    tmp = fgetl(fid);
                    break;
                end
            end
            fclose(fid);
        else
            if(j == num_files+1)
                tmp = ['% <b><i>Download the whole package. (' num2str(num_mat_files) ' .m files) </i></b>'];
            else % sub-directories
                tmp = ['% <b><i>Go to the ' remove_dir_from_file_name(mat_files{j}) ' package. </i></b>'];
            end
        end
        S{ctr+5} = ['<td>' tmp(3:end)  '</td>'];
        S{ctr+6}  = '</tr>';
        ctr=ctr+7;
    end    
    S = [S' R_end']';
    
%     S_IS = S 
%     save_S = matlab_html_lib_file
    
    matlab_html_lib_file = fullfile(sub_dirs{i}, ['matlab_' remove_dir_from_file_name( sub_dirs{i} ) '_utils.html']);
    savecellfile(S, matlab_html_lib_file); % save html file
    
    cur_dir = pwd; % zip the files into one file
    cd_to_dir = sub_dirs{i}
    eval(['cd ''' sub_dirs{i} '''']);
    zip_str = ['tar cvzf ' remove_dir_from_file_name( sub_dirs{i} ) '.tgz' ' *.m *.cpp *.h *.mex*'] %  fullfile(sub_dirs{i},'*.m')]
    system(zip_str); % need to check that it works
    eval(['cd ' cur_dir]);
end % loop on sub-directories


% Update number of files (not working for now)
ctr = strfind_cell(SS, 'Download the entire library');
SS{ctr} = ['<td><b><i>Download the entire library. (excluding big data files. ' num2str(num_dirs) ...
    ' libs, ' num2str(sum(num_mat_files_vec)) ' files) </i></b></td>'];
%savecellfile(vec2column(SS), matlab_master_html_libs_file); % save master html file


