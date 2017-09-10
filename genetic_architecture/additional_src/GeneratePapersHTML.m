% Create a web-page for genetic-architecture papers 
function GeneratePapersHTML()

AssignGeneralConstants;
papers_dir = '~eliana/public_html/Papers';
matlab_html_template_file = fullfile(matlab_libs_root_dir, 'matlab_utilities_template.html');

F = GetFileNames(fullfile(papers_dir, '*.pdf'));

R = loadcellfile(matlab_html_template_file); % load template file
start_ind = strfind_cell(lower(R), '<tr>');
R_start = R(1:start_ind(1)); R_end = R(start_ind(end):end);
S = R_start; ctr = length(S)+1;
X = loadcellfile(fullfile(papers_dir, 'genetic_architecture_papers.txt'));

[common_names I J] = intersect(F, X(:,1)); 

for j=1:length(F) % num_files+1 % loop over matlab files and add a line per file
    %        j_is = j
    S{ctr} = '<tr>';
    S{ctr+1} = ['<td><a href="' F{j} '" name="' F{j}(1:end-2) ...
        '" id="' F{j}(1:end-2) '">' F{j} '</a></td>'];
    [A B] = ismember(j, I);
    if(A) % found intersection
        S{ctr+1} = [S{ctr+1} ' <td> ' X{J(B),2}  ' </td> '];
    end
    ctr=ctr+2;
end


papers_file = fullfile(papers_dir, 'papers_index.html');
%papers_file = 'papers_index.html';
S = [S' R_end']';

S = strrep_cell(S, 'Matlab Template Utilities', 'Papers for Genetic Architecture Project ');
S = strrep_cell(S, 'see documentation and examples within the functions', ...
    'short description follows every file ');

savecellfile(S, papers_file); % save html file

