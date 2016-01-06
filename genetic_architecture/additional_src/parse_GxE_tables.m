% Convert GxE tables to latex
function parse_GxE_tables()

AssignGeneralConstants;
tables_dir = '../../common_disease_model/docs/pnas/genetic_interactions/tables/';
trait_vec = {'raffinose', 'grapejuice', 'galactose', 'exudate'};

for i=1:length(trait_vec)
    tables_file = fullfile(tables_dir, ['GxE_glucose_' trait_vec{i} '_tab.txt']);
    R = loadcellfile(tables_file); % load data
    for j=1:size(R,1) % conver to percent
        for k=1:size(R,2)-1
            R{j,k} = 100*R{j,k};
        end
    end
        
    R = [{'$f (\%)$', '$d (\%)$', '$\mu (\%)$', ...
        '$V_P (\%)$', '$V_G (\%)$', '$V_E (\%)$', '$V_{GxE} (\%)$', '$\pi_{phantom} (\%)$'}' R']'; 
    for kk=2:size(R,1)
        for j=3:size(R,2)
            R{kk,j} = sprintf('%.1f', R{kk,j});
        end
    end
%    R = num2str_cell(R, 2);     
    
    R_latex = latex(R, 2, precision);
    R_latex = mat2cell(R_latex, ones(size(R_latex,1),1));
    R_latex{1} = strrep(R_latex{1},  '|c', '|r'); % align to right  
    tab_header = {'\begin{table}[h!] % table summarizing all diseases', ...
        '\begin{center}', '{\scriptsize', ...
        ['\caption[Phantom heritability due to GxE: Glucose vs. ' trait_vec{i} ']' ...
        '{Phantom heritability due to GenexEnvironment interaction for Glucose and ' trait_vec{i}  ...
        ' as a function of allele frequency f and environment frequency d. }'], ...
        R_latex{1}};
    tab_footer = {'}',  ...
        '\label{table:fig1_data_quantitative_traits}', ...
        '\end{center}', '\end{table}', '', '\newpage', ''};
    % S_latex = [tab_header  R_latex(2:end)' tab_footer]';
    R_latex = split_latex_table(R_latex(2:end), tab_header, ['\end{tabular}' tab_footer], 60);
    savecellfile(R_latex, strdiff(tables_file, '_tab'), [], 1);
end

