% Save family tree information in pedigree format
% 
% Input:
% family_tree - tree representing the family
% names_vec - names of nodes (if empty use numberings)
% gender_vec - vector of genders: Male (1), Female (0), Unknown (-1)
% traits_vec - vector of traits for each individual
% pedigree_outfile - where to save the pedigree
% pedigree_format - 'Morgan' or 'Linkage'
%
function save_pedigree(family_tree, names_vec, gender_vec, traits_vec, ...
    pedigree_outfile, pedigree_format)

AssignGeneralConstants();

n = length(family_tree); % number of nodes
if(~exist('pedigree_format', 'var') || isempty(pedigree_format))
    pedigree_format = 'Morgan';
end
pedigree_sources = get_graph_sources(family_tree);

switch pedigree_format % choose format
    case 'Morgan' % Use the Morgan format: individual, father, mother, sex, observed, other
        R = cell(n, 4);
        %        S{1} =
        S = cell(1,4); S{1,1} = '********'; % file seperator
        for i=1:n
            R{i,1} = i;
            par = parents(family_tree, i);
            if(isempty(par)) % pedigree sources - they don't have parents
                R{i,2} = 0; R{i,3} = 0;
            else
                R{i,2} = par(1); R{i,3} = par(2); % set parents
            end
            R{i,4} = family_tree(i,i) + 1; % specify gender
        end % loop on individuals 
        R = [S' R']';
    case 'Linkage' % Use the Linkage format: family, individual, father, mother, sex, observed, other
        
    case 'pedigraph'
        R = cell(n, 5);
        for i=1:n
            R{i,1} = names_vec{i}; % i;
            par = parents(family_tree, i);
            if(isempty(par)) % pedigree sources - they don't have parents
                R{i,2} = 0; R{i,3} = 0;
            else
                R{i,2} = names_vec{par(1)}; R{i,3} = names_vec{par(2)}; % set parents
            end
            switch gender_vec(i)
                case MALE
                    R{i,4} = 'M';
                case FEMALE
                    R{i,4} = 'F';
                otherwise
                    R{i,4} = '?';
            end
            %            R{i,4} = family_tree(i,i) + 1; % specify gender
            R{i,5} = num2str(traits_vec(i));
        end % loop on individuals       
end % switch pedigree format 
savecellfile(R, pedigree_outfile);
