% function AssignGeneticArchitectureConstants

catalog_type = 'MIT-TURK'; % 'NHGRI'; % 'MIT-TURK';
DISCOVERY = 1; REPLICATION = 2; COMBINED = 3;

special_traits = {'Height', 'Menarche (age at onset)', 'Breast Cancer'}

% , 'Body mass index', 'Crohn''s disease'}; %  ... %   }; % debug only height
%%    'Breast Cancer', 'Type 1 diabetes', 'Type 2 diabetes', ...    % 'Crohn''s disease', 'Height', ...
%%    'Lipids', 'Lipid LDL', 'Lipid HDL', 'Lipid TriGly'}; % , ...
%    'Menarche (age at onset)', 'HbF'};
MIN_POWER = 0.05; MIN_POWER_LIBERAL = 0.0001; % truncate power correction. Minimal allowed power
MAX_CORRECTION = 20; % 1 over MIN_POWER

correction_mode = 'floor'; % 'exact'; % 'floor'; % 'round', 'ceil'

% Set different power corrections:
correction_type = {'theoretical', 'pop-gen'}; % 'theoretical' % 'empirical'; 'pop-gen', 'Park'

% Indices: 11 - Park with variable beta, 3 - pop. gen. [0.01,0.5] (common), fitted s
correction_inds = [11 3]; % indices in the correction strings (should be computed automatically rather than hard-coded)

lambda_type = {{'MZ-twins'}, {'DZ-twins', 'sibs', 'parent-offspring'}, ...
    {'grandparent-grandchild', 'half-sibs', 'uncle/ant'}, ...
    {'1st-cousins', 'graet-grandparent-great-grand-child'}, ...
    {'1st-cousing-once-removed', 'great-great-grandfother'}, {'2nd-cousings'}}; % different family relations 
