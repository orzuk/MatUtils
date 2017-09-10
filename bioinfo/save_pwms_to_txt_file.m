% Saves a set of pwms in a convenient .txt format 
%
% Input: 
% pwms - a cell array of pwms, in a 4*n format
% output_file - a name of a text file to output
% fasta format - a different format to save the pwms (only name, and a '>' mark)
%
function Dummy = save_pwms_to_txt_file(pwms, output_file, fasta_format, varargin)

if(ischar(pwms)) % enable input from file 
    load(pwms); 
end
if(~exist('fasta_format', 'var') || isempty(fasta_format)) % default is no fasta format 
    fasta_format = 0; 
end
n = size(pwms, 1); % number of pwms
% q = {}; 
if(~fasta_format)
    q = cell(1, n*8);
else
    q = cell(1, n*5);
end
ctr = 1;
for i=1:n
    if(~fasta_format)
        q{ctr} = pwms{i,1};
        q{ctr+1} = pwms{i,3};
        q{ctr+2}  = pwms{i,4};
        ctr = ctr+2;
    else
        q{ctr} = ['>' pwms{i,1}];
    end
    for j=1:4
        q{ctr+j} = num2str_delim(pwms{i,2}(j,:), '\t'); % seperate by tabs
    end
    if(~fasta_format)
        q{ctr+5} = ' ';
        ctr = ctr+6;   
    else
        ctr = ctr+5;
    end
end
q = q';
savecellfile(q, output_file); 
Dummy = [];


