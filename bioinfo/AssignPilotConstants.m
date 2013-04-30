INTERNAL = 0; TOPHAT = 1;
solexa_cerebellum_dir = '/seq/orzuk/24mammals/data/rna_seq/human_mouse_cerebellum_pilot/TopHat';
solexa_lung_dir = '/seq/rinnscratch/mguttman/RNASeq/LungFB/TopHat'; % Only Human!!!
rna_matrix_output_file = ...
    fullfile(solexa_cerebellum_dir, 'solexa_matrix_rna.mat');
user = 'NIDA'; % NIDA
agilent_prenormalization_flag = 'OLD'; % OLD - previous normalization. NEW - new normalization


switch user
    case 'OR'
        switch agilent_prenormalization_flag
            case 'OLD'
                agilent_dir = '../../mammals_expression/Agilent/pilot/';
            case 'NEW'
                agilent_dir = '../../mammals_expression/Agilent/pilot/new_protocol';
        end
    case 'NIDA'
        switch agilent_prenormalization_flag
            case 'OLD'
                agilent_dir = '~/agilent';
            case 'NEW'
                agilent_dir = '~/agilent/new_protocol';
                
        end
end

agilent_pilot_array_date = '';
barcode_mapping_file = 'barcode_mapping.txt';
