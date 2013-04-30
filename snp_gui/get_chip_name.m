function chip_name=get_chip_name(chip_name)


if strcmp(chip_name,'Hind50k ( from Affy 100K)')
    chip_name='Hind';
elseif strcmp(chip_name,'Xba50k ( from Affy 100K)')
    chip_name='Xba';
    %     else
    %         errordlg('Invalid chip','Chip Error','modal');
    %         return;
end