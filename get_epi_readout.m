function get_epi_readout
% calculates the EPI readout time - relevant for fieldmap based distortion
%  correction

fprintf('Please select one functional image of which you would like to know the TOTAL EPI READOUT TIME\n\n')

file = spm_select;
dinfo = dicominfo(file);
bandperpixphenc = dinfo.Private_0019_1028; % bandwith per pixel phase encode
tert = 1000*(round((1/bandperpixphenc),5)); % total EPI readout time

fprintf('The TOTAL EPI READOUT TIME is: %s \n\n', num2str(tert));
fprintf('You can get the short and long TE information from fieldmap sequence info sheet!\n') ;

end