function concatenateXLSfiles(dirName)

%% concatenating all XLS files in DIR into a single XLS file 

% find all the xls files in dir 
files=dir(fullfile(dirName, '*.xls'));

% get the names of the files
filesNames={files.name}';

% create array to store all xls data
concatenatedXLS = [];

% append all xls files
for i=1:numel(filesNames)
    concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
end

% save xls file with data 
filename = strcat(files(1).name(1:15)," - concatenated");
fulldirectory = strcat(dirName,'\',filename,'.xls');
writetable(concatenatedXLS,fulldirectory);

% %% concatenating PPR
% 
% % find all the files in dir that contain the string "PPR" and have a
% % .xls extension
% pprfiles=dir(fullfile(dirName, '*PPR*.xls'));
% 
% % get the names of the files
% pprfilesNames={pprfiles.name}';
% 
% % create array to store all xls data
% concatenatedXLSppr = [];
% 
% % append all xls files
% for i=1:numel(pprfilesNames)
%     concatenatedXLSppr = [concatenatedXLSppr; readtable(pprfiles(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(pprfiles(1).name(1:15)," - concatenated PPR");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLSppr,fulldirectory);


% %% concatenating AMPA NMDA ratio
% 
% % find all the files in dir that contain the string "PPR" and have a
% % .xls extension
% pprfiles=dir(fullfile(dirName, '*ratio*.xls'));
% 
% % get the names of the files
% pprfilesNames={pprfiles.name}';
% 
% % create array to store all xls data
% concatenatedXLSratio = [];
% 
% % append all xls files
% for i=1:numel(pprfilesNames)
%     concatenatedXLSratio = [concatenatedXLSratio; readtable(pprfiles(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(pprfiles(1).name(1:15)," - concatenated AMPA NMDA ratio");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLSratio,fulldirectory);


end
