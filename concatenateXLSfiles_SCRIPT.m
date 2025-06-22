dirName='/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/Data Analysis/2025-02-05 analysis window 15 ms';

% function concatenateXLSfiles(dirName)

addpath(dirName)
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% %% concatenating all XLS files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(files(1).name(1:15)," - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);




%% concatenating all files with a particular suffix in DIR into a single XLS file 

% find all the xls files in dir 
% files=dir(fullfile(dirName, 'DLS*', '*sCRACM_40x_5x5_ROI_by_ROI.xls'));
files=dir(fullfile(dirName, '*sCRACM_40x_5x5_ROI_by_ROI.xls'));

% get the names of the files
filesNames={files.name}';

% create array to store all xls data
concatenatedXLS = [];

% append all xls files
for i=1:numel(filesNames)
    concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
end

% save xls file with data 
filename = strcat(analysisDate, " - LED - sCRACM_40x_5x5_ROI_by_ROI - concatenated.xls");
% filename = strcat(analysisDate, " - sCRACM_40x_5x5_ROI_by_ROI - concatenated");
fulldirectory = fullfile(dirName,filename);
% fulldirectory = strcat(dirName,'\',filename,'.xls');
writetable(concatenatedXLS,fulldirectory);
% writematrix(concatenatedXLS,fulldirectory);



% %% concatenating all 'psc_vs_light - cell.xls' files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*PA_FP_stim_DualSite_rewrite_MEAN.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(analysisDate, "PA_FP_stim_DualSite_rewrite_MEAN - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);

 
% %% concatenating all 'psc_vs_light - cell.xls' files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*first_psc_kinetics - cell.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(files(1).name(1:15)," - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);


% %% concatenating all 'first_psc_kinetics - cell.xls' files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*first_psc_kinetics - cell.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(files(1).name(1:15)," first_psc_kinetics - cell - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);


% %% concatenating all confocal analysis.xls files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*confocal analysis.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(analysisDate," - confocal analysis - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);


% %% concatenating all cell_avgs XLS files in DIR into a single XLS file 
% 
% % find all the xls files in dir 
% files=dir(fullfile(dirName, '*cell_avgs.xls'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % create array to store all xls data
% concatenatedXLS = [];
% 
% % append all xls files
% for i=1:numel(filesNames)
%     concatenatedXLS = [concatenatedXLS; readtable(files(i).name)];
% end
% 
% % save xls file with data 
% filename = strcat(analysisDate," - cell_avgs - concatenated");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedXLS,fulldirectory);

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



