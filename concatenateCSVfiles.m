function concatenateCSVfiles(dirName)

addpath(dirName)

%% concatenating all CSV files in DIR into an XLS file

% find all the CSV files in dir 
files=dir(fullfile(dirName, '*.csv'));

% get the names of the files
filesNames={files.name}';

% create array to store all CSV data
concatenatedCSV = [];

% append all CSV files
for i=1:numel(filesNames)
    concatenatedCSV = [concatenatedCSV; readtable(files(i).name)];
end

% save XLS file with data 
filename = strcat(files(1).name(1:15)," - concatenated");
fulldirectory = strcat(dirName,'\',filename,'.xls');
writetable(concatenatedCSV,fulldirectory);


% %% concatenating light vs firing rate files
% 
% % find all the files in dir that contain the string "firing" and have a
% % .csv extension
% firingfiles=dir(fullfile(dirName, '*firing*.csv'));
% 
% % get the names of the files
% firingfilesNames={firingfiles.name}';
% 
% % create array to store all csv data
% concatenatedCSVfiring = [];
% 
% % append all csv files
% for i=1:numel(firingfilesNames)
%     concatenatedCSVfiring = [concatenatedCSVfiring; readtable(firingfiles(i).name)];
% end
% 
% % save csv file with data 
% filename = strcat(firingfiles(1).name(1:15)," - concatenated firing");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVfiring,fulldirectory);


% %% concatenating excitability files
% 
% excitabilityfiles=dir(fullfile(dirName, '*excitability*.csv'));
% excitabilityfilesNames={excitabilityfiles.name}';
% concatenatedCSVexcitability = [];
% 
% for i=1:numel(excitabilityfilesNames)
%     concatenatedCSVexcitability = [concatenatedCSVexcitability; readtable(excitabilityfiles(i).name)];
% end
% 
% filename = strcat(excitabilityfiles(1).name(1:15)," - concatenated excitability");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVexcitability,fulldirectory);


% %% concatenating CC sag files
% 
% CCsagfiles=dir(fullfile(dirName, '*CC sag*.csv'));
% CCsagfilesNames={CCsagfiles.name}';
% concatenatedCSVCCsag = [];
% 
% for i=1:numel(CCsagfilesNames)
%     concatenatedCSVCCsag = [concatenatedCSVCCsag; readtable(CCsagfiles(i).name)];
% end
% 
% filename = strcat(CCsagfiles(1).name(1:15)," - concatenated CC sag");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVCCsag,fulldirectory);
% 
% 
% %% concatenating VC sag files
% 
% VCsagfiles=dir(fullfile(dirName, '*VC sag*.csv'));
% VCsagfilesNames={VCsagfiles.name}';
% concatenatedCSVVCsag = [];
% 
% for i=1:numel(VCsagfilesNames)
%     concatenatedCSVVCsag = [concatenatedCSVVCsag; readtable(VCsagfiles(i).name)];
% end
% 
% filename = strcat(VCsagfiles(1).name(1:15)," - concatenated VC sag");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVVCsag,fulldirectory);
% 
% 
% %% concatenating Rs sag files
% 
% rsfiles=dir(fullfile(dirName, '*Rs*.csv'));
% rsfilesNames={rsfiles.name}';
% concatenatedCSVrs = [];
% 
% for i=1:numel(rsfilesNames)
%     concatenatedCSVrs = [concatenatedCSVrs; readtable(rsfiles(i).name)];
% end
% 
% filename = strcat(rsfiles(1).name(1:15)," - concatenated Rs");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVrs,fulldirectory);
% 
% 
% %% concatenating normmono files
% 
% normmonofiles=dir(fullfile(dirName, '*Baseline*.csv'));
% normmonofilesNames={normmonofiles.name}';
% concatenatedCSVnormmono = [];
% 
% for i=1:numel(normmonofilesNames)
%     concatenatedCSVnormmono = [concatenatedCSVnormmono; readtable(normmonofiles(i).name)];
% end
% 
% filename = strcat(normmonofiles(1).name(1:15)," - concatenated normmono");
% fulldirectory = strcat(dirName,'\',filename,'.xls');
% writetable(concatenatedCSVnormmono,fulldirectory);


%% encouraging statements

disp('I did it =3')   

end