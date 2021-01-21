function mouseNumber = newWSfileFromDir(dirName, varargin)

    % EX: m011 = WSfileFromDir('D:\Temp\Ephys\20190426')
    
    addpath(dirName)

    % set defaults for optional inputs 
    optargs = {17 20};  % For files named like m052_2020-08-14_0001-0018
%     optargs = {16 19};  % For files named like m52_2020-08-14_0001-0018
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % Place optional args in memorable variable names
    [startOfSweepNumberInFileName,endOfSweppNumberInFileName] = optargs{:};      

    files=dir(fullfile(dirName, '*.h5'));
    fileNames={files.name}';
    
    mouseNumber = struct;
    
    for i=1:length(fileNames)
        field = strcat('s',fileNames{i}(startOfSweepNumberInFileName:endOfSweppNumberInFileName));      
        mouseNumber.(field) = newWSfile(fileNames{i});
    end
end