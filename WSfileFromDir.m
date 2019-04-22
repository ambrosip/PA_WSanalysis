function mouseNumber = WSfileFromDir(dirName, varargin)

    % set defaults for optional inputs 
    optargs = {17 20};
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
        mouseNumber.(field) = WSfile(fileNames{i});
    end
end

