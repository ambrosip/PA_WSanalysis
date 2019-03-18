function WSfileFromFolder(dirName)

    files=dir(fullfile(dirName, '*.h5'));
    fileNames={files.name}';
    
    
    
    for i=1:length(fileNames)
        strcat(fileNames{i}(1:5),fileNames{i}(end-11:end-8)) = WSfile(fileNames{i});       
    end

end
