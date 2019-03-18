function mouseNumber = WSfileFromDir(dirName)

    files=dir(fullfile(dirName, '*.h5'));
    fileNames={files.name}';
    
    mouseNumber = struct;
    
    for i=1:length(fileNames)
        field = strcat('s',fileNames{i}(17:20));        
        mouseNumber.(field) = WSfile(fileNames{i});
    end
end

