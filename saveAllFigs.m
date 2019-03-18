function saveAllFigs(dirName)

    FolderName = dirName;   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      saveas(FigHandle,fullfile(FolderName, [FigName '.tiff']));
    end

end