function saveAllFigs(dirName)

    % EX:
    % saveAllFigs('R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data
    % summaries\2019-04-29')

    % saves all figures in destination folder set by dirName (remember to
    % add single quotes around the directory) and uses the figure name as
    % file name.

    FolderName = dirName;   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName = FigList(iFig).Name;
      set(0, 'CurrentFigure', FigHandle);
      saveas(FigHandle,fullfile(FolderName, [FigName '.tiff']));
    end

end