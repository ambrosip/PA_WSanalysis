function saveAllFigs(varargin)


    % optional arguments
    % set defaults for optional inputs 
%     optargs = {'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
%     optargs = {'D:\Temp\From MATLAB'};
%     optargs = {'E:\From MATLAB'};
    optargs = {'D:\CORONAVIRUS DATA\2021-09-21 ipsc kinetics'};
 
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [dirName] = optargs{:};

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
      
      % forces matlab to save fig as a vector
      FigHandle.Renderer = 'painters';
      
      % actually saves a vector file
      saveas(FigHandle,fullfile(FolderName, [FigName '.svg']));
    end
    
    disp('I did it')

end