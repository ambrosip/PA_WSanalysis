function niceplot(obj, varargin)

% optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
    optargs = {15 3 10 23 -85 50};
    optargs(1:numvarargs) = varargin;
    [lightOnsetTime, lightDuration, xmin, xmax, ymin, ymax] = optargs{:};
    
    % finding sweep numbers from file name
    if length(obj.file) == 28
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));        
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;     
    
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),') - niceplot'));
        [x,y] = obj.xy(sweepNumber, 1);
        plot(x,y,'k','LineWidth',1.5);
        axis([xmin xmax ymin ymax+10])
%         xlabel('Time (s)');
%         ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        set(gca,'Visible','off')
        rectangle('Position', [lightOnsetTime 55 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    end 
end