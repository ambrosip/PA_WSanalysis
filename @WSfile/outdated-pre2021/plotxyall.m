function plotxyall(obj, varargin)

    % EX: plotxyall(m293.s0061);
    % EX: plotxyall(m293.s0061,1);
    % EX: plotxyall(m293.s0061,1,-inf,inf,-1000,500);    
    
    % optional arguments: axis range
    numvarargs = length(varargin);
    optargs = {1 -inf inf -inf inf};
    optargs(1:numvarargs) = varargin;
    [channel, xmin, xmax, ymin, ymax] = optargs{:};
    
    % finding sweep numbers from file name
    if length(obj.file) == 28
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;
    
    % plotting one figure with all sweeps
    figure('name', strcat(obj.file,' (all)')); % naming figure file
    hold on;
    for sweepNumber = allSweeps
        [x,y] = obj.xy(sweepNumber, channel);
        plot(x,y);
    end
	hold off;
    axis([xmin xmax ymin ymax])
    xlabel('Time (s)');
    ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
    title([obj.file ' (all)'],'Interpreter','none');
    movegui('northwest');
    
    % plotting individual figures for all sweeps
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),')'));
        [x,y] = obj.xy(sweepNumber, channel);
        plot(x,y)
        axis([xmin xmax ymin ymax])
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        movegui('north');
    end

end


