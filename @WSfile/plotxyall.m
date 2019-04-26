function plotxyall(obj, varargin)
    
    % optional arguments: axis range
    numvarargs = length(varargin);
    optargs = {-inf inf -inf inf};
    optargs(1:numvarargs) = varargin;
    [xmin, xmax, ymin, ymax] = optargs{:};
    
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
    figure;
    hold on;
    for sweepNumber = allSweeps
        [x,y] = obj.xy(sweepNumber);
        plot(x,y)
    end
	hold off;
    axis([xmin xmax ymin ymax])
    xlabel('Time (s)');
    ylabel(obj.header.Acquisition.AnalogChannelUnits(1));
    title([obj.file ' (all)'],'Interpreter','none');
    
    % plotting individual figures for all sweeps
    for sweepNumber = allSweeps
        figure;
        [x,y] = obj.xy(sweepNumber);
        plot(x,y)
        axis([xmin xmax ymin ymax])
        xlabel('Time (s)');
        ylabel(obj.header.Acquisition.AnalogChannelUnits(1));
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
    end
    

end


