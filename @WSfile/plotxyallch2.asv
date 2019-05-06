function plotxyallch2(obj, varargin)

    % plots all channels and all sweeps combined in single figure and in
    % separate figures per sweep.

    % EX: plotxyallch2(m293.s0061);
    % EX: plotxyallch2(m293.s0061,-inf,inf,-1000,500);    
    
    % optional arguments: axis range for channel 1 (monitor CC or VC)
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
    
    % finding total number of active channels from header
    totalActiveChannels = obj.header.Acquisition.ActiveChannelIndexFromChannelIndex(end);
    
    % plotting one figure with all channels and all sweeps
    figure('name', strcat(obj.file,' (all) - all Channels joined')); % naming figure file
    
    % separate code for channel 1 to apply axis range and figure title
    % select sweeps at odd indexes only
    for sweepNumber = allSweeps(1:2:end)  
        hold on;
        [x1,y1] = obj.xy(sweepNumber, 1);
        [x2,y2] = obj.xy(sweepNumber+1, 1);
        x2=x2+15;
        % join data from two sweeps
        x = [x1;x2];
        y = [y1;y2];
        subplot(totalActiveChannels,1,1)
        plot(x,y);
        axis([xmin xmax ymin ymax])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) joined'],'Interpreter','none');
    end
    hold off;    
    
    % code for all other channels
    for channel = 2:totalActiveChannels
        for sweepNumber = allSweeps(1:2:end)
        hold on;
        [x1,y1] = obj.xy(sweepNumber, channel);
        [x2,y2] = obj.xy(sweepNumber+1, channel);
        x2=x2+15;
        x = [x1;x2];
        y = [y1;y2];  
        subplot(totalActiveChannels,1,channel)
        plot(x,y);
        yminhere = min(y)-5;
        ymaxhere = max(y)+5;
        axis([xmin xmax yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(channel), ' (', obj.header.Acquisition.AnalogChannelUnits(channel), ')'));
        end
        hold off;
    end
    movegui('northwest');
    
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps(1:2:end)  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),') - all Channels joined'));
        [x1,y1] = obj.xy(sweepNumber, 1);
        [x2,y2] = obj.xy(sweepNumber+1, 1);
        x2=x2+15;
        % join data from two sweeps
        x = [x1;x2];
        y = [y1;y2];
        subplot(totalActiveChannels,1,1)
        plot(x,y);
        axis([xmin xmax ymin ymax])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (' num2str(sweepNumber) ') joined'],'Interpreter','none');
        
        for channel = 2:totalActiveChannels
            [x1,y1] = obj.xy(sweepNumber, channel);
            [x2,y2] = obj.xy(sweepNumber+1, channel);
            x2=x2+15;
            x = [x1;x2];
            y = [y1;y2];       
            subplot(totalActiveChannels,1,channel)
            plot(x,y);
            yminhere = min(y)-5;
            ymaxhere = max(y)+5;
            axis([xmin xmax yminhere ymaxhere])
            xlabel('Time (s)');
            ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(channel), ' (', obj.header.Acquisition.AnalogChannelUnits(channel), ')'));
        end
        movegui('north');
    end    
    
%     distFig()
end