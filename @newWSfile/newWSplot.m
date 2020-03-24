function newWSplot(obj, varargin)

    % optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
    optargs = {-inf inf -inf inf};
    optargs(1:numvarargs) = varargin;
    [xmin, xmax, ymin, ymax] = optargs{:};
    
    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj)
    
    % forcing total number of active channels to THREE
    totalActiveChannels = 3;
    
    % plotting one figure with all channels and all sweeps
    figure('name', strcat(obj.file,' (all) - all Channels')); % naming figure file
    
    % separate code for channel 1 to apply axis range and figure title
    for sweepNumber = allSweeps  
        hold on;
        [x,y] = obj.xy(sweepNumber, 1);
        subplot(totalActiveChannels,1,1)
        plot(x,y);
        axis([xmin xmax ymin ymax])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([obj.file ' (all)'],'Interpreter','none');
    end
    hold off;    
    
    % code for all other channels
    for channel = 2:totalActiveChannels
        for sweepNumber = allSweeps
        hold on;
        [x,y] = obj.xy(sweepNumber, channel);     
        subplot(totalActiveChannels,1,channel)
        plot(x,y);
        yminhere = min(y)-5;
        ymaxhere = max(y)+5;
        axis([xmin xmax yminhere ymaxhere])
        ylabel(strcat(obj.header.AIChannelNames(channel), ' (', obj.header.AIChannelUnits(channel), ')'));
        end
        xlabel('Time (s)');
        hold off;
    end
    movegui('northwest');
    
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),') - all Channels'));
        [x,y] = obj.xy(sweepNumber, 1);
        subplot(totalActiveChannels,1,1)
        plot(x,y);
        axis([xmin xmax ymin ymax])
        ylabel(strcat(obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        
        for channel = 2:totalActiveChannels
            [x,y] = obj.xy(sweepNumber, channel);     
            subplot(totalActiveChannels,1,channel)
            plot(x,y);
            yminhere = min(y)-5;
            ymaxhere = max(y)+5;
            axis([xmin xmax yminhere ymaxhere])
            ylabel(strcat(obj.header.AIChannelNames(channel), ' (', obj.header.AIChannelUnits(channel), ')'));
        end
        xlabel('Time (s)');
        movegui('north');
    end 
end