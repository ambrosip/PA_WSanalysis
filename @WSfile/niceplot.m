function niceplot(obj, varargin)

% optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
%     optargs = {15 3 10 23 -85 50};
    optargs = {-200 200 1 20 15 3 10 23};
    optargs(1:numvarargs) = varargin;
    [ymin, ymax, filterORnot, GaussianFilterWindow, lightOnsetTime, lightDuration, xmin, xmax] = optargs{:};
    
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
        yfiltered = smoothdata(y,'gaussian',GaussianFilterWindow);
%         plot(x,y,'k','LineWidth',1.5);
        plot(x,y,'k','LineWidth',1);
        axis([xmin xmax ymin ymax+10])
%         xlabel('Time (s)');
%         ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        set(gca,'Visible','off')
        
        % adding light stim
%         rectangle('Position', [lightOnsetTime 55 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%         rectangle('Position', [lightOnsetTime (ymax-((ymax-ymin)/20)) lightDuration ymax], 'FaceColor', [0.4660, 0.6740, 0.1880], 'LineStyle', 'none')
%         line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0.4660, 0.6740, 0.1880],'LineWidth',10)
        line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)

        % adding scale bar
        line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    
        if filterORnot
            figure('name', strcat(obj.file,' (',num2str(sweepNumber),') - niceplot filtered'));
            plot(x,yfiltered,'k','LineWidth',1);
            axis([xmin xmax ymin ymax+10])
            title([obj.file ' (' num2str(sweepNumber) ") filtered"],'Interpreter','none');
            set(gca,'Visible','off')
%             line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0.4660, 0.6740, 0.1880],'LineWidth',10)
            line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
            line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
            line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
            text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
            text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))        
        end    
    end             
end
