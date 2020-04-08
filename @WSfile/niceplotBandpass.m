function niceplotBandpass(obj, varargin)

    % optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
%     optargs = {'black' 100 2000 -200 100 15 3 10 23};   % for ON data
    optargs = {'black' 100 2000 -20 10 15 3 10 23};   % for ON data
%     optargs = {'black' 100 1000 -200 150 15 3 10 23};   % for ON data
%     optargs = {'black' 2 4000 -500 100 15 3 10 23};   % for mIPSC data
    optargs(1:numvarargs) = varargin;
    [colorName, highpassThreshold, lowpassThreshold, ymin, ymax, lightOnsetTime, lightDuration, xmin, xmax] = optargs{:};
    
    % finding sweep numbers from file name
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    samplingFrequency = obj.header.Acquisition.SampleRate;
    
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),") - niceplotBandpass ", colorName));
        [x,y] = obj.xy(sweepNumber, 1);
        yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
        
        % choosing plot color based on user input        
        if colorName == 'bg'
            plot(x,yFiltered,'Color', [0/255 158/255 115/255],'LineWidth',1);   % bluish green 
        elseif colorName == 'v'
            plot(x,yFiltered,'Color', [213/255 94/255 0/255],'LineWidth',1);   % vermillion
        elseif colorName == 'rp'
            plot(x,yFiltered,'Color', [204/255 121/255 167/255],'LineWidth',1);   % reddish purple
        elseif colorName == 'sb'
            plot(x,yFiltered,'Color', [86/255 180/255 233/255],'LineWidth',1);   % sky blue
        elseif colorName == 'o'
            plot(x,yFiltered,'Color', [230/255 159/255 0/255],'LineWidth',1);   % orange
        else 
            plot(x,yFiltered,'k','LineWidth',1);   % black
        end

        axis([xmin xmax ymin ymax+10])
%         xlabel('Time (s)');
%         ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        set(gca,'Visible','off')
        
        % adding light stim
        % green line
%         line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0.4660, 0.6740, 0.1880],'LineWidth',10)
%         line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0, 200/255, 0],'LineWidth',10)

        % blue line
%         line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
        line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0, 114/255, 178/255],'LineWidth',10)

        % adding scale bar
        line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    
        
    end             
end