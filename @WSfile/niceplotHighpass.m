function niceplotHighpass(obj, varargin)

    % highpass(data,0.5,10000,'Steepness', 0.5,'StopbandAttenuation', 10)

    % optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
    optargs = {'k' 0.5 0.5 10 -500 100 10 20};          % for mIPSC data
    optargs(1:numvarargs) = varargin;
    [colorName, passbandFrequency, steepness, stopbandAttenuation, ymin, ymax, xmin, xmax] = optargs{:};
    
    % finding sweep numbers from file name
    [~, ~, allSweeps] = getSweepNumbers(obj);

    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    samplingFrequency = obj.header.Acquisition.SampleRate;
    
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),") - niceplotHighpass ", colorName));
        [x,y] = obj.xy(sweepNumber, 1);
        yFiltered = highpass(y, passbandFrequency, samplingFrequency, 'Steepness', steepness, 'StopbandAttenuation', stopbandAttenuation);
        
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

        % adding scale bar
        line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))           
    
    end             
end


