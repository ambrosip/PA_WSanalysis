function niceplot(obj, varargin)

    % optional arguments: axis range for channel 1 (monitor CC or VC)
    numvarargs = length(varargin);
%     optargs = {'k' -100 50 0 20 15 3 10 23};     % for WC data
%     optargs = {'k' -200 200 0 20 15 3 10 23};    % for ON data
%     optargs = {'k' -600 100 0 20 15 3 10 23};    % for mIPSC data
    optargs = {'k' -500 100 1 10 15 3 10 23};    % for filtered mIPSC data
    optargs(1:numvarargs) = varargin;
    [colorName, ymin, ymax, filterORnot, GaussianFilterWindow, lightOnsetTime, lightDuration, xmin, xmax] = optargs{:};
    
    % finding sweep numbers from file name
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);
        
    % plotting individual figures for all sweeps    
    for sweepNumber = allSweeps  
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),") - niceplot ", colorName));
        [x,y] = obj.xy(sweepNumber, 1);
        
        % choosing plot color based on user input        
        if colorName == 'bg'
            plot(x,y,'Color', [0/255 158/255 115/255],'LineWidth',1);   % bluish green 
        elseif colorName == 'v'
            plot(x,y,'Color', [213/255 94/255 0/255],'LineWidth',1);   % vermillion
        elseif colorName == 'rp'
            plot(x,y,'Color', [204/255 121/255 167/255],'LineWidth',1);   % reddish purple
        elseif colorName == 'sb'
            plot(x,y,'Color', [86/255 180/255 233/255],'LineWidth',1);   % sky blue
        elseif colorName == 'o'
            plot(x,y,'Color', [230/255 159/255 0/255],'LineWidth',1);   % orange
        else 
            plot(x,y,'k','LineWidth',1);   % black
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
%         line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0, 114/255, 178/255],'LineWidth',10)

        % adding scale bar
        line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))

        % most likely useless code now that I have niceplotBandpass
        if filterORnot
            yfiltered = smoothdata(y,'gaussian',GaussianFilterWindow);
            figure('name', strcat(obj.file,' (',num2str(sweepNumber),') - niceplot filtered'));
            plot(x,yfiltered,'k','LineWidth',1);
            axis([xmin xmax ymin ymax+10])
            title([obj.file ' (' num2str(sweepNumber) ") filtered"],'Interpreter','none');
            set(gca,'Visible','off')
            % adding light stim line
%             line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0.4660, 0.6740, 0.1880],'LineWidth',10)
%             line([lightOnsetTime,lightOnsetTime+lightDuration],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
            line([xmin,xmin+(xmax-xmin)/13],[ymin,ymin],'Color','k')
            line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
            text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/13)," s"))
            text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))        
        end    
    end             
end