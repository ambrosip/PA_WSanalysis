function psallON(obj, varargin)

    % optional arguments
    % set defaults for optional inputs 
    %======= NOTE CAPPING OF SIGNAL AT 12 HZ: 1/12 = 0.083
    optargs = {10 10 3 15 0.083 12 inf 10 1 1};
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [window, MinPeakHeight, LightDur, LightOnsetTime, MinPeakDistance, ymaxhist, ymax, NumPeaksBeforeAfter, DelayScaleFactor, LightExtensionFactor] = optargs{:};

    % finding sweep numbers from file name
    if length(obj.file) == 28
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;
    
    % plotting one figure per sweep that shows firing rate histogram for the whole sweep
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - histogram')); % naming figure file
        subplot(2,1,1)
        [x,y] = obj.xy(sweepNumber, 1);
y = smoothdata(y,'gaussian',window);
        [yupper, ylower] = envelope(y);     
        [pks,locs,w,p] = findpeaks(yupper,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        
        sweepDuration = obj.header.Acquisition.Duration;
        sweepTime=0;
        inverseISIperSecBin=[];
                   
        while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
            indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
            locsDuringSweepTime = locs;
            locsDuringSweepTime(indicesToDelete) = [];
            inverseISIperSecBin=[inverseISIperSecBin 1/mean(diff(locsDuringSweepTime))];
            sweepTime=sweepTime+1;
        end

        bar(0:29,inverseISIperSecBin,1);
        axis([-inf inf 0 ymaxhist]);
            
        xlabel('Bins (1 s long)');
%         ylabel('Firing Frequency (Mean 1/ISI)');
        ylabel('1/ISI (Hz)');
        title([' (' num2str(sweepNumber) ') - hist'],'Interpreter','none');
%         title([obj.file ' (' num2str(sweepNumber) ') - histogram'],'Interpreter','none');
%         set(gcf,'Position',[1 1 560 420]) % default figure size in Matlab
        
        subplot(2,1,2)
        [x,y] = obj.xy(sweepNumber, 2);     
        plot(x,y);
        yminhere = min(y)-5;
        ymaxhere = max(y)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
        
        set(gcf,'Position',[1 1 280 420])
        movegui('northwest');
        
%         inverseISIperSecBin
%         locs
        
    end
    
    
    % plotting one figure per sweep that shows Ch1, Ch2, peaks, and 1/ISI    
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - peaks')); % naming figure file
        [x,y] = obj.xy(sweepNumber, 1);
y = smoothdata(y,'gaussian',window);
        [yupper, ylower] = envelope(y);          
        [pks,locs,w,p] = findpeaks(yupper,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        
        % store locs and pks values in new variables that will be editted
        % acordingly
        locsLight = locs;
        locsBaseline1 = locs;
        locsBaseline2 = locs;
        pksLight = pks;
        pksBaseline1 = pks;
        pksBaseline2 = pks;
        
        % find indices at which the value of locs is between the light
        % onset and offset (with optional extension factor to look for
        % lingering light effect)
        indicesLight = find(locs<LightOnsetTime | locs>(LightOnsetTime+LightDur*LightExtensionFactor));
        
        % deletes all peaks that happen before or after the light pulse
        locsLight(indicesLight) = [];
        pksLight(indicesLight) = [];
        
        % find indices for baseline prior to light stim and delete unwanted
        % peaks
        indicesBaseline1 = find(locs>=LightOnsetTime);
        locsBaseline1(indicesBaseline1) = [];
        pksBaseline1(indicesBaseline1) = [];
        locsBaseline1(1:length(locsBaseline1)-NumPeaksBeforeAfter)=[];
        pksBaseline1(1:length(pksBaseline1)-NumPeaksBeforeAfter)=[];
        
        % find indices for baseline after light stim and delete unwanted
        % peaks
        indicesBaseline2 = find(locs<=(LightOnsetTime+LightDur*DelayScaleFactor));
        locsBaseline2(indicesBaseline2) = [];
        pksBaseline2(indicesBaseline2) = [];
        locsBaseline2(NumPeaksBeforeAfter+1:length(locsBaseline2))=[];
        pksBaseline2(NumPeaksBeforeAfter+1:length(pksBaseline2))=[];
        
        % calculate interspike interval (ISI) and firing frequency as 1/ISI
        % (Hz)
        inverseISIforLight = 1/mean(diff(locsLight));
        inverseISIforBaseline1 = 1/mean(diff(locsBaseline1));
        inverseISIforBaseline2 = 1/mean(diff(locsBaseline2));
     
        % Ch1 subplot (recorded data)
        subplot(2,1,1);
        % OLD: subplot(2,lastSweepNumber-firstSweepNumber+1,sweepNumber-firstSweepNumber+1);
        plot(x,yupper)
        hold on;
        plot(locsLight,pksLight,'o','color','red');
        plot(locsBaseline1,pksBaseline1,'o','color','blue');
        plot(locsBaseline2,pksBaseline2,'o','color','blue');
%         rectangle('Position',[LightOnsetTime,-85,LightDur,145],'FaceColor', [0 0 1 0.05],'LineStyle','none');
        hold off;
%         ymax = 1.1*max(yupper);
        axis([LightOnsetTime-4 LightOnsetTime+LightDur+4 -inf ymax])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ' (env)'));
        title([obj.file ' (' num2str(sweepNumber) ') - peaks'],'Interpreter','none');
%         text(LightOnsetTime + LightDur/4,max(yupper),num2str(round(inverseISIforLight,2)),'color','red');
%         text(LightOnsetTime - 3,max(yupper),num2str(round(inverseISIforBaseline1,2)),'color','blue');
%         text(LightOnsetTime + LightDur + 2,max(yupper),num2str(round(inverseISIforBaseline2,2)),'color','blue');
        text(LightOnsetTime + LightDur/4,0.9*ymax,num2str(round(inverseISIforLight,2)),'color','red');
        text(LightOnsetTime - 3,0.9*ymax,num2str(round(inverseISIforBaseline1,2)),'color','blue');
        text(LightOnsetTime + LightDur + 2,0.9*ymax,num2str(round(inverseISIforBaseline2,2)),'color','blue');
        
        % Ch2 subplot (light stim command)
        subplot(2,1,2);
        [x,y] = obj.xy(sweepNumber, 2);
        plot(x,y)
        yminhere = min(y)-5;
        ymaxhere = max(y)+5;
        axis([LightOnsetTime-4 LightOnsetTime+LightDur+4 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        
%         set(gcf,'Position',[100 200 1750 375]) % for 3 plots        
%         locsLight
    end
    
    
    % plotting one figure per sweep that shows whole sweep
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - raw')); % naming figure file
        
        subplot(2,1,1)
        [x,y] = obj.xy(sweepNumber, 1);
y = smoothdata(y,'gaussian',window);
        plot(x,y);
        axis([-inf inf -inf ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
      
        subplot(2,1,2)
        [x,y] = obj.xy(sweepNumber, 2);     
        plot(x,y);
        yminhere = min(y)-5;
        ymaxhere = max(y)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
        
        set(gcf,'Position',[1 1 280 420])
        movegui('northeast');
    end
    
end