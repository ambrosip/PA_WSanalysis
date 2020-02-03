function lightvsfiringONbandpass(obj, varargin)

    % optional arguments
    % set defaults for optional inputs 
    optargs = {100 400 10 30 3 15 0.01 50 1 1 'E:\From MATLAB'};
    
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    
    % place optional args in memorable variable names
    [highpassThreshold, lowpassThreshold, MinPeakHeight, ymax, LightDur, LightOnsetTime, MinPeakDistance, ymaxhist, DelayScaleFactor, LightExtensionFactor, savefileto] = optargs{:};

    data = [];
    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    samplingFrequency = obj.header.Acquisition.SampleRate;

    % finding sweep numbers from file name
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);
    

    for sweepNumber = allSweeps
        
        [x,y] = obj.xy(sweepNumber, 1);
        y = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
        [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        [xch2,ych2] = obj.xy(sweepNumber, 2);  
        
        sweepDuration = obj.header.Acquisition.Duration;
        sweepTime=0;
        inverseISIperSecBin=[];
        
        % plotting one figure per sweep that shows firing rate histogram for the whole sweep       
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - histogram')); % naming figure file
        subplot(2,1,1)
         
        while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
            indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
            locsDuringSweepTime = locs;
            locsDuringSweepTime(indicesToDelete) = [];
            inverseISIperSecBin=[inverseISIperSecBin 1/mean(diff(locsDuringSweepTime))];
            sweepTime=sweepTime+1;
        end

        bar(0:length(inverseISIperSecBin)-1,inverseISIperSecBin,1);
        axis([-inf inf 0 ymaxhist]);            
        xlabel('Bins (1 s long)');
        ylabel('1/ISI (Hz)');
        title([' (' num2str(sweepNumber) ') - hist'],'Interpreter','none');        
        
        subplot(2,1,2)   
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        set(gcf,'Position',[1 1 280 420])
        movegui('northwest');   
            
        % plotting one figure per sweep that shows Ch1, Ch2, peaks, and 1/ISI     
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - peaks')); % naming figure file
        
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
        indicesBaseline1 = find(locs<LightOnsetTime-LightDur | locs>=LightOnsetTime);
        locsBaseline1(indicesBaseline1) = [];
        pksBaseline1(indicesBaseline1) = [];
        
        % find indices for baseline after light stim and delete unwanted
        % peaks
        indicesBaseline2 = find(locs<(LightOnsetTime+LightDur*LightExtensionFactor) | locs>(LightOnsetTime+(LightDur*LightExtensionFactor)+LightDur));
        locsBaseline2(indicesBaseline2) = [];
        pksBaseline2(indicesBaseline2) = [];
        
        % calculate interspike interval (ISI) and firing frequency as 1/ISI
        % (Hz)
        ISIforLight = (diff(locsLight));
        ISIforBaseline1 = (diff(locsBaseline1));
        ISIforBaseline2 = (diff(locsBaseline2));
        
        % calculate interspike interval (ISI) and firing frequency as 1/ISI
        % (Hz)
        avgInverseISIforLight = 1/mean(diff(locsLight));
        avgInverseISIforBaseline1 = 1/mean(diff(locsBaseline1));
        avgInverseISIforBaseline2 = 1/mean(diff(locsBaseline2));
     
        % Ch1 subplot (recorded data)
        subplot(2,1,1);
        plot(x,y)
        hold on;
        plot(locsLight,pksLight,'o','color','red');
        plot(locsBaseline1,pksBaseline1,'o','color','blue');
        plot(locsBaseline2,pksBaseline2,'o','color','blue');
        hold off;
        axis([LightOnsetTime-4 LightOnsetTime+LightDur+4 -inf ymax])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ' (env)'));
        title([obj.file ' (' num2str(sweepNumber) ') - peaks'],'Interpreter','none');
        text(LightOnsetTime + LightDur/4,0.9*ymax,num2str(round(avgInverseISIforLight,2)),'color','red');
        text(LightOnsetTime - 3,0.9*ymax,num2str(round(avgInverseISIforBaseline1,2)),'color','blue');
        text(LightOnsetTime + LightDur + 2,0.9*ymax,num2str(round(avgInverseISIforBaseline2,2)),'color','blue');
        
        % Ch2 subplot (light stim command)
        subplot(2,1,2);
        plot(xch2,ych2)
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([LightOnsetTime-4 LightOnsetTime+LightDur+4 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));  
                   
        placeholder1 = [];
        placeholder2 = [];
        placeholder3 = [];
        placeholder4 = [];
        placeholder5 = [];
        placeholder6 = [];
        
                if isempty(ISIforBaseline1)
                    placeholder1 = NaN;
                    placeholder2 = NaN;
                else
                    placeholder1 = 1/ISIforBaseline1(1);
                    placeholder2 = 1/ISIforBaseline1(end);
                end
                
               if isempty(ISIforLight)
                    placeholder3 = NaN;
                    placeholder4 = NaN;
                else
                    placeholder3 = 1/ISIforLight(1);
                    placeholder4 = 1/ISIforLight(end);
                 end
                
                if isempty(ISIforBaseline2)
                    placeholder5 = NaN;
                    placeholder6 = NaN;
                else
                    placeholder5 = 1/ISIforBaseline2(1);
                    placeholder6 = 1/ISIforBaseline2(end);
                end
                
        data = [data; mouseNumber, experimentDate, sweepNumber, ...
                avgInverseISIforBaseline1, avgInverseISIforLight, avgInverseISIforBaseline2, ...
                length(pksBaseline1), length(pksLight), length(pksBaseline2),...
                placeholder1, placeholder2,...
                placeholder3, placeholder4,...
                placeholder5, placeholder6];
                        
        % plotting one figure per sweep that shows whole sweep    
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - raw')); % naming figure file    
        
        subplot(2,1,1)
        plot(x,y);
        axis([-inf inf -inf ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
      
        subplot(2,1,2) 
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
        
        set(gcf,'Position',[1 1 280 420])
        movegui('northeast');

        % save csv file with data 
        filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - firing vs light cell attached");
        fulldirectory = strcat(savefileto,'\',filename,'.csv');        
        dataInCellFormat = {};
        dataInCellFormat = num2cell(data);
        labeledData = cell2table(dataInCellFormat,'VariableNames',...
            {'mouse', 'date', 'sweep',...
            'invISIpre', 'invISIlight', 'invISIpost',...
            'APsPre', 'APsLight', 'APsPost',...
            'firstInvISIpre','lastInvISIpre',...
            'firstInvISIlight','lastInvISIlight',...
            'firstInvISIpost','lastInvISIpost'});
        writetable(labeledData,fulldirectory);
        disp('I saved it')
        disp('Change directory if you want this saved elsewhere!')     
    
end