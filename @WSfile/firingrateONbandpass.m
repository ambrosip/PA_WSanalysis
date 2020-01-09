function firingrateONbandpass(obj, varargin)

    %% BEFORE USING THIS FUNCTION
    
    % Add the folder where your h5 files to MATLAB's path
    
    % Load the h5 files using the fuction WSfileFromDir. Name the object
    % using the mouse number and make sure there are 3 numbers. If the
    % mouse number is m23, name it as m023.   
    % EX: m155=WSfileFromDir('E:\Priscilla\Ephys\20190808 dms chr2 habit')
    % TIP: don't load the files directly from the server because that takes
    % a million years. 
    
    % Choose which files you want to analyze - you will have to tell matlab
    % what is the first sweep number of the file you need analyzed.
    
    % Run the function using the optional arguments (or not) - feel free to 
    % edit the default values of the optional arguments to speed up your analysis
    % EX: firingrate(m155.s0201)
    % EX: firingrate(m155.s0201,10,1,5) - this one sets an optional
    % threshold and start-end points for firing rate analysis

    %% BEHOLD THE CODE
    
    % optional arguments
    % set defaults for optional inputs 
%     optargs = {100 500 15 0 30 12 -100 100 0.01 'D:\Temp\From MATLAB 2020'};
%     optargs = {100 500 3 0 30 12 -30 30 0.08 'D:\Temp\From MATLAB 2020'};
    optargs = {100 400 5 0 30 12 -30 30 0.08 'D:\Temp\From MATLAB 2020'};
%     optargs = {100 500 0.4 0 30 12 -1 1 0.05 'D:\Temp\From MATLAB 2020'};
    
    % NOTE that setting MinPeakDistance to 0.08 caps the firing rate at 12
    % Hz (1/12 = 0.083), so change that if you expect your cell to be
    % firing at higher frequencies - this was added to dismiss artifact in
    % VClamp.
    
    % overwrite defaults with values spclose ecified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [highpassThreshold, lowpassThreshold, MinPeakHeight, timeToStartCountingAPs, timeToStopCountingAPs, ymaxhist, ymin, ymax, MinPeakDistance, savefileto] = optargs{:};
    
    data = [];
    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    samplingFrequency = obj.header.Acquisition.SampleRate;

    % finding sweep numbers from file name
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    for sweepNumber = allSweeps
        
    % plotting one figure per sweep that shows firing rate histogram (1s bins) for the whole sweep
    
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - histogram bandpass')); % naming figure file
        [x,y] = obj.xy(sweepNumber, 1);
        yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
        [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        
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
        ylabel('1/ISI (Hz)');
        title([' (' num2str(sweepNumber) ') - hist'],'Interpreter','none');              
        set(gcf,'Position',[1 1 280 420])
        movegui('northwest');
        
        
    % plotting one figure per sweep that shows whole sweep filtered
    
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - bandpass filtered')); % naming figure file
        plot(x,yFiltered);
        axis([-inf inf ymin ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - filtered'],'Interpreter','none');     
        set(gcf,'Position',[1 1 280 420])
        movegui('northeast');
        
        
    % plotting one figure per sweep that shows whole sweep raw
    
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - raw')); % naming figure file
        plot(x,y);
        axis([-inf inf -inf inf]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
        % adding line at 0 pA
        line([0 30],[0,0],'Color','black','LineStyle','--')
        set(gcf,'Position',[1 1 280 420])
        movegui('southeast');
        
        
    % plotting one figure per sweep that shows found peaks 
    
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - peaks bandpass filtered')); % naming figure file
       
        % store locs and pks values in new variables that will be editted
        % accordingly
        locsToPlot = locs;
        pksToPlot = pks;
        
        % find indices at which the value of locs is beyond the desired
        % time interval
        indicesToDelete = find(locs<timeToStartCountingAPs | locs>timeToStopCountingAPs);
        
        % deletes all peaks that happen before or after the desired
        % interval
        locsToPlot(indicesToDelete) = [];
        pksToPlot(indicesToDelete) = [];
                
        % calculate interspike interval (ISI) and firing frequency as 1/ISI
        % (Hz)
        ISItoPlot = (diff(locsToPlot));
        avgInverseISI = 1/mean(diff(locsToPlot));
     
        % plot
        plot(x,yFiltered)
        hold on;
        plot(locsToPlot,pksToPlot,'o','color','red');
        hold off;
        axis([timeToStartCountingAPs-1 timeToStopCountingAPs+1 ymin ymax])
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([obj.file ' (' num2str(sweepNumber) ') - peaks bandpass filtered'],'Interpreter','none');
        text(timeToStartCountingAPs+(timeToStopCountingAPs-timeToStartCountingAPs)/2,-75,num2str(round(avgInverseISI,2)),'color','red');

        
    % creating csv file 
    
        placeholder1 = [];
        placeholder2 = [];
        
        if isempty(ISItoPlot)
            placeholder1 = NaN;
            placeholder2 = NaN;
        else
            placeholder1 = 1/ISItoPlot(1);
            placeholder2 = 1/ISItoPlot(end);
        end

        timeInterval = timeToStopCountingAPs-timeToStartCountingAPs;
        nAPsperTime = length(pksToPlot)/timeInterval;
               
        data = [data; 
                mouseNumber, ...
                experimentDate, ...
                sweepNumber, ...
                timeInterval,...
                avgInverseISI, ...
                nAPsperTime,...
                length(pksToPlot),...
                highpassThreshold,...
                lowpassThreshold,...
                placeholder1, placeholder2];        
       
    end
    
    
    % save csv file with data
    
    filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - firing");
    fulldirectory = strcat(savefileto,'\',filename,'.csv');        
    dataInCellFormat = {};
    dataInCellFormat = num2cell(data);
    labeledData = cell2table(dataInCellFormat,'VariableNames',...
        {'mouse',...
        'date',...
        'sweep',...
        'timeInterval',...
        'avgInvISI',...
        'nAPsperTime',...
        'numberOfAPs',...
        'highpassThreshold',...
        'lowpassThreshold',...
        'firstInvISI','lastInvISI'});
    writetable(labeledData,fulldirectory);
    disp('I saved the csv file!')
    disp('Change directory if you want this saved elsewhere!') 

    saveAllFigs(savefileto);
    disp('I saved the figures!')

end