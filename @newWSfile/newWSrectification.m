function newWSrectification(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);
optargs = {0.015 10 2.9 2.95 3.2 -250 250 0.1 'D:\CORONAVIRUS DATA\From MATLAB'}; 
optargs(1:numvarargs) = varargin;
[analysisTimeWindow, threshold, baselineOnsetTime, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
sampleRate = obj.header.AcquisitionSampleRate;
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);

% matrixes that will be filled
allPeakCurrent = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
data = [];
allRs = [];
allOnsetLatencyInMilliSeconds = [];
allPeakLatencyInMilliSeconds = []; 
dataPerSweepCh1abs = [];
allHoldingV = [];
allPeakCurrentAbs = [];

% checking for incomplete sweeps for obj
% and not analyzing incomplete sweeps (to
% avoid this error: "Dimensions of matrices being concatenated are not
% consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
% ych2]" and other concatenations.       
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
    lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 2;
    allSweeps = firstSweepNumber:lastSweepNumber;
end

    
%% getting all the data for obj
for sweepNumber = allSweeps

    % getting light stim data from ch2
    [xch2,ych2] = obj.xy(sweepNumber, 2);     
    dataPerSweepCh2 = [dataPerSweepCh2, ych2];

    % finding info about light stim        
    lightPulseStartList = find(diff(ych2>1)>0);
    lightPulseEndList = find(diff(ych2<1)>0);
    lightPulseOnsetDataPoint = lightPulseStartList(1);
    lightPulseOnsetTime = lightPulseStartList(1)/sampleRate;                           % in seconds
    stimDur = (lightPulseEndList(end)-lightPulseStartList(end))/sampleRate;            % duration of the light pulse (s)
    afterLightDataPoint = (lightPulseOnsetTime+analysisTimeWindow)*sampleRate;

    % getting data from ch1
    [x,y] = obj.xy(sweepNumber, 1);
    
    % getting holding voltage from ch3
    % TL started holding the cell at -60 mV
    [xch3,ych3] = obj.xy(sweepNumber, 3);
    holdingV = -60 + round(ych3(baselineOnsetDataPoint));

    % calculating series resistance
    rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
    rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
    dCurrent = rsTransientCurrent-rsBaselineCurrent;
    dVoltage = -5;
    seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm      

    % normalizing data based on mean current 10 ms before baselineOnsetDataPoint
    ynorm = y-mean(y(baselineOnsetDataPoint-10:baselineOnsetDataPoint));
    
    % getting absolute values to calculate info about peaks irregardless of
    % their direction (inward or outward)
    yabs = abs(ynorm);

    % calculating latency 
    % rationale for finding light evoked response latency: I am looking
    % for a monotonic change of the signal for at least "x" data points 
    % (x=threshold), and I am selecting the first occurence of this monotonic change.

    % the expected direction of the monotonic change is determined
    % by the optarg "inwardORoutward". If the value is 1, the function
    % will look for a monotonic increase (outward current), and if the
    % value is -1, the function will look for a monotonic decrease
    % (inward current).

    % y(lightOnsetDataPoint:lightOffDataPoint) is the signal during light pulse

    % zeros(threshold,1)+(inwardORoutward) is a vector with "threshold" 
    % number of rows containing -1 or 1, depending on the value of the
    % optarg "inwardORoutward"

    % function diff calculates difference between data points

    % function sign only keeps info about decrease (-1) or increase (1) in signal

    % function conv2 performs convolution, looking for a sequence of
    % "threshold" points of value equal to the value stored in
    % "inwardORoutward" (-1 or 1). I think conv2 collapses all found 
    % elements into a single element of value "threshold" and moves 
    % forward in its search. So if threshold=20, inwardORoutward=-1,
    % and there are 40 elements of value -1 in sequence (one after the
    % other), conv2 will spit out a vector with 2 elements, both of
    % value 20.

    % function find looks for the index of the elements equal to threshold in the output of the convolution

    % function min looks for the min index (aka the first data point
    % that marks the start of a monotonic change in signal for "x"
    % points (x=threshold).
    
    % UPDATE for this particular code: I am using yabs, so all peaks will
    % look outward (aka inwardORoutward = 1).
    inwardORoutward = 1;
    onsetLatencyInDataPoints = min(find(conv2(sign(diff(yabs(lightPulseOnsetDataPoint:afterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
    
    % 1000 multiplication is conversion from seconds to milliseconds
    onsetLatencyInMilliSeconds = 1000*onsetLatencyInDataPoints/sampleRate;
    
    % to avoid errors due to empty matrices
    if isempty(onsetLatencyInMilliSeconds)
        onsetLatencyInMilliSeconds = NaN;
    end

    % storing normalized data
    dataPerSweepCh1 = [dataPerSweepCh1, ynorm];
    dataPerSweepCh1abs = [dataPerSweepCh1abs, yabs];

    % storing relevant info
    % first: find current peaks using yabs (using the absolute values of y
    % allows me to find peak inward and outward currents)
    peakCurrentAbs = max(yabs(lightPulseOnsetDataPoint:afterLightDataPoint));
    % now find the location (index) of the max current
    peakCurrentLocs = find(yabs(lightPulseOnsetDataPoint:afterLightDataPoint) == peakCurrentAbs);  % data points
    % sometimes there will be more than one location where yabs=max(yabs),
    % so let's keep only the first one:
    peakCurrentLoc = peakCurrentLocs(1,1);                                                          % 1st data point
    % now find the peak value in y (we'll get a negative value for inward
    % current and a positive value for outward currents). Note that the
    % index found above is much smaller than the actual index needed to
    % find the correct peak inside ynorm. In a previous attempt, I used
    % this code "peakCurrent = ynorm(peakCurrentLoc)" and it gave me
    % bizarre results. You have to add "lightPulseOnsetDataPoint" to 
    % "peakCurrentLoc" to get the right index!!
    peakCurrent = ynorm(lightPulseOnsetDataPoint + peakCurrentLoc);
    peakLatencyInMilliSeconds = 1000*peakCurrentLoc/sampleRate;               

    % create table that will be saved as csv file      
    data = [data;...
        mouseNumber, experimentDate, sweepNumber, seriesResistance,...
        holdingV, stimDur, peakCurrent,...
        onsetLatencyInMilliSeconds, peakLatencyInMilliSeconds];

    % create lists that will be used to plot sweep x Rs, sweep x latency
    % and IV curves (current x voltage)
    allRs = [allRs, seriesResistance];
    allOnsetLatencyInMilliSeconds = [allOnsetLatencyInMilliSeconds, onsetLatencyInMilliSeconds];
    allPeakLatencyInMilliSeconds = [allPeakLatencyInMilliSeconds, peakLatencyInMilliSeconds];  
    allPeakCurrent = [allPeakCurrent, peakCurrent];
    allPeakCurrentAbs = [allPeakCurrentAbs, peakCurrentAbs];
    allHoldingV = [allHoldingV, holdingV];
   
end


%% plotting

% dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
% calculates the average column. mean(dataPerSweepCh1) calculates the
% wrong thing: the average row.   

% plotting found peaks, raw data and averages
    figure('name', strcat(obj.file,' - AMPA rectification')); 
    hold on;
    % Talia's protocol tests 6 holding currents (i=6) 3 times.
    for i=1:6
        hold on;
        % plot raw data from each sweep
        plot(x, dataPerSweepCh1(:,i),'Color', [0 0 0 0.2]);
        plot(x, dataPerSweepCh1(:,i+6),'Color', [0 0 0 0.2]);
        plot(x, dataPerSweepCh1(:,i+12),'Color', [0 0 0 0.2]);
        % plot mean current per holding step
        plot(x, mean([dataPerSweepCh1(:,i),...
            dataPerSweepCh1(:,i+6),...
            dataPerSweepCh1(:,i+12)],2), 'Color',[0 0 0],'LineWidth',1.5);
    end
    plot(lightPulseOnsetTime + allPeakLatencyInMilliSeconds/1000, allPeakCurrent, 'o', 'color', 'red');
    rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    hold off;
    axis([xmin xmax ymin ymax]);
    xlabel('Time (s)');
    ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
    title([obj.file ' AMPA rectification'],'Interpreter','none');
    movegui('northwest');
    
% plotting current at -60, 0 and +40 mV - niceplot
    figure('name', strcat(obj.file,' - AMPA rectification niceplot')); 
    hold on;
    % Talia's protocol tests 6 holding currents (i=6) 3 times.
    for i=[1 4 6]
        hold on;
        % plot raw data from each sweep
        plot(x, dataPerSweepCh1(:,i),'Color', [0 0 0 0.2]);
        plot(x, dataPerSweepCh1(:,i+6),'Color', [0 0 0 0.2]);
        plot(x, dataPerSweepCh1(:,i+12),'Color', [0 0 0 0.2]);
        % plot mean current per holding step
        plot(x, mean([dataPerSweepCh1(:,i),...
            dataPerSweepCh1(:,i+6),...
            dataPerSweepCh1(:,i+12)],2), 'Color',[0 0 0],'LineWidth',1.5);
    end
    rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    axis([xmin xmax ymin ymax]);
    set(gca,'Visible','off')
    % adding scale bar
    line([xmax-(xmax-xmin)/10,xmax],[ymax,ymax],'Color','k')
    line([xmax,xmax],[ymax,ymax-((ymax-ymin)/10)],'Color','k')
    text(xmax-(xmax-xmin)/10,ymax-((ymax-ymin)/25),strcat(num2str(1000*(xmax-xmin)/10)," ms"))
    text(xmax-(xmax-xmin)/8,ymax-((ymax-ymin)/10),strcat(num2str((ymax-ymin)/10)," ",obj.header.AIChannelUnits(1)))     
    hold off;
    movegui('southeast');

% plotting IV curve
    figure('name', strcat(obj.file,' - AMPA rectification IV'));     
    plot(allHoldingV, allPeakCurrent, 'o')
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
    xlabel('mV');
    ylabel('pA');
    axis([-80 60 ymin ymax]);
    title([obj.file ' AMPA rectification IV curve'],'Interpreter','none');
    
%  % plotting found peaks, raw data and averages - used for troubleshooting
%     figure('name', strcat(obj.file,' - AMPA rectification ABS')); 
%     hold on;
%     % Talia's protocol tests 6 holding currents (i=6) 3 times.
%     for i=1:6
%         hold on;
%         % plot raw data from each sweep
%         plot(x, dataPerSweepCh1abs(:,i),'Color', [0 0 0 0.1]);
%         plot(x, dataPerSweepCh1abs(:,i+6),'Color', [0 0 0 0.1]);
%         plot(x, dataPerSweepCh1abs(:,i+12),'Color', [0 0 0 0.1]);
%         % plot mean current per holding step
%         plot(x, mean([dataPerSweepCh1abs(:,i),...
%             dataPerSweepCh1abs(:,i+6),...
%             dataPerSweepCh1abs(:,i+12)],2), 'Color',[0 0 0],'LineWidth',1.5);
%     end
%     plot(lightPulseOnsetTime + allPeakLatencyInMilliSeconds/1000, allPeakCurrentAbs, 'o', 'color', 'red');
%     rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%     hold off;
%     axis([xmin xmax ymin ymax]);
%     xlabel('Time (s)');
%     ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%     title([obj.file ' AMPA rectification'],'Interpreter','none');
    
% % plotting niceplot
%     % naming figure file
%     figure('name', strcat(obj.file,' (all) - niceplot - Ampa rectification')); 
%     hold on;
%     plot(x, dataPerSweepCh1,'Color', [0 0 0 0.05]);
%     plot(x, dataPerSweepCh1Ampa,'Color', [213/255 94/255 0/255 0.05]);
%     plot(x, mean(dataPerSweepCh1,2),'Color',[0 0 0],'LineWidth',1.5); 
%     plot(x, mean(dataPerSweepCh1Ampa,2),'Color',[213/255 94/255 0/255],'LineWidth',1.5);
%     rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%     axis([xmin xmax ymin ymax]);
%     xlabel('Time (s)');
%     ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%     set(gca,'Visible','off')
%     % adding scale bar
%     line([xmax-(xmax-xmin)/10,xmax],[ymax,ymax],'Color','k')
%     line([xmax,xmax],[ymax,ymax-((ymax-ymin)/10)],'Color','k')
%     text(xmax-(xmax-xmin)/10,ymax-((ymax-ymin)/25),strcat(num2str(1000*(xmax-xmin)/10)," ms"))
%     text(xmax-(xmax-xmin)/8,ymax-((ymax-ymin)/10),strcat(num2str((ymax-ymin)/10)," ",obj.header.AIChannelUnits(1)))     
%     hold off;
    
% % plot data and found peaks
%     % naming figure file
%     figure('name', strcat(obj.file,' - Baseline Subtracted - Total and Ampa')); 
%     hold on;   
%     % plotting currents
%     plot(x, dataPerSweepCh1(:,:),'color',[0 0 0 0.2]);
%     plot(x, dataPerSweepCh1Ampa(:,:),'color',[213/255 94/255 0/255 0.2]); 
%     plot(lightPulseOnsetTime + data(:,8)/1000, data(:,6),'o','color','red');
%     plot(lightPulseOnsetTime + dataAmpa(:,8)/1000, dataAmpa(:,6),'o','color','red');
%     % plotting light stim
%     rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none');  
%     axis([xmin xmax ymin ymax]);
%     xlabel('Time (s)');
%     ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%     title([obj.file ' Total and AMPA'],'Interpreter','none');
%     hold off;
%     movegui('northwest');
    
% plot rs 
    figure('name', strcat(obj.file,' (all) - rs - AMPA rectification')); % naming figure file
    plot(allSweeps, allRs,'-o');
    % plot lines marking 30% increase and 30% decrese in Rs compared to first
    % test pulse
    line([allSweeps(1) allSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
    line([allSweeps(1) allSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
    axis([allSweeps(1) inf 0 60])
    ylabel('Rs (M\Omega)');
    xlabel('Sweeps');
    title([obj.file ' rs'],'Interpreter','none');
    movegui('northeast');
          
    
% save xls file with all the data
    filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - AMPA rectification",'(',num2str(inwardORoutward),')');
    fulldirectory = strcat(savefileto,'\',filename,'.xls');        
    dataInCellFormat = {};
    dataInCellFormat = num2cell(data);
    labeledData = cell2table(dataInCellFormat,'VariableNames',...
        {'mouse', 'date', 'sweep', 'seriesResistance',...
        'holdingVinMilliV', 'lightStimDur', 'peakCurrent',...
        'onsetLatencyInMilliSeconds', 'peakLatencyInMilliSeconds'});
    writetable(labeledData,fulldirectory);
    disp('I saved it')

end