% Code adapted from newWSnormmono and lightvsfiringONbandpassAUTO
% Use this code for paired pulse ratio (PPR) analysis

function newWSppr(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);
optargs = {-1 0.01 5 0.2 0.2 0.9 -500 50 0.1 'D:\CORONAVIRUS DATA\From MATLAB'};  
optargs(1:numvarargs) = varargin;
[inwardORoutward, analysisTimeWindow, threshold, baselineOnsetTime, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
totalActiveChannels = 2;
sampleRate = obj.header.AcquisitionSampleRate;
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);

% matrixes that will be filled
allFirstLightEvokedCurrent = [];
allSecondLightEvokedCurrent = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
data = [];
allRs = [];
allFirstLightEvokedResponseLatencyInMilliSeconds = []; 
allSecondLightEvokedResponseLatencyInMilliSeconds = []; 
allFirstPeakLatencyInMilliSeconds = [];  
allSecondPeakLatencyInMilliSeconds = []; 

% checking for incomplete sweeps and not analyzing incomplete sweeps (to
% avoid this error: "Dimensions of matrices being concatenated are not
% consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
% ych2]" and other concatenations.       
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
    lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 2;
    allSweeps = firstSweepNumber:lastSweepNumber;
end
    
% getting all the data
for sweepNumber = allSweeps

    % getting light stim data from ch2
    [xch2,ych2] = obj.xy(sweepNumber, 2);     
    dataPerSweepCh2 = [dataPerSweepCh2, ych2];

    % finding info about light stim        
    lightPulseStartList = find(diff(ych2>1)>0);
    lightPulseEndList = find(diff(ych2<1)>0);
    firstLightPulseOnsetDataPoint = lightPulseStartList(1);
    secondLightPulseOnsetDataPoint = lightPulseStartList(2);
    firstLightPulseOnsetTime = lightPulseStartList(1)/sampleRate;                           % in seconds
    secondLightPulseOnsetTime = lightPulseStartList(2)/sampleRate;
    StimDur = (lightPulseEndList(end)-lightPulseStartList(end))/sampleRate;                 % duration of each light pulse in the train (s)
    StimInterval = (lightPulseStartList(2)-lightPulseStartList(1))/sampleRate - StimDur;    % interval between each pulse (s)
    StimFreq = 1/StimInterval;                                                              % frequency of the light stim (Hz)
    firstAfterLightDataPoint = (firstLightPulseOnsetTime+analysisTimeWindow)*sampleRate;
    secondAfterLightDataPoint = (secondLightPulseOnsetTime+analysisTimeWindow)*sampleRate;

    [x,y] = obj.xy(sweepNumber, 1);

    % calculating series resistance
    rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
    rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
    dCurrent = rsTransientCurrent-rsBaselineCurrent;
    dVoltage = -5;
    seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm      

    % normalizing data based on mean current 10 ms before baselineOnsetDataPoint
    y = y-mean(y(baselineOnsetDataPoint-10:baselineOnsetDataPoint));

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
    firstLightEvokedResponseLatencyInDataPoints = min(find(conv2(sign(diff(y(firstLightPulseOnsetDataPoint:firstAfterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
    secondLightEvokedResponseLatencyInDataPoints = min(find(conv2(sign(diff(y(secondLightPulseOnsetDataPoint:secondAfterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));

    % 1000 multiplication is conversion from seconds to milliseconds
    firstLightEvokedResponseLatencyInMilliSeconds = 1000*firstLightEvokedResponseLatencyInDataPoints/sampleRate;
    secondLightEvokedResponseLatencyInMilliSeconds = 1000*secondLightEvokedResponseLatencyInDataPoints/sampleRate;

    % to avoid errors due to empty matrices
    if isempty(firstLightEvokedResponseLatencyInMilliSeconds)
        firstLightEvokedResponseLatencyInMilliSeconds = NaN;
    end

    if isempty(secondLightEvokedResponseLatencyInMilliSeconds)
        secondLightEvokedResponseLatencyInMilliSeconds = NaN;
    end

    % storing normalized data
    dataPerSweepCh1 = [dataPerSweepCh1, y];

    % storing relevant info
    firstLightEvokedCurrent = min(y(firstLightPulseOnsetDataPoint:firstAfterLightDataPoint));
    firstPeakLatencyDataPoint = find(y(firstLightPulseOnsetDataPoint:firstAfterLightDataPoint)==firstLightEvokedCurrent);
    firstPeakLatencyInMilliSeconds = 1000*firstPeakLatencyDataPoint(1,1)/sampleRate;         
    secondLightEvokedCurrent = min(y(secondLightPulseOnsetDataPoint:secondAfterLightDataPoint));
    secondPeakLatencyDataPoint = find(y(secondLightPulseOnsetDataPoint:secondAfterLightDataPoint)==secondLightEvokedCurrent);
    secondPeakLatencyInMilliSeconds = 1000*secondPeakLatencyDataPoint(1,1)/sampleRate;        
    ppr = secondLightEvokedCurrent/firstLightEvokedCurrent;          

    % not sure I need this?
    allFirstLightEvokedCurrent = [allFirstLightEvokedCurrent, firstLightEvokedCurrent];  
    allSecondLightEvokedCurrent = [allSecondLightEvokedCurrent, secondLightEvokedCurrent]; 

    % create table that will be saved as csv file      
    data = [data; mouseNumber, experimentDate, sweepNumber,...
        seriesResistance, StimDur, StimInterval,...
        firstLightPulseOnsetTime, firstLightEvokedCurrent,...
        firstLightEvokedResponseLatencyInMilliSeconds, firstPeakLatencyInMilliSeconds,...
        secondLightPulseOnsetTime, secondLightEvokedCurrent,...
        secondLightEvokedResponseLatencyInMilliSeconds, secondPeakLatencyInMilliSeconds,...
        ppr];

    % create lists that will be used to plot Rs and latency x sweep
    allRs = [allRs, seriesResistance];
    allFirstLightEvokedResponseLatencyInMilliSeconds = [allFirstLightEvokedResponseLatencyInMilliSeconds, firstLightEvokedResponseLatencyInMilliSeconds];
    allSecondLightEvokedResponseLatencyInMilliSeconds = [allSecondLightEvokedResponseLatencyInMilliSeconds, secondLightEvokedResponseLatencyInMilliSeconds];
    allFirstPeakLatencyInMilliSeconds = [allFirstPeakLatencyInMilliSeconds, firstPeakLatencyInMilliSeconds];  
    allSecondPeakLatencyInMilliSeconds = [allSecondPeakLatencyInMilliSeconds, secondPeakLatencyInMilliSeconds]; 

end

% plotting 
% dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
% calculates the average column. mean(dataPerSweepCh1) calculates the
% wrong thing: the average row.   

% plot 1 figure with raw data for each inter stimulus interval:
    % Talia's protocol tests 6 inter stimulus intervals (i=1:6) 3 times (z=1:3).
    % for all tested stimulus intervals (i=1:6)
    for i=1:6
        % naming figure file
        figure('name', strcat(obj.file,' - Baseline Subtracted - PPR interval ', num2str(1000*data(i,6)))); 
        hold on;

        % for all replications of this stimulus interval (i=1:3)
        for z=1:3
            % plotting currents
            plot(x, dataPerSweepCh1(:,i+6*(z-1)));

            % plotting light stim
            rectangle('Position', [data(i+6*(z-1),7) 25 data(i+6*(z-1),5) 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
            rectangle('Position', [data(i+6*(z-1),11) 25 data(i+6*(z-1),5) 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')

            % plotting min current within LightPulseOnsetTime and LightPulseOnsetTime + analysisTimeWindow
            line([data(i+6*(z-1),7) data(i+6*(z-1),7) + analysisTimeWindow], [data(i+6*(z-1),8), data(i+6*(z-1),8)],'Color','black','LineStyle','--')
            line([data(i+6*(z-1),11) data(i+6*(z-1),11) + analysisTimeWindow], [data(i+6*(z-1),12), data(i+6*(z-1),12)],'Color','black','LineStyle','--')

            % mark peaks using latency
            plot((data(i+6*(z-1),7) + data(i+6*(z-1),10)/1000), data(i+6*(z-1),8),'o','color','red');
            plot((data(i+6*(z-1),11) + data(i+6*(z-1),14)/1000), data(i+6*(z-1),12),'o','color','red');
        end
        axis([xmin xmax ymin ymax])
        hold off;
    end

% plot ppr per stim frequency    
    figure('name', strcat(obj.file,' - PPR')); 
    hold on;
    plot(1000*data(:,6), data(:,15),'o');
    for i=1:6
        plot(1000*data(i,6), mean([data(i,15); data(i+6,15); data(i+12,15)]),'*','color','red')
    end
    line([0 525],[1 1],'Color','black','LineStyle','--')
    axis([0 525 0 2])
    xlabel('Paired Pulse Interval (ms)');
    ylabel('Paired Pulse Ratio');
    title([obj.file ' (all) - PPR'],'Interpreter','none');
    hold off;
    movegui('northwest');
 
%     % plotting figure with only Ch1 niceplot
%         figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted')); % naming figure file
%         hold on;
%         plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
%         plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
%         rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%         axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
%         xlabel('Time (s)');
%         ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%         title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
%         set(gca,'Visible','off')
%         
%         xminScale = lightOnsetTime-0.01;
%         xmaxScale = lightOnsetTime+0.1;
%         
%         % adding scale bar
%         line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
%         line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
%         text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
%         text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.AIChannelUnits(1)))                      
%         hold off;
%         
%     
%     % plotting figure with only Ch1 niceplot - AVG only
%         figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted - AVG only')); % naming figure file
%         hold on;
%         plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
%         rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%         axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
%         xlabel('Time (s)');
%         ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%         title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
%         set(gca,'Visible','off')
%         
%         xminScale = lightOnsetTime-0.01;
%         xmaxScale = lightOnsetTime+0.1;
%         
%         % adding scale bar
%         line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
%         line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
%         text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
%         text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.AIChannelUnits(1)))                      
%         hold off;
 
% save xls file with data 
    filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - PPR",'(',num2str(inwardORoutward),')');
    fulldirectory = strcat(savefileto,'\',filename,'.xls');        
    dataInCellFormat = {};
    dataInCellFormat = num2cell(data);
    labeledData = cell2table(dataInCellFormat,'VariableNames',...
        {'mouse', 'date', 'sweep',...
        'seriesResistance', 'stimDurInSeconds', 'stimIntervalInSeconds',...
        'firstLightPulseOnsetTime', 'firstLightEvokedCurrent',...
        'firstLightEvokedResponseLatencyInMilliSeconds', 'firstPeakLatencyInMilliSeconds',...
        'secondLightPulseOnsetTime', 'secondLightEvokedCurrent',...
        'secondLightEvokedResponseLatencyInMilliSeconds', 'secondPeakLatencyInMilliSeconds',...
        'ppr'});
    writetable(labeledData,fulldirectory);
    disp('I saved it')
    disp('Change directory if you want this saved elsewhere!')     

% plot rs
    figure('name', strcat(obj.file,' (all) - rs')); % naming figure file
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

% plot light response ONSET latency  
    figure('name', strcat(obj.file," (all) - onset latency ", '(',num2str(inwardORoutward),')')); 
    hold on;
    plot(allSweeps, allFirstLightEvokedResponseLatencyInMilliSeconds,'-^');
    plot(allSweeps, allSecondLightEvokedResponseLatencyInMilliSeconds,'-v');
    line([allSweeps(1) allSweeps(end)],[5, 5],'Color','black','LineStyle','--')
    line([allSweeps(1) allSweeps(end)],[1, 1],'Color','red','LineStyle','--')
    axis([allSweeps(1) inf 0 20])
    ylabel('Onset Latency (ms)');
    xlabel('Sweeps');
    title([obj.file ' onset latency ' '(' num2str(inwardORoutward) ')'],'Interpreter','none');
    hold off;
    movegui('southeast');        

% plot light response PEAK latency
    figure('name', strcat(obj.file," (all) - peak latency ", '(',num2str(inwardORoutward),')')); 
    hold on;
    plot(allSweeps, allFirstPeakLatencyInMilliSeconds,'-^');
    plot(allSweeps, allSecondPeakLatencyInMilliSeconds,'-v');
    line([allSweeps(1) allSweeps(end)],[5, 5],'Color','black','LineStyle','--')
    line([allSweeps(1) allSweeps(end)],[1, 1],'Color','red','LineStyle','--')
    axis([allSweeps(1) inf 0 20])
    ylabel('Peak Latency (ms)');
    xlabel('Sweeps');
    title([obj.file ' peak latency ' '(' num2str(inwardORoutward) ')'],'Interpreter','none');
    hold off;
    movegui('south');          
     
end