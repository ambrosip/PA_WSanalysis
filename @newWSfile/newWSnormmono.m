function newWSnormmono(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);

% USED AFter 2020-03-20 (due to long latency of chrimson-evoked light responses)
% Adjusted for Talia's protocol settings
% optargs = {-1 0.02 15 0.29 0.3 0.005 0.24 0.44 -10000 500 0.1 'D:\CORONAVIRUS DATA\From MATLAB'}; 
optargs = {-1 0.02 2 0.99 1 0.001 0.95 1.05 -100 500 0.1 'E:\TEMP'};     % For AJ
optargs(1:numvarargs) = varargin;
[inwardORoutward, analysisTimeWindow, threshold, baselineOnsetTime, lightOnsetTime, lightDuration, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
totalActiveChannels = 2;
sampleRate = obj.header.AcquisitionSampleRate;
lightOnsetDataPoint = lightOnsetTime*sampleRate;
lightOffDataPoint = (lightOnsetTime+lightDuration)*sampleRate;
afterLightDataPoint = (lightOnsetTime+analysisTimeWindow)*sampleRate;
round(afterLightDataPoint);
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
mouseNumber = getMouseNumber(obj);
% experimentDate = getExperimentDate(obj);

% matrixes that will be filled
lightEvokedCurrents = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
seriesResistance = [];
data = [];
allRs = [];
lightEvokedResponseLatency = [];
allLightEvokedResponseLatencyInMilliSeconds = []; 

    % checking for incomplete sweeps and not analyzing incomplete sweeps (to
    % avoid this error: "Dimensions of matrices being concatenated are not
    % consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
    % ych2]" and other concatenations.   
    
    if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
        lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 2;
        allSweeps = firstSweepNumber:lastSweepNumber;
    end 
    
    
    % plotting one figure with all sweeps and average data 
    figure('name', strcat(obj.file,' (all) - 2 channels - Baseline Subtracted')); % naming figure file
    
    % code for channel 2
    for channel = 2:totalActiveChannels
        for sweepNumber = allSweeps
        hold on;
        [xch2,ych2] = obj.xy(sweepNumber, channel);     
        subplot(totalActiveChannels,1,channel)
        plot(xch2,ych2,'k');
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([xmin xmax yminhere ymaxhere])
        ylabel(strcat(obj.header.AIChannelNames(channel), ' (', obj.header.AIChannelUnits(channel), ')'));
        dataPerSweepCh2 = [dataPerSweepCh2, ych2];
        end
        xlabel('Time (s)');
        hold off;
    end
    
    % separate code for channel 1 to apply axis range and figure title
    for sweepNumber = allSweeps  
        hold on;
        [x,y] = obj.xy(sweepNumber, 1);
        
        % calculating series resistance
        rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
        rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
        dCurrent = rsTransientCurrent-rsBaselineCurrent;
        %%%%%%%%%%%%%%%%%%%% IM NOT SURE ABOUT THE dV!!!!!!
        dVoltage = -2.5;
        seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
        
        % normalizing data based on mean baseline current between t=1.99 (or whatever was the input time) and t=2
        y = y-mean(y(baselineOnsetDataPoint:lightOnsetDataPoint));
        
        subplot(totalActiveChannels,1,1)
        plot(x,y,'Color', [0.75, 0.75, 0.75]);
        axis([xmin xmax ymin ymax])
%         xlabel('Time (s)');
        ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([obj.file ' (all) - 2 channels - Baseline Subtracted'],'Interpreter','none');

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
        lightEvokedResponseLatencyInDataPoints = min(find(conv2(sign(diff(y(lightOnsetDataPoint:afterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
        
        % 1000 multiplication is conversion from seconds to milliseconds
        lightEvokedResponseLatencyInMilliSeconds = 1000*lightEvokedResponseLatencyInDataPoints/sampleRate;
        
        if isempty(lightEvokedResponseLatencyInMilliSeconds)
            lightEvokedResponseLatencyInMilliSeconds = NaN;
        end

        lightEvokedCurrent = min(y(lightOnsetDataPoint:afterLightDataPoint));
        
        dataPerSweepCh1 = [dataPerSweepCh1, y];
        lightEvokedCurrents = [lightEvokedCurrents, min(y(lightOnsetDataPoint:afterLightDataPoint))];  
              
        % UPDATE 2022-11-02: I removed experimentDate from data and from
        % the csv code to accomodate AJ's data
        data = [data; mouseNumber, sweepNumber, seriesResistance, lightEvokedCurrent, lightEvokedResponseLatencyInMilliSeconds];
        allRs = [allRs, seriesResistance];
        allLightEvokedResponseLatencyInMilliSeconds = [allLightEvokedResponseLatencyInMilliSeconds, lightEvokedResponseLatencyInMilliSeconds];
        
        
    end

    % plotting average for channel 1
    % dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
    % calculates the average column. mean(dataPerSweepCh1) calculates the
    % wrong thing: the average row.
    subplot(totalActiveChannels,1,1)
    plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
    rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%     line(lightOnsetDataPoint:lightOffDataPoint,100:100,'Color','blue','LineWidth',3)
    axis([xmin xmax ymin ymax])
    hold off;    
    movegui('northwest');
      
    
    % plotting figure with only Ch1 niceplot
        figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted')); % naming figure file
        hold on;
        plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
        plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
        rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
        axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
        xlabel('Time (s)');
        ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
        set(gca,'Visible','off')
        
        xminScale = lightOnsetTime-0.01;
        xmaxScale = lightOnsetTime+0.1;
        
        % adding scale bar
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
        text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.AIChannelUnits(1)))                      
        hold off;
        
    
    % plotting figure with only Ch1 niceplot - AVG only
        figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted - AVG only')); % naming figure file
        hold on;
        plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
        rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
        axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
        xlabel('Time (s)');
        ylabel(strcat("Baseline subtracted ", obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
        set(gca,'Visible','off')
        
        xminScale = lightOnsetTime-0.01;
        xmaxScale = lightOnsetTime+0.1;
        
        % adding scale bar
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
        text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.AIChannelUnits(1)))                      
        hold off;
                   

 % save csv file with data 
        filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - Baseline Subtracted light evoked currents",'(',num2str(inwardORoutward),')');
        fulldirectory = strcat(savefileto,'\',filename,'.csv');        
        dataInCellFormat = {};
        dataInCellFormat = num2cell(data);
        % UPDATE 2022-11-02: I removed experimentDate from data and from
        % the csv code (labeledData) to accomodate AJ's data
        labeledData = cell2table(dataInCellFormat,'VariableNames',...
            {'mouse', 'sweep', 'seriesResistance', 'lightEvokedCurrent', 'responseLatency'});
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
%         axis([allSweeps(1) allSweeps(end) 0 60])
        axis([allSweeps(1) inf 0 60])
        ylabel('Rs (M\Omega)');
        xlabel('Sweeps');
        title([obj.file ' rs'],'Interpreter','none');
        movegui('northeast');
        
 % plot light response latency
        figure('name', strcat(obj.file," (all) - latency ", '(',num2str(inwardORoutward),')')); % naming figure file
        plot(allSweeps, allLightEvokedResponseLatencyInMilliSeconds,'-o');
        line([allSweeps(1) allSweeps(end)],[5, 5],'Color','black','LineStyle','--')
        line([allSweeps(1) allSweeps(end)],[1, 1],'Color','red','LineStyle','--')
%         axis([allSweeps(1) allSweeps(end) 0 10])
%         axis([allSweeps(1) inf 0 10])
        axis([allSweeps(1) inf 0 20])
        ylabel('Response Latency (ms)');
        xlabel('Sweeps');
        title([obj.file ' latency ' '(' num2str(inwardORoutward) ')'],'Interpreter','none');
        movegui('southeast');        
        

end