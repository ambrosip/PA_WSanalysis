function newExcitability(obj,varargin)

% defaults for optional args
optargs = {10 50 1.25 -10 0.001 'D:\CORONAVIRUS DATA\From MATLAB'};

% overwrite defaults with values specified in varargin
numvarargs = length(varargin);
optargs(1:numvarargs) = varargin;

% place optional args in memorable variable names
[ymaxAP, ymaxISI, midCurrentStepTimePoint, MinPeakHeight, MinPeakDistance, savefileto] = optargs{:};

% get info from file
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
samplingFrequency = obj.header.AcquisitionSampleRate;

% create matrix that will be filled
dataPerCurrentStep = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
dataPerSweepCh3 = [];

% analyze each sweep
for sweepNumber = allSweeps
    
    % collect data from ch1 (voltage) and ch3 (current)
    [x,y] = obj.xy(sweepNumber, 1);         % voltage
    [xch2,ych2] = obj.xy(sweepNumber, 2);   % light
    [xch3,ych3] = obj.xy(sweepNumber, 3);   % current
    
    % find all action potentials in sweep
    [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    
    % find data point index where current step is at its max or min
    midCurrentStepDataPoint = midCurrentStepTimePoint*samplingFrequency;
    
    % round down value of current step
    currentStep = floor(ych3(midCurrentStepDataPoint));
    
    % keep all action potentials (APs) during current step (1 - 1.5 s), ignore the rest 
    locsCurrentStep = locs;
    indicesNotCurrentStep = find(locs<1 | locs>1.5);
    locsCurrentStep(indicesNotCurrentStep) = [];    
    numberAP = length(locsCurrentStep);
    
    % calculate inter-spike interval (ISI) if there are enough APs
    if length(locsCurrentStep) > 1
        allISI = diff(locsCurrentStep);
        firstISIinv = 1/allISI(1);
        lastISIinv = 1/allISI(end);
    else
        allISI = 0;
        firstISIinv = 0;
        lastISIinv = 0;
    end
    
    % collect all data into one matrix
    dataPerCurrentStep = [dataPerCurrentStep; experimentDate, mouseNumber, sweepNumber, currentStep, numberAP, firstISIinv, lastISIinv];
    dataPerSweepCh1 = [dataPerSweepCh1, x, y];
    dataPerSweepCh2 = [dataPerSweepCh2, xch2, ych2];
    dataPerSweepCh3 = [dataPerSweepCh3, xch3, ych3];
    
end

% plot number of APs per current step
figure('name', strcat(obj.file, ' (all) - AP per pA'));
plot(dataPerCurrentStep(:,4),dataPerCurrentStep(:,5),'-o')
ylabel('# APs');
xlabel('Current Step (pA)')
title([strcat(obj.file, ' (all) - AP per pA')],'Interpreter','none');
axis([-inf inf 0 ymaxAP])
movegui('northeast');

% plot first and last ISI per current step
figure('name', strcat(obj.file, ' (all) - first and last ISI per pA'));
hold on;
plot(dataPerCurrentStep(:,4),dataPerCurrentStep(:,6),'-^');
plot(dataPerCurrentStep(:,4),dataPerCurrentStep(:,7),'-v');
hold off;
axis([-inf inf 0 ymaxISI])
ylabel('1/ISI (Hz)');
xlabel('Current Step (pA)')
title([strcat(obj.file, ' (all) - first and last 1/ISI')],'Interpreter','none');
movegui('northwest');

% plot example raw data - voltage vs time
figure('name', strcat(obj.file, ' (subset) - raw'));
subplot(3,1,1)
    hold on;
    plot(dataPerSweepCh1(:,1),dataPerSweepCh1(:,2));
    plot(dataPerSweepCh1(:,end-1),dataPerSweepCh1(:,end));
    axis([-inf, inf, -120, 40])
    hold off;
    ylabel(strcat(obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
    title([strcat(obj.file, ' - subset raw')],'Interpreter','none');
subplot(3,1,2)
    hold on;
    plot(dataPerSweepCh3(:,1),dataPerSweepCh3(:,2));
    plot(dataPerSweepCh3(:,end-1),dataPerSweepCh3(:,end));
    axis([-inf, inf, -150, 250])
    hold off;
    ylabel(strcat(obj.header.AIChannelNames(3), ' (', obj.header.AIChannelUnits(3), ')'));
    xlabel('Time (s)');
subplot(3,1,3)
    hold on;
    plot(dataPerSweepCh2(:,1),dataPerSweepCh2(:,2));
    plot(dataPerSweepCh2(:,end-1),dataPerSweepCh2(:,end));
    axis([-inf, inf, -5, 10])
    hold off;
    ylabel(strcat(obj.header.AIChannelNames(2), ' (', obj.header.AIChannelUnits(2), ')'));
    xlabel('Time (s)');
set(gcf, 'Position', [680,558,560,630]);
movegui('north');

size(dataPerSweepCh1,2)
length(dataPerSweepCh1(:,1))
length(dataPerSweepCh1(:,end))

% % plot all current steps - raw data - voltage vs time
% figure('name', strcat(obj.file, ' (all steps) - raw'));
% for i = 1:length(dataPerCurrentStep(:,4)) 
%     subplot(3,1,1) 
%         hold on;
%         plot(dataPerSweepCh1(:,2*i-1),dataPerSweepCh1(:,2*i));
%         axis([-inf, inf, -120, 40])
%         ylabel(strcat(obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
%         title([strcat(obj.file, ' - all steps')],'Interpreter','none');
%     subplot(3,1,2)
%         hold on;
%         plot(dataPerSweepCh3(:,2*i-1),dataPerSweepCh3(:,2*i));
%         axis([-inf, inf, -150, 250])
%         ylabel(strcat(obj.header.AIChannelNames(3), ' (', obj.header.AIChannelUnits(3), ')'));
%         xlabel('Time (s)');
%     subplot(3,1,3)
%         hold on;
%         plot(dataPerSweepCh2(:,2*i-1),dataPerSweepCh2(:,2*i));
%         axis([-inf, inf, -5, 10])
%         ylabel(strcat(obj.header.AIChannelNames(2), ' (', obj.header.AIChannelUnits(2), ')'));
%         xlabel('Time (s)');
% end
% hold off;
% set(gcf, 'Position', [680,558,560,630]);
% movegui('southwest');



% plot multiple figures - one for every current step - voltage vs time
for i = 1:length(dataPerCurrentStep(:,4)) 
    figure('name', strcat(obj.file, " - raw (", num2str(dataPerCurrentStep(i,4)), ')'));
    subplot(3,1,1) 
        plot(dataPerSweepCh1(:,2*i-1),dataPerSweepCh1(:,2*i));
        axis([-inf, inf, -120, 40])
        ylabel(strcat(obj.header.AIChannelNames(1), ' (', obj.header.AIChannelUnits(1), ')'));
        title([strcat(obj.file, ' - raw (', num2str(dataPerCurrentStep(i,4)), ')')],'Interpreter','none');
    subplot(3,1,2)
        plot(dataPerSweepCh3(:,2*i-1),dataPerSweepCh3(:,2*i));
        axis([-inf, inf, -150, 250])
        ylabel(strcat(obj.header.AIChannelNames(3), ' (', obj.header.AIChannelUnits(3), ')'));
        xlabel('Time (s)');
    subplot(3,1,3)
        plot(dataPerSweepCh2(:,2*i-1),dataPerSweepCh2(:,2*i));
        axis([-inf, inf, -5, 10])
        ylabel(strcat(obj.header.AIChannelNames(2), ' (', obj.header.AIChannelUnits(2), ')'));
        xlabel('Time (s)');
        set(gcf, 'Position', [680,558,560,630]);
        movegui('southwest');
end



% save csv file with dataPerCurrentStep
filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - excitability");
fulldirectory = strcat(savefileto,'\',filename,'.csv');
dataPerCurrentStepInCellFormat = {};
dataPerCurrentStepInCellFormat = num2cell(dataPerCurrentStep);
labeledData = cell2table(dataPerCurrentStepInCellFormat,'VariableNames',{'date', 'mouse', 'sweepNumber', 'currentStep', 'numberAP', 'firstISIinv', 'lastISIinv'});
writetable(labeledData,fulldirectory);        
disp('change directory if you want this saved elsewhere!') 

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

