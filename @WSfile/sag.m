function [dataPerCurrentStep, dataPerSweepCh1, dataPerSweepCh2] = sag(obj,varargin);

    % optional arguments
    % set defaults for optional inputs 
%     optargs = {'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
    optargs = {'D:\Temp\From MATLAB'};
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [savefileto] = optargs{:};

    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    dataPerCurrentStep = [];
    dataPerSweepCh1 = [];
    dataPerSweepCh2 = [];
    allSweeps(1)=[];
    data = [];
    
    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);

    for sweepNumber = allSweeps
        
        [x,y] = obj.xy(sweepNumber, 1);
        [xch2,ych2] = obj.xy(sweepNumber, 2);
        
        sagPeak = min(y);
        steadyState = y(1.45*obj.header.Acquisition.SampleRate); % index 1.45*10000 corresponds to x=1.45 s
        currentStep = round(min(ych2));
        sagRatio = sagPeak/steadyState;    
        
        dataPerCurrentStep = [dataPerCurrentStep; sweepNumber, currentStep, sagRatio, sagPeak, steadyState];
        data = [data; mouseNumber, experimentDate, sweepNumber, currentStep, sagRatio, sagPeak, steadyState];
        dataPerSweepCh1 = [dataPerSweepCh1, x, y];
        dataPerSweepCh2 = [dataPerSweepCh2, xch2, ych2];
        
    end
    
    figure('name', strcat(obj.file, ' - sag ratio per pA'));
    plot(dataPerCurrentStep(:,2),dataPerCurrentStep(:,3),'-o')
    ylabel('Sag Ratio');
    xlabel('Current Step (pA)')
    title([strcat(obj.file, ' - sag ratio per pA')],'Interpreter','none');
    axis([-inf inf 0 2])
    set(gca, 'XDir','reverse') % reverses x axis 
    
    figure('name', strcat(obj.file, ' (-150 pA step)'));
    subplot(2,1,1)
    hold on;
    plot(dataPerSweepCh1(:,5),dataPerSweepCh1(:,6));
    line([0, 3],[dataPerCurrentStep(3,4), dataPerCurrentStep(3,4)],'Color','red','LineStyle','--')
%         before adding "sweeps" to dataPerCurrentStep, this line was like
%         this:
%         line([0, 3],[dataPerCurrentStep(3,3), dataPerCurrentStep(3,3)],'Color','red','LineStyle','--')
    plot(1.45, dataPerSweepCh1(1.45*obj.header.Acquisition.SampleRate,6),'o','color','red');
    axis([-inf, inf, -120, 40])
    ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
%     xlabel('Time (s)')
    title([strcat(obj.file, ' (-150 pA step)')],'Interpreter','none');
    
    subplot(2,1,2)
    plot(dataPerSweepCh2(:,5),dataPerSweepCh2(:,6));
    axis([-inf, inf, -200, 50])
    ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
    xlabel('Time (s)');
    
    % save csv file with dataPerCurrentStep 
        filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - CC sag");
        fulldirectory = strcat(savefileto,'\',filename,'.csv');
%         csvwrite(fulldirectory,dataPerCurrentStep); 
%         csvwrite works but does not let me add strings to label the variables
%         to add labels, I have to use writetable. To use writetable, I
%         have to convert my array into a cell array...
        dataPerCurrentStepInCellFormat = {};
        dataPerCurrentStepInCellFormat = num2cell(data);
        labeledData = cell2table(dataPerCurrentStepInCellFormat,'VariableNames',{'mouse', 'date', 'sweepNumber', 'currentStep', 'sagRatio', 'sagPeakCurrent', 'steadyStateCurrent'})
        writetable(labeledData,fulldirectory);        
        disp('I saved it, ur welcome love')
        disp('Change directory if you want this saved elsewhere!') 

end