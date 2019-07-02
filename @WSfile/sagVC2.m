function [data, dataPerSweepCh1, dataPerSweepCh2] = sagVC2(obj,varargin)

    % optional arguments
    % set defaults for optional inputs 
    optargs = {'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [savefileto] = optargs{:};

    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);
    
    data = [];
    dataPerSweepCh1 = [];
    dataPerSweepCh2 = [];

    for sweepNumber = allSweeps
        [x,y] = obj.xy(sweepNumber, 1);
        [xch2,ych2] = obj.xy(sweepNumber, 2);
        
        ysubset = y;
        % VERY RELEVANT: remove elements from the end of the list FIRST,
        % and elements from the beginning of the list LATER, otherwise the
        % ranges set do not match the size of the list!!!
        ysubset(3.5*obj.header.Acquisition.SampleRate:end)=[]; % remove all points after hyperpolarizing step
        ysubset(1:2.0005*obj.header.Acquisition.SampleRate)=[]; % remove all points before hyperpolarizing step        
        leakCurrent = max(ysubset);
        baselineCurrent = y(1.9*obj.header.Acquisition.SampleRate);
%         sagPeak = min(ysubset);
        sagPeak = y(3.5*obj.header.Acquisition.SampleRate); % change 3.5 value when cell dies before ih step protocol is completed
                
        sagRatio = (sagPeak-baselineCurrent)/(leakCurrent-baselineCurrent);    
        sagCurrent = -(sagPeak-leakCurrent);
        
        data = [data; sagRatio, sagCurrent, leakCurrent, baselineCurrent, sagPeak];
        dataPerSweepCh1 = [dataPerSweepCh1, x, y];
        dataPerSweepCh2 = [dataPerSweepCh2, xch2, ych2];
    end    
    
    figure('name', strcat(obj.file, ' - Ih'));
    subplot(2,1,1)
    hold on;
    plot(x,y);
    axis([1.5, 4, -2000, 500])
    line([1.5, 4],[leakCurrent, leakCurrent],'Color','red','LineStyle','--')
    plot(3.5, sagPeak,'o','color','red'); % change 3.5 value when cell dies before ih step protocol is completed
%     text(4.5,250,strcat("Sag Ratio = ", num2str(round(sagRatio,2))),'color','black');
    text(2.1,0,strcat('Sag Current = ', num2str(round(sagCurrent,2)), ' pA'),'color','black');
%     text(LightOnsetTime + LightDur/4,-75,num2str(round(inverseISIforLight,2)),'color','red');
    hold off;
    ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
%     xlabel('Time (s)')
    title([strcat(obj.file, ' - Ih')],'Interpreter','none');    
    
    subplot(2,1,2)
    hold on;
    plot(xch2,ych2);
    axis([1.5, 4, min(ych2)-50, max(ych2)+50])
    ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
    xlabel('Time (s)');
    
    % save csv file with data 
        filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - VC sag");
        fulldirectory = strcat(savefileto,'\',filename,'.csv');
        dataInCellFormat = {};
        dataInCellFormat = num2cell(data);
        labeledData = cell2table(dataInCellFormat,'VariableNames',{'sagRatio', 'sagCurrent', 'leakCurrent', 'baselineCurrent', 'sagPeak'});
        writetable(labeledData,fulldirectory);
        disp('I saved it, ur welcome love')
        disp('Change directory if you want this saved elsewhere!')

end