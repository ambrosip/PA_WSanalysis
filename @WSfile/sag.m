function [dataPerCurrentStep, dataPerSweepCh1, dataPerSweepCh2] = sag(obj);

    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    dataPerCurrentStep = [];
    dataPerSweepCh1 = [];
    dataPerSweepCh2 = [];
    allSweeps(1)=[];

    for sweepNumber = allSweeps
        
        [x,y] = obj.xy(sweepNumber, 1);
        [xch2,ych2] = obj.xy(sweepNumber, 2);
        
        sagPeak = min(y);
        steadyState = y(1.45*obj.header.Acquisition.SampleRate); % index 1.45*10000 corresponds to x=1.45 s
        currentStep = min(ych2);
        sagRatio = sagPeak/steadyState;    
        
        dataPerCurrentStep = [dataPerCurrentStep; currentStep, sagRatio, sagPeak, steadyState];
        dataPerSweepCh1 = [dataPerSweepCh1, x, y];
        dataPerSweepCh2 = [dataPerSweepCh2, xch2, ych2];

    end
    
    figure('name', strcat(obj.file, ' - sag ratio per pA'));
    plot(dataPerCurrentStep(:,1),dataPerCurrentStep(:,2),'-o')
    ylabel('Sag Ratio');
    xlabel('Current Step (pA)')
    title([strcat(obj.file, ' - sag ratio per pA')],'Interpreter','none');
    axis([-inf inf 0 5])
    set(gca, 'XDir','reverse') % reverses x axis 
    
    figure('name', strcat(obj.file, ' (-150 pA step)'));
    subplot(2,1,1)
    hold on;
    plot(dataPerSweepCh1(:,5),dataPerSweepCh1(:,6));
    line([0, 3],[dataPerCurrentStep(3,3), dataPerCurrentStep(3,3)],'Color','red','LineStyle','--')
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
    

end