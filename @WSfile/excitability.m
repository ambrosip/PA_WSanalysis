function [dataPerCurrentStep, dataPerSweepCh1, dataPerSweepCh2] = excitability(obj, varargin);

    % optional arguments
    % set defaults for optional inputs 
    optargs = {-10 0 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [MinPeakHeight, MinPeakDistance, savefileto] = optargs{:};

    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    dataPerCurrentStep = [];
    dataPerSweepCh1 = [];
    dataPerSweepCh2 = [];

    for sweepNumber = allSweeps
        [x,y] = obj.xy(sweepNumber, 1);
        
                       
        if length(dataPerSweepCh1) <= length(x)
            [xch2,ych2] = obj.xy(sweepNumber, 2);

            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
            currentStep = max(ych2);

            locsCurrentStep = locs;
        %     pksCurrentStep = pks;

            indicesNotCurrentStep = find(locs<1 | locs>1.5);
            locsCurrentStep(indicesNotCurrentStep) = [];
        %     pksCurrentStep(indicesNotCurrentStep) = [];

            numberAP = length(locsCurrentStep);

            if length(locsCurrentStep) > 1       
                allISI = diff(locsCurrentStep);
                firstISIinv = 1/allISI(1);
                lastISIinv = 1/allISI(end);
            else
                allISI = 0;
                firstISIinv = 0;
                lastISIinv = 0;
            end         

            dataPerCurrentStep = [dataPerCurrentStep; currentStep, numberAP, firstISIinv, lastISIinv];
            dataPerSweepCh1 = [dataPerSweepCh1, x, y];
            dataPerSweepCh2 = [dataPerSweepCh2, xch2, ych2];
            
        else
            disp('Last sweep was incomplete');
        end

    end
    
    figure('name', strcat(obj.file, ' (all) - AP per pA'));
    plot(dataPerCurrentStep(:,1),dataPerCurrentStep(:,2),'-o')
    ylabel('# APs');
    xlabel('Current Step (pA)')
%     title([strcat(objectArray(1).file(1:15)," - #AP/pA: ",num2str(allSweeps(1)),'-',num2str(allSweeps(end)))],'Interpreter','none');
    title([strcat(obj.file, ' (all) - AP per pA')],'Interpreter','none');
    axis([-inf inf 0 12])

    figure('name', strcat(obj.file, ' (all) - first and last ISI per pA'));
    hold on;
    plot(dataPerCurrentStep(:,1),dataPerCurrentStep(:,3),'-^');
    plot(dataPerCurrentStep(:,1),dataPerCurrentStep(:,4),'-v');
    hold off;
    axis([-inf inf 0 50])
    ylabel('1/ISI (Hz)');
    xlabel('Current Step (pA)')
    title([strcat(obj.file, ' (all) - first and last 1/ISI')],'Interpreter','none');
    
    figure('name', strcat(obj.file, ' (all) - raw'));
    subplot(2,1,1)
    hold on;
    plot(dataPerSweepCh1(:,3),dataPerSweepCh1(:,4));
    plot(dataPerSweepCh1(:,end-1),dataPerSweepCh1(:,end));
    axis([-inf, inf, -85, 40])
    hold off;
    ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
%     xlabel('Time (s)')
    title([strcat(obj.file, ' - subset raw')],'Interpreter','none');
    
    subplot(2,1,2)
    hold on;
    plot(dataPerSweepCh2(:,3),dataPerSweepCh2(:,4));
    plot(dataPerSweepCh2(:,end-1),dataPerSweepCh2(:,end));
    axis([-inf, inf, min(ych2)-50, max(ych2)+50])
    hold off;
    ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
    xlabel('Time (s)');
     
    % save csv file with dataPerCurrentStep 
    filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - excitability");
    fulldirectory = strcat(savefileto,'\',filename,'.csv');
    dataPerCurrentStepInCellFormat = {};
    dataPerCurrentStepInCellFormat = num2cell(dataPerCurrentStep);
    labeledData = cell2table(dataPerCurrentStepInCellFormat,'VariableNames',{'currentStep', 'numberAP', 'firstISIinv', 'lastISIinv'});
    writetable(labeledData,fulldirectory);        
    disp('change directory if you want this saved elsewhere!')    

end




