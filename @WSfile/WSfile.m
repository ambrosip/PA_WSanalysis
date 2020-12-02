% EX:   s184=WSfile('m282_2019-02-12_0184-0188.h5')  load file
%       s184.plotxy(184)                             plot raw sweep 184
%       s184.plotxdy(184)                            plot deriv sweep 184
%       s184.plotxy([184,185])                       plot sweeps 184 & 185 in same figure
%       s184.plotxy([184:187])                       plot sweeps 184-187 in same figure
%       data.ff(109,0,7,10)                          prints firing frequency from sweep 109, threshold 0 mV, from 7-10 s     

%       s96.plotps(99,0,5)                           plots sweep 99, threshold 0, light onset at 5 s

%       arrayfun(@(x) s184.plotxy(x),[184:187])      plot sweeps 184-187 in different figures

% figure;s117.plotpsall(0,15,3);figure;s117.ploth1all(0,10);figure;s117.ploth2all(0,10)

% To hear the data as a sound
% First load the data              file = WSfile(...)
% Then load the sweep data         [x,y] = file.xy(sweepNumber)
% Then you play by using           soundsc(y, sampleRate)


% % % % % firing
% % % % % figure;gates.s0200.plotxy(200, -inf, inf)         1 sweep
% % % % % 
% % % % % Ih
% % % % % figure;gates.s0174.plotxy(174:182,3, 8, -60, 10)	9 swweeps (first + 8)
% % % % % 
% % % % % excitability
% % % % % figure;gates.s0232.plotxy(232:236, 0, 4, -80, 40)	5 sweeps (first + 4)

classdef WSfile
   
    properties
        header   
        sweeps
        file
    end
    
    methods
        plotxyall(self, varargin)
        plotxdyenvall(self, varargin)
        plotxyallch(self, varargin)
        psallch(self, varargin)
        psallON(self, varargin)
        psallONdy(self, varargin)
        plotxyallch2(self, varargin)    % merges 2 sweeps at a time
        plotxyall2ch(self, varargin)    % plots only TWO channels
        plotxyall2ch2(self, varargin)    % plots only TWO channels & merges 2 sweeps
        phaseplot(self)
        [allRs, allSweeps] = rs(self,varargin)    % plots rs
        [seriesResistance, firstSweepNumber, mouseNumber, experimentDate] = calculateRs(obj)
        [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj)
        [dataPerCurrentStep, dataPerSweepCh1, dataPerSweepCh2] = excitability(self, varargin)
        [dataPerCurrentStep, dataPerSweepCh1, dataPerSweepCh2] = sag(self,varargin)
        [data, dataPerSweepCh1, dataPerSweepCh2] = sagVC(self, varargin)
        [lightEvokedCurrents,dataPerSweepCh1,dataPerSweepCh2] = mono(self, varargin)
        [data, dataPerSweepCh1, dataPerSweepCh2] = sagVC2(self,varargin)
        niceplot(self, varargin)
        lightvsfiring(self, varargin)
        lightvsfiringON(self, varargin)
        [lightEvokedCurrents,dataPerSweepCh1,dataPerSweepCh2,seriesResistance] = normmono(self, varargin) % baseline subtracted
        [mouseNumber] = getMouseNumber(self)
        [experimentDate] = getExperimentDate(self)
        normmonofig(self, varargin)
        firingrate(self, varargin)
        firingrateON(self, varargin)
        plotfft(self,sweepNumber,highpassThreshold,lowpassThreshold)
        firingrateONbandpass(self, varargin)
        bandpassplot(self, varargin)
        filteredplot(self, varargin)
        lightvsfiringONbandpass(self, varargin)
        niceplotBandpass(self, varargin)
        lightvsfiringONbandpassAUTO(self, varargin)
        niceplotHighpass(self, varargin)
        niceplotBandpassNew(obj, varargin)
        lightvsfiringONbandpassAUTOhist(obj)
        lightvsfiringONbandpassAUTOhistZOOM(obj)

        
        function obj = WSfile(fileName)
            % Need to put single quotes around file name
            % Ex: WSfile('m282_2019-02-12_0089-0090.h5')
            
            % Load file with sweeps and header using load function from WS
            loadedFile = ws.loadDataFile(fileName);
            
            % Assign header from loaded file to header property
            obj.header = loadedFile.header;
            
            % sweeps is what is left after you remove header
            obj.sweeps = rmfield(loadedFile,'header');
            
            obj.file = fileName;
        end
        
        
        function loadedSweep = sweep(self, sweepNumber)
            % Ex sweepNumber: 90 for sweep_0090  
            sweepName = strcat('sweep_', num2str(sweepNumber, '%04.f')); %maybe change f (float) to integer
            loadedSweep = getfield(self.sweeps, sweepName);
            % MATLAB suggests using new notation: loadedSweep = loadedFile.(sweepName)            
        end

              
        function [x,y] = xy(self, sweepNumber, channel) 
            
            % optional arguments: channel
            if nargin < 3
                selectedChannel = 1;
            else
                selectedChannel = channel;
            end
            
            % Raw y values from sweep
            % Ex sweepNumber: 90 for sweep_0090 
            samplingFrequency = self.header.Acquisition.SampleRate;
            sweepData = self.sweep(sweepNumber);
            
%             % NOTE that you're getting data from Channel 1 - add new
%             % variable if you want to choose channel
%             rawY = sweepData.analogScans(:,1);
            rawY = sweepData.analogScans(:,selectedChannel);
            
%             % Need to scale y values according to Multiclamp gain - only
%             % relevant for early data acquired before reconnecting
%             % communication between amplifier and digidata.
%             y = rawY*self.header.Acquisition.AnalogChannelScales(1)*100;
            y = rawY;
            
            sweepDuration = size(y,1)/samplingFrequency;
            x = linspace(0,sweepDuration,size(y,1))';
        end
        
        
        % varargin: variable argument input (axes range)
        function plotxy(self, sweepNumbers,varargin)
%             figure('Name',[self.file ' (' num2str(sweepNumbers) ')'],'NumberTitle','off');
            hold on;
            for sweepNumber = sweepNumbers
                [x,y] = self.xy(sweepNumber);
                plot(x,y)
            end
            hold off;
            
            % only want 4 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 4
                error('myfuns:somefun2Alt:TooManyInputs', ...
                    'requires at most 4 optional inputs');
            end
            
            % set defaults for optional inputs - auto axes
%             optargs = {-inf inf -80 40};
            optargs = {11 20 -80 40};
            
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            % Place optional args in memorable variable names
            [xmin, xmax, ymin, ymax] = optargs{:};
            axis([xmin xmax ymin ymax])
            
            xlabel('Time (s)');
            ylabel('Voltage (mV) or Current (pA)');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
            
        end
        
        
        function [x,dy] = xdy(self, sweepNumber)
            % Derivative of y values from sweep
            samplingFrequency = self.header.Acquisition.SampleRate;
            sweepData = self.sweep(sweepNumber);
            rawY = sweepData.analogScans(:,1);
%             y = rawY*self.header.Acquisition.AnalogChannelScales(1)*100;
%             scaling no longer applies!
            y = rawY;
            dy = diff(y);
            
            duration = size(y,1)/samplingFrequency;
            x = linspace(0,duration,size(y,1))';
            x(end) = [];
        end
        
        
        function plotxdy(self, sweepNumber)
            figure
            [x,y] = self.xdy(sweepNumber);
            plot(x,y)
        end
        
        
        function [pks,locs,w,p] = peaks(self, sweepNumber, MinPeakHeight)
            % MinPeakHeight is the threshold value above which function
            % will find peaks.
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',0.005); %% added 'MinPeakDistance',0.005
            % pks   peak y value
            % locs  peak x value
            % w     peak half-width
            % p     peak amplitude
        end
        
        
        function plotp(self,sweepNumber,MinPeakHeight)
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',0.005); %% added 'MinPeakDistance',0.005
            
            self.plotxy(sweepNumber);
            hold on;
            plot(locs,pks,'o');
            hold off;
        end
        
        
        % plots selected peaks: 3 APs prior to light pulse, all APs during
        % light pulse, and 3 APs after light pulse. Displays firing
        % frequency as 1/ISI (ISI: Inter Spike Interval).
        function plotps(self,sweepNumber,varargin)
                       
            
            % set defaults for optional inputs
%             original inputs
%             optargs = {0 15 1 5 1.3 1.1}; 
%           inputs being used for data summary on 2019-04-01
            optargs = {0 15 1 10 1 1};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [MinPeakHeight, LightOnsetTime, LightDur, NumPeaksBeforeAfter, DelayScaleFactor, LightExtensionFactor] = optargs{:};            
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',0.005); %% added 'MinPeakDistance',0.005
            
            % Save locs and pks values in different variable names that 
            % will be editted differently
            locsLight = locs;
            locsBaseline1 = locs;
            locsBaseline2 = locs;
            
            pksLight = pks;
            pksBaseline1 = pks;
            pksBaseline2 = pks;
            
            % indices are integers. Need to find indices at which the value
            % of locs is between LightOnsetTime and LightOnsetTime+1.25. Remember that locs stores the
            % time stamp at which a peak was detected
            indicesLight = find(locs<LightOnsetTime | locs>(LightOnsetTime+LightDur*LightExtensionFactor));
            % deletes all peaks that happen before or after the light pulse
            % UPDATE: get data beyond light pulse (0.25 s after light pulse
            % ended)
            locsLight(indicesLight) = [];
            pksLight(indicesLight) = [];
            
            indicesBaseline1 = find(locs>=LightOnsetTime);
            locsBaseline1(indicesBaseline1) = [];
            pksBaseline1(indicesBaseline1) = [];
            locsBaseline1(1:length(locsBaseline1)-NumPeaksBeforeAfter)=[];
            pksBaseline1(1:length(pksBaseline1)-NumPeaksBeforeAfter)=[];
            
            indicesBaseline2 = find(locs<=(LightOnsetTime+LightDur*DelayScaleFactor));
            locsBaseline2(indicesBaseline2) = [];
            pksBaseline2(indicesBaseline2) = [];
            locsBaseline2(NumPeaksBeforeAfter+1:length(locsBaseline2))=[];
            pksBaseline2(NumPeaksBeforeAfter+1:length(pksBaseline2))=[];
            
            inverseISIforLight = 1/mean(diff(locsLight));
            inverseISIforBaseline1 = 1/mean(diff(locsBaseline1));
            inverseISIforBaseline2 = 1/mean(diff(locsBaseline2));
                        
            self.plotxy(sweepNumber);
            hold on;
            plot(locsLight,pksLight,'o','color','red');
            plot(locsBaseline1,pksBaseline1,'o','color','blue');
            plot(locsBaseline2,pksBaseline2,'o','color','blue');
            rectangle('Position',[LightOnsetTime,-85,LightDur,145],'FaceColor', [0 0 1 0.05],'LineStyle','none');
            hold off;
            axis([LightOnsetTime-4 LightOnsetTime+LightDur+5 -85 45]);
            xlabel('Time (s)');
            ylabel('Voltage (mV)');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
            text(LightOnsetTime-3,-75,'1/ISI (Hz)');
            text(LightOnsetTime+0.25,-75,num2str(round(inverseISIforLight,2)),'color','red');
            text(LightOnsetTime-1,-75,num2str(round(inverseISIforBaseline1,2)),'color','blue');
            text(LightOnsetTime+LightDur+3,-75,num2str(round(inverseISIforBaseline2,2)),'color','blue');
        end
        
        
%         function plotpsall(self,MinPeakHeight,LightOnsetTime)
%             % This converts struct into array so I can use indexing to
%             % access each sweep. sweepNamesFromFileStruct(1) gives me data
%             % from first sweep. The data stored is the same as the data
%             % stored in loadedSweep from function sweep(self, sweepNumber).
%             % sweepNamesFromFileStruct(1:end) lists the data of all sweeps
%             % from the file.
%             sweepNamesFromFileStruct = struct2array(self.sweeps);
%             numberOfSweepsInFile = size(sweepNamesFromFileStruct,2);
%             arrayfun(@(x) self.xy(x),[1:end])
%             
%         end
        

        function plotpsall(self,MinPeakHeight,LightOnsetTime,varargin)
            
            % set defaults for optional inputs 
            optargs = {1 10};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [LightDur, NumPeaksBeforeAfter] = optargs{:};
            
            firstSweepNumber = str2num(self.file(end-11:end-8));
            lastSweepNumber = str2num(self.file(end-6:end-3));
            plotIndex = 1;            
            for sweepNumber = firstSweepNumber:lastSweepNumber
                subplot(1,lastSweepNumber-firstSweepNumber+1,plotIndex);
                hold on;
                self.plotps(sweepNumber,MinPeakHeight,LightOnsetTime,LightDur,NumPeaksBeforeAfter);
                plotIndex=plotIndex+1;
                set(gcf,'Position',[100 200 1750 375]) % for 3 plots
%                 set(gcf,'Position',[100 200 2200 375])  % for 5 plots
            end
            hold off;
        end    
            
            
%             while firstSweepNumber < firstSweepNumber + size(struct2array(self.sweeps),2)-1
% %                 subplotColumnIndex=1;
% %                 subplot(1,size(struct2array(self.sweeps),2),subplotColumnIndex);
%                 self.plotps(firstSweepNumber,MinPeakHeight,LightOnsetTime);
%                 firstSweepNumber = firstSweepNumber+1;
% %                 subplotColumnIndex = subplotColumnIndex+1;
%             end
         
        
        function firingFrequency = ff(self, sweepNumber, MinPeakHeight, timeStart, timeEnd)
            [pks,locs,w,p] = self.peaks(sweepNumber, MinPeakHeight);
            
            % l = locs(locs>=timeStart && locs<= timeEnd) ;
            % numPeaks = length(l);
            
            numPeaks = sum((locs>=timeStart & locs<=timeEnd));
            firingFrequency = numPeaks/(timeEnd-timeStart);          
        end
        
        
        % plots histogram of firing frequency based on 1/ISI
        function ploth1(self,sweepNumber,MinPeakHeight,varargin)

            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight);

            sweepDuration = self.header.Acquisition.Duration;
            sweepTime=0;
            inverseISIperSecBin=[];
            
            while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
                indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
                locsDuringSweepTime = locs;
                locsDuringSweepTime(indicesToDelete) = [];
                inverseISIperSecBin=[inverseISIperSecBin 1/mean(diff(locsDuringSweepTime))];
%                 inverseISIperSecBin=[inverseISIperSecBin; sweepTime 1/mean(diff(locsDuringSweepTime))];
%                 locsDuringSweepTime
%                 diff(locsDuringSweepTime)
%                 1/mean(diff(locsDuringSweepTime))
%                 inverseISIperSecBin
                sweepTime=sweepTime+1;
            end
            
%             histogram(inverseISIperSecBin);
            bar(1:30,inverseISIperSecBin,1);
            axis([-inf inf 0 ymax]);
            
            xlabel('Bins (1 s long)');
            ylabel('Firing Frequency (Mean 1/ISI)');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        end
        
        
        % plots histogram of firing frequency based number of action
        % potentials per second
        function ploth2(self,sweepNumber,MinPeakHeight,varargin)

            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight);

            histogram(locs,0:1:30)
            axis([-inf inf 0 ymax]);
            
            xlabel('Bins (1 s long)');
            ylabel('Number of Action Potentials');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        end
        
        function ploth1all(self,MinPeakHeight,varargin)
            
            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            firstSweepNumber = str2num(self.file(end-11:end-8));
            lastSweepNumber = str2num(self.file(end-6:end-3));
            plotIndex = 1;            
            for sweepNumber = firstSweepNumber:lastSweepNumber
                subplot(1,lastSweepNumber-firstSweepNumber+1,plotIndex);
                hold on;
                self.ploth1(sweepNumber,MinPeakHeight,ymax);
                plotIndex=plotIndex+1;
                set(gcf,'Position',[100 200 1750 375]) % for 3 plots
%                 set(gcf,'Position',[100 200 2200 375])  % for 5 plots
            end
            hold off;
        end
        
        function ploth2all(self,MinPeakHeight,varargin)
            
            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            firstSweepNumber = str2num(self.file(end-11:end-8));
            lastSweepNumber = str2num(self.file(end-6:end-3));
            plotIndex = 1;            
            for sweepNumber = firstSweepNumber:lastSweepNumber
                subplot(1,lastSweepNumber-firstSweepNumber+1,plotIndex);
                hold on;
                self.ploth2(sweepNumber,MinPeakHeight,ymax);
                plotIndex=plotIndex+1;
                set(gcf,'Position',[100 200 1750 375]) % for 3 plots
%                 set(gcf,'Position',[100 200 2200 375])  % for 5 plots
            end
            hold off;
        end
        
        
        % finds minima (findpeaks applied to inverse of data) in cell
        % attached recordings (ON). Minima are more informative of action
        % potentials than maxima (peaks).
        % Displays firing frequency as 1/ISI (ISI: Inter Spike Interval).
        function plotpsON(self,sweepNumber,varargin)
            
            % set defaults for optional inputs
            optargs = {0 15 1 10 1 1};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [MinPeakHeight, LightOnsetTime, LightDur, NumPeaksBeforeAfter, DelayScaleFactor, LightExtensionFactor] = optargs{:};            
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(-y,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',0.005); 
            
            % Save locs and pks values in different variable names that 
            % will be editted differently
            locsLight = locs;
            locsBaseline1 = locs;
            locsBaseline2 = locs;
            
            pksLight = pks;
            pksBaseline1 = pks;
            pksBaseline2 = pks;
            
            % indices are integers. Need to find indices at which the value
            % of locs is between LightOnsetTime and LightOnsetTime+1.25. Remember that locs stores the
            % time stamp at which a peak was detected
            indicesLight = find(locs<LightOnsetTime | locs>(LightOnsetTime+LightDur*LightExtensionFactor));
            % deletes all peaks that happen before or after the light pulse
            locsLight(indicesLight) = [];
            pksLight(indicesLight) = [];
            
            indicesBaseline1 = find(locs>=LightOnsetTime);
            locsBaseline1(indicesBaseline1) = [];
            pksBaseline1(indicesBaseline1) = [];
            locsBaseline1(1:length(locsBaseline1)-NumPeaksBeforeAfter)=[];
            pksBaseline1(1:length(pksBaseline1)-NumPeaksBeforeAfter)=[];
            
            indicesBaseline2 = find(locs<=(LightOnsetTime+LightDur*DelayScaleFactor));
            locsBaseline2(indicesBaseline2) = [];
            pksBaseline2(indicesBaseline2) = [];
            locsBaseline2(NumPeaksBeforeAfter+1:length(locsBaseline2))=[];
            pksBaseline2(NumPeaksBeforeAfter+1:length(pksBaseline2))=[];
            
            inverseISIforLight = 1/mean(diff(locsLight));
            inverseISIforBaseline1 = 1/mean(diff(locsBaseline1));
            inverseISIforBaseline2 = 1/mean(diff(locsBaseline2));
                        
            self.plotxy(sweepNumber);
            hold on;
            plot(locsLight,-pksLight,'o','color','red');
            plot(locsBaseline1,-pksBaseline1,'o','color','blue');
            plot(locsBaseline2,-pksBaseline2,'o','color','blue');
            rectangle('Position',[LightOnsetTime,-85,LightDur,145],'FaceColor', [0 0 1 0.05],'LineStyle','none');
            hold off;
            axis([LightOnsetTime-4 LightOnsetTime+LightDur+5 -85 45]);
            xlabel('Time (s)');
            ylabel('Voltage (mV)');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
            text(LightOnsetTime-3,-75,'1/ISI (Hz)');
            text(LightOnsetTime+0.25,-75,num2str(round(inverseISIforLight,2)),'color','red');
            text(LightOnsetTime-1,-75,num2str(round(inverseISIforBaseline1,2)),'color','blue');
            text(LightOnsetTime+LightDur+3,-75,num2str(round(inverseISIforBaseline2,2)),'color','blue');
        end
        
        
        function plotpsONall(self,MinPeakHeight,LightOnsetTime,varargin)
            
            % set defaults for optional inputs 
            optargs = {1 10};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [LightDur, NumPeaksBeforeAfter] = optargs{:};
            
            firstSweepNumber = str2num(self.file(end-11:end-8));
            lastSweepNumber = str2num(self.file(end-6:end-3));
            plotIndex = 1;            
            for sweepNumber = firstSweepNumber:lastSweepNumber
                subplot(1,lastSweepNumber-firstSweepNumber+1,plotIndex);
                hold on;
                self.plotpsON(sweepNumber,MinPeakHeight,LightOnsetTime,LightDur,NumPeaksBeforeAfter);
                plotIndex=plotIndex+1;
                set(gcf,'Position',[100 200 1750 375]) % for 3 plots
%                 set(gcf,'Position',[100 200 2200 375])  % for 5 plots
            end
            hold off;
        end
        
        
        % plots histogram of firing frequency based on 1/ISI
        function ploth1ON(self,sweepNumber,MinPeakHeight,varargin)

            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(-y,x,'MinPeakHeight',MinPeakHeight);

            sweepDuration = self.header.Acquisition.Duration;
            sweepTime=0;
            inverseISIperSecBin=[];
            
            while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
                indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
                locsDuringSweepTime = locs;
                locsDuringSweepTime(indicesToDelete) = [];
                inverseISIperSecBin=[inverseISIperSecBin 1/mean(diff(locsDuringSweepTime))];

                sweepTime=sweepTime+1;
            end
            
            bar(1:30,inverseISIperSecBin,1);
            axis([-inf inf 0 ymax]);
            
            xlabel('Bins (1 s long)');
            ylabel('Firing Frequency (Mean 1/ISI)');
            title([self.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        end
        
        function ploth1ONall(self,MinPeakHeight,varargin)
            
            % set defaults for optional inputs 
            optargs = {6};
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            numvarargs = length(varargin);
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [ymax] = optargs{:};
            
            firstSweepNumber = str2num(self.file(end-11:end-8));
            lastSweepNumber = str2num(self.file(end-6:end-3));
            plotIndex = 1;            
            for sweepNumber = firstSweepNumber:lastSweepNumber
                subplot(1,lastSweepNumber-firstSweepNumber+1,plotIndex);
                hold on;
                self.ploth1ON(sweepNumber,MinPeakHeight,ymax);
                plotIndex=plotIndex+1;
                set(gcf,'Position',[100 200 1750 375]) % for 3 plots
%                 set(gcf,'Position',[100 200 2200 375])  % for 5 plots
            end
            hold off;
        end
        
        
    end
end