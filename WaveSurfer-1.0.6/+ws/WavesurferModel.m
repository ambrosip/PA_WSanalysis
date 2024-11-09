classdef WavesurferModel < ws.Model & ws.EventBroadcaster
    % The main Wavesurfer model object.
    
    properties (Constant = true, Transient=true)
        NFastProtocols = 6
    end

    properties (Dependent = true)
        AllDeviceNames
        PrimaryDeviceName
        IsPrimaryDeviceAPXIDevice
%         AvailableReferenceClockSources
%         ReferenceClockSource
%         ReferenceClockRate
        NDIOTerminalsPerDevice
        NPFITerminalsPerDevice
        NCountersPerDevice
        NAITerminalsPerDevice
        AITerminalIDsOnEachDevice
        NAOTerminalsPerDevice
        AllChannelNames
        IsAIChannelTerminalOvercommitted
        IsAOChannelTerminalOvercommitted
        IsDIChannelTerminalOvercommitted
        IsDOChannelTerminalOvercommitted
        IsProcessingIncomingCommand
        IsAIChannelActive
        AIChannelScales
        AIChannelUnits
        IsDIChannelActive
        IsAIChannelMarkedForDeletion
        IsDIChannelMarkedForDeletion
        AcquisitionSampleRate  % Hz
        ExpectedSweepScanCount
        IsXSpanSlavedToAcquistionDuration
          % if true, the x span for all the scopes is set to the acquisiton
          % sweep duration
        IsXSpanSlavedToAcquistionDurationSettable
          % true iff IsXSpanSlavedToAcquistionDuration is currently
          % settable
        XSpan  % s, the span of time showed in the signal display
        NAIChannels
        NDIChannels
        AOChannelScales
          % A row vector of scale factors to convert each channel from native units to volts on the coax.
          % This is implicitly in units of ChannelUnits per volt (see below)        
        AOChannelUnits
          % An cellstring row vector that describes the real-world units 
          % for each analog stimulus channel.          
        NAOChannels
        NDOChannels
        IsDOChannelTimed
        DOChannelStateIfUntimed
        StimulationSampleRate  % Hz
        IsAOChannelMarkedForDeletion
        IsDOChannelMarkedForDeletion
        %IsTestPulsing
        DoSubtractBaselineInTestPulseView
        TestPulseYLimits
        TestPulseDuration
        IsAutoYInTestPulseView
        IsAutoYRepeatingInTestPulseView
        DoTrodeUpdateBeforeRun
        %IsElectrodeMarkedForTestPulse  % logical, nElectrodes x 1
        IsElectrodeMarkedForRemoval  % logical, nElectrodes x 1
        TestPulseElectrodeIndex  % either a scalar or [] if no electrode is currently selected in the test pulser
        NextRunAbsoluteFileName
        UserClassName
        IsUserClassNameValid
        DoesTheUserObjectMatchTheUserClassName
        TheUserObject
        ElectrodeCount
        IsInControlOfSoftpanelModeAndGains        
        DidLastElectrodeUpdateWork
        AreSoftpanelsEnabled
        IsDoTrodeUpdateBeforeRunSensible
        %TestPulseElectrodeNames
        TestPulseElectrodesCount
        %TestPulseElectrodeAmplitude
        %TestPulseElectrodeName
        IsStimulationEnabled
        IsLoggingEnabled
        IsDisplayEnabled
        DataFileLocation
        DataFileBaseName
        IsOKToOverwriteDataFile
        NPlots
        NextSweepIndex
        DisplayUpdateRate
        AIScalingCoefficients
        IsAIChannelDisplayed
        IsDIChannelDisplayed
        AreColorsNormal
        XOffset
        YLimitsPerAIChannel
        IsGridOn
        DoShowZoomButtons
        DoColorTraces
        AreYLimitsLockedTightToDataForAIChannel
        ChannelIndexWithinTypeFromPlotIndex
        IsAnalogFromPlotIndex
        ChannelIndexFromPlotIndex
        CacheInputChannelIndexFromInputChannelIndex
        PlotHeightFromPlotIndex
        PlotIndexFromChannelIndex  % 1 x nChannels
        AIChannelTerminalNames
        DIChannelTerminalNames
        AOChannelTerminalNames
        DOChannelTerminalNames
        DoRepeatStimulusSequence
        DoIncludeDateInDataFileName
        DoIncludeSessionIndexInDataFileName
        SessionIndex
        %StimulusLibrary
        CurrentRunAbsoluteFileName
        PlotHeightFromAIChannelIndex
        PlotHeightFromDIChannelIndex
        RowIndexFromAIChannelIndex
        RowIndexFromDIChannelIndex
        NRunsCompleted  % number of runs *completed* (not stopped or aborted) since WS was started
        AIChannelDeviceNames
        AOChannelDeviceNames
        DIChannelDeviceNames
        DOChannelDeviceNames
        AIChannelTerminalIDs
        AOChannelTerminalIDs
        DIChannelTerminalIDs
        DOChannelTerminalIDs
        IsPerformingRun
        IsPerformingSweep
        DataCacheDurationWhenContinuous
        IsInputChannelInCacheFromInputChannelIndex
        WidthOfPlotsInPixels
        XDataForDisplay
        YDataForDisplay
    end
    
    properties (Access=protected)
        PrimaryDeviceName_ = ''   % an empty string represents "no device specified"
        %ReferenceClockSource_ = '100MHzTimebase' 
    end

    properties (Access=protected, Transient=true)
        AllDeviceNames_ = cell(1,0)   % transient b/c we want to probe the hardware on startup each time to get this        
        IsPrimaryDeviceAPXIDevice_ = false  % transient b/c depends on hardware
        %AvailableReferenceClockSources_ = cell(1,0)
        %OnboardClockReferenceClockRate_ = []  % Hz
        
        % The terminal counts are transient b/c e.g. "Dev1" could refer to a different board on protocol
        % file load than it did when the protocol file was saved.  Further, the set
        % of available device could be completely different.
        NDIOTerminalsPerDevice_ = zeros(1,0)  
        NPFITerminalsPerDevice_ = zeros(1,0)
        NCountersPerDevice_ = zeros(1,0)
        NAITerminalsPerDevice_ = zeros(1,0)
        AITerminalIDsOnEachDevice_ = cell(1,0)
        NAOTerminalsPerDevice_ = zeros(1,0)

        IsAIChannelTerminalOvercommitted_ = false(1,0)        
        IsAOChannelTerminalOvercommitted_ = false(1,0)        
        
        IsDIChannelTerminalOvercommitted_ = false(1,0)        
        IsDOChannelTerminalOvercommitted_ = false(1,0)        
        
        NRunsCompleted_ = 0  % number of runs *completed* (not stopped or aborted) since WS was started
        DoUsePreferences_ = true
        TheBigTimer_ = []
        AllowTimerCallback_ = true ;
        CurrentProfileName_
        ProfileNames_
        LastProtocolFilePath_ = ''
    end   

    properties (Dependent = true)
        HasUserSpecifiedProtocolFileName
        AbsoluteProtocolFileName
        HasUserSpecifiedUserSettingsFileName
        AbsoluteUserSettingsFileName
        %FastProtocols
        IndexOfSelectedFastProtocol   % Invariant: Always a scalar real double, and an integer between 1 and NFastProtocols (never empty)
        %Acquisition
        %Stimulation
        %Triggering
        %Display
        %Logging
        %UserCodeManager
        %Ephys
        SweepDuration  % the sweep duration, in s
        AreSweepsFiniteDuration  % boolean scalar, whether the current acquisition mode is sweep-based.
        AreSweepsContinuous  
          % boolean scalar, whether the current acquisition mode is continuous.  Invariant: self.AreSweepsContinuous == ~self.AreSweepsFiniteDuration
        NSweepsPerRun
        SweepDurationIfFinite
        NSweepsCompletedInThisRun    % Current number of completed sweeps while the run is running (range of 0 to NSweepsPerRun).
        IsYokedToScanImage
        NTimesDataAvailableCalledSinceRunStart
        ClockAtRunStart  
          % We want this written to the data file header, but not persisted in
          % the .cfg file.  Having this property publically-gettable, and having
          % ClockAtRunStart_ transient, achieves this.
        State
        VersionString
          % VersionString property exists so that the version is written to the data file
          % header.  The version string that gets written to the protocol file is stored
          % outside the serialization of the WavesurferModel, and is gotten directly from
          % ws.versionString().          
        IsITheOneTrueWavesurferModel  % deprecated, same as IsAwake
        IsAwake
        %WarningLog
        %LayoutForAllWindows
       
        AIChannelNames
        DIChannelNames
        AOChannelNames
        DOChannelNames
        
        % Triggering settings
        TriggerCount
        CounterTriggerCount
        ExternalTriggerCount
        AcquisitionTriggerIndex  % this is an index into Schemes
        StimulationUsesAcquisitionTrigger  % boolean
        StimulationTriggerIndex  % this is an index into Schemes
        DoUsePreferences
        DoesProtocolNeedSave
        
        IsWavesurferMainFigureVisible
        IsGeneralSettingsFigureVisible
        IsChannelsFigureVisible
        IsStimulusLibraryFigureVisible
        IsStimulusPreviewFigureVisible
        IsTriggersFigureVisible
        IsUserCodeManagerFigureVisible
        IsElectrodeManagerFigureVisible
        IsTestPulserFigureVisible        
        IsFastProtocolsFigureVisible

        MainFigurePosition
        GeneralSettingsFigurePosition
        ChannelsFigurePosition
        StimulusLibraryFigurePosition
        StimulusPreviewFigurePosition
        TriggersFigurePosition
        UserCodeManagerFigurePosition
        ElectrodeManagerFigurePosition
        TestPulserFigurePosition        
        FastProtocolsFigurePosition        
        
        CurrentProfileName
        ProfileNames
        LastProtocolFilePath
    end
   
%     properties (Access=protected, Constant = true, Transient=true)
%         CommunicationFolderName_ = tempdir()
%         
%         WSCommandFileName_ = 'ws_command.txt'  % Commands *to* WS
%         WSResponseFileName_ = 'ws_response.txt'  % Responses *from* WS
%         
%         SICommandFileName_ = 'si_command.txt'  % Commands *to* SI
%         SIResponseFileName_ = 'si_response.txt'  % Responses *from* SI
%     end
%     
%     properties (Access=protected, Dependent = true)
%         WSCommandFilePath_
%         WSResponseFilePath_
%         SICommandFilePath_
%         SIResponseFilePath_
%     end
    
    properties (Access=protected)
        % Saved to protocol file
        Triggering_
        Acquisition_
        Stimulation_
        Display_
        Ephys_
        UserCodeManager_
        %IsYokedToScanImage_ = false
        %ExecutingScanImageCommandNow_ = false;
        AreSweepsFiniteDuration_ = true
        NSweepsPerRun_ = 1
        SweepDurationIfFinite_ = 1  % s
        
        %LayoutForAllWindows_  % Yeah, this is view-related, but it's persisted, so it belongs in the model
        
        IsGeneralSettingsFigureVisible_ = false
        IsChannelsFigureVisible_ = false
        IsStimulusLibraryFigureVisible_ = false
        IsStimulusPreviewFigureVisible_ = false
        IsTriggersFigureVisible_ = false
        IsUserCodeManagerFigureVisible_ = false
        IsElectrodeManagerFigureVisible_ = false
        IsTestPulserFigureVisible_ = false
        
        MainFigurePosition_ = []
        GeneralSettingsFigurePosition_ = []
        ChannelsFigurePosition_ = []
        StimulusLibraryFigurePosition_ = []
        StimulusPreviewFigurePosition_ = []
        TriggersFigurePosition_ = []
        UserCodeManagerFigurePosition_ = []
        ElectrodeManagerFigurePosition_ = []
        TestPulserFigurePosition_ = []
        
        % Saved to .usr file
        FastProtocols_ = cell(1,0)
        
        % Not saved to either protocol or .usr file
        Logging_
    end

    properties (Access=protected, Transient=true)
        Looper_
        Refiller_        
%         IPCPublisher_
%         LooperIPCSubscriber_
%         RefillerIPCSubscriber_
%         LooperIPCRequester_
%         RefillerIPCRequester_
        HasUserSpecifiedProtocolFileName_ = false
        AbsoluteProtocolFileName_ = ''
        HasUserSpecifiedUserSettingsFileName_ = false
        DoesProtocolNeedSave_ = false 
        AbsoluteUserSettingsFileName_ = ''
        IndexOfSelectedFastProtocol_ = 1  % Invariant: Always a scalar real double, and an integer between 1 and NFastProtocols (never empty)
        State_ = 'uninitialized'
        %Subsystems_
        t_ = 0  % During a sweep, the time stamp of the scan *just after* the most recent scan
        %NScansAcquiredSoFarThisSweep_
        FromRunStartTicId_
        FromSweepStartTicId_
        TimeOfLastWillPerformSweep_
        %TimeOfLastSamplesAcquired_
        NTimesDataAvailableCalledSinceRunStart_ = 0
        %PollingTimer_
        %MinimumPollingDt_
        TimeOfLastPollInSweep_
        ClockAtRunStart_
        %DoContinuePolling_
        %DidLooperCompleteSweep_
        %DidRefillerCompleteEpisodes_
        DidAnySweepFailToCompleteSoFar_
        %WasRunStopped_
        %WasRunStoppedInLooper_
        %WasRunStoppedInRefiller_
        %WasExceptionThrown_
        %ThrownException_
        NSweepsCompletedInThisRun_ = 0
        AreAllSweepsCompleted_
        IsAwake_
        IsHeaded_
        DesiredNScansPerUpdate_        
        NScansPerUpdate_        
        SamplesBuffer_
        IsPerformingRun_ = false
        IsPerformingSweep_ = false
        %IsDeeplyIntoPerformingSweep_ = false
        %TimeInSweep_  % wall clock time since the start of the sweep, updated each time scans are acquired        
        %DoLogWarnings_ = false
        UnmatchedLogWarningStartCount_ = 0  % technically, the number of starts without a corresponding stop
        WarningCount_ = 0
        WarningLog_ = MException.empty(0,1)   % N.B.: a col vector
        %LayoutForAllWindows_ = []   % this should eventually get migrated into the persistent state, but don't want to deal with that now
        DrawnowTicId_
        TimeOfLastDrawnow_
        %DidLooperCompleteSweep_
        %SICommandPollTimer_
        %SICommandFileExistenceChecker_
        IsFastProtocolsFigureVisible_ = false
        FastProtocolsFigurePosition_ = []
        CommandClient_
        CommandServer_
    end
    
    events
        Update
        UpdateMain
        UpdateGeneral
        UpdateChannels
        UpdateTriggering
        UpdateStimulusLibrary
        UpdateStimulusPreview
        UpdateFastProtocols
        UpdateLogging
        UpdateElectrodeManager
        UpdateUserCodeManager
        UpdateForNewData
        UpdateIsYokedToScanImage
        %WillSetState
        DidSetState
        DidCompleteSweep
        UpdateDigitalOutputStateIfUntimed
        DidChangeNumberOfInputChannels
        RequestLayoutForAllWindows
        LayoutAllWindows
        DidSetAcquisitionSampleRate
        DidSetStimulationSampleRate
        UpdateElectrodes
        UpdateTestPulser
        RaiseDialogOnException
        DidMaybeChangeProtocol
        DidMaybeSetUserClassName
        DidSetSingleFigureVisibility
        UpdateDoIncludeSessionIndexInDataFileName
        TPDidSetIsInputChannelActive
        TPUpdateTrace
        
        EMDidSetIsInputChannelActive
        EMDidSetIsDigitalOutputTimed
        EMDidChangeNumberOfInputChannels
        EMDidChangeNumberOfOutputChannels
        
        %UpdateDisplay
        DidSetUpdateRate
        DidSetXSpan
        DidSetXOffset
        DidSetYAxisLimits
        UpdateTraces
        UpdateAfterDataAdded
        
        UpdateReadiness        
    end
    
    properties (Dependent = true, SetAccess=immutable, Transient=true)
        IsReady  % true <=> figures are showing the normal (as opposed to waiting) cursor
    end
    
    properties (Access = protected, Transient=true)
        DegreeOfReadiness_ = 1
    end

    methods
        function self = WavesurferModel(isAwake, isHeaded, doUsePreferences)
            %self@ws.Model();
            
            if ~exist('isAwake','var') || isempty(isAwake) ,
                isAwake = false ;
            end                       
            if isAwake ,
                if ~exist('isHeaded','var') || isempty(isHeaded) ,
                    isHeaded = false ;
                end
            else
                isHeaded = false ;                
            end         
            if isAwake ,
                if ~exist('doUsePreferences','var') || isempty(doUsePreferences) ,
                    doUsePreferences = true ;
                end
            else
                doUsePreferences = false ;                
            end                     
            
            %doRunInDebugMode = true ;
            %dbstop('if','error') ;
            
            self.IsAwake_ = isAwake ;
            self.IsHeaded_ = isHeaded ;
            self.DoUsePreferences_ = doUsePreferences ;
            
            % We only set up the sockets if we are the one true
            % WavesurferModel, and not some blasted pretender!
            % ("Pretenders" are created when we load a protocol from disk,
            % for instance.)
            if isAwake ,
                % Create the looper, refiller
                self.Looper_ = ws.Looper() ;
                self.Refiller_ = ws.Refiller() ;                                

%                 % Create the timer that runs during a run
%                 self.TheBigTimer_ = timer('ExecutionMode', 'fixedRate', ...
%                                           'BusyMode', 'drop', ...
%                                           'Period', 0.1, ...
%                                           'TimerFcn', @(timerObject, event)(self.handleTimerTick())) ;
                
                % Get the list of all device names, and cache it in our own
                % state
                self.AllDeviceNames_ = ws.getAllDeviceNamesFromHardware() ;
                self.syncDeviceResourceCountsFromDeviceNames_() ;
            end  % if isITheOneTrueWavesurfer
            
            % Initialize the fast protocols
            self.FastProtocols_ = cell(1,self.NFastProtocols) ;
            for i=1:self.NFastProtocols ,
                self.FastProtocols_{i} = ws.FastProtocol();
            end
            self.IndexOfSelectedFastProtocol_ = 1;
            
            % Create all subsystems
            self.Acquisition_ = ws.Acquisition() ;  % Acq subsystem doesn't need to know its parent, now.
            self.Stimulation_ = ws.Stimulation() ;  % Stim subsystem doesn't need to know its parent, now.
            self.Display_ = ws.Display() ;  % Stim subsystem doesn't need to know its parent, now.
            self.Triggering_ = ws.Triggering() ;  % Triggering subsystem doesn't need to know its parent, now.
            self.UserCodeManager_ = ws.UserCodeManager() ;
            self.Logging_ = ws.Logging() ;
            self.Ephys_ = ws.Ephys() ;
            
            % Create a list for methods to iterate when excercising the
            % subsystem API without needing to know all of the property
            % names.  Ephys must come before Acquisition to ensure the
            % right channels are enabled and the smart-electrode associated
            % gains are right, and before Display and Logging so that the
            % data values are correct.
            %self.Subsystems_ = {self.Ephys_, self.Acquisition_, self.Stimulation_, self.Display_, self.Triggering_, self.Logging_, self.UserCodeManager_};
            
            % The object is now initialized, but not very useful until a
            % device is specified.
            self.setState_('no_device') ;
            
            % Finally, set the device name to the first device name, if
            % there is one (and if we are the one true wavesurfer object)
            if isAwake ,
                % Set the device name to the first device
                allDeviceNames = self.AllDeviceNames ;
                if ~isempty(allDeviceNames) ,
                    self.PrimaryDeviceName = allDeviceNames{1} ;
                end
            end
            
            % Have to call this at object creation to set the override
            % correctly.
            self.overrideOrReleaseStimulusMapDurationAsNeeded_();            
            
            % Create a command connector (which will
            % remain disabled for now)
            % This needs to be done before loading the preferences, since it's the
            % CommandClient_ that takes care of notifying SI of preference file loading,
            % depending on whether we are yoked.
            self.CommandServer_ = ws.CommandServer(self) ;
            self.CommandClient_ = ws.CommandClient(self) ;
            if isAwake ,
                self.CommandServer_.IsEnabled = true ;
            end            
            
            % Set the protocol file name
            if isAwake ,
                %self.addStarterChannelsAndStimulusLibrary() ;
                self.loadProfileNameAndNames_() ;  % populates self.CurrentProfileName_                
                self.loadPreferences_(self.CurrentProfileName_) ;
                
                lastProtocolFilePath = self.LastProtocolFilePath_ ;  % just loaded from disk, typically 
                lastProtocolFileFolderPath = fileparts(lastProtocolFilePath) ;
                if isempty(lastProtocolFileFolderPath) || ~exist(lastProtocolFileFolderPath, 'dir') ,
                    protocolFileFolderPath = pwd() ;
                else
                    protocolFileFolderPath = lastProtocolFileFolderPath ;
                end
%                 index = 0 ;
%                 didFindGoodFileName = false ;
%                 while ~didFindGoodFileName ,
%                     index = index + 1 ;
%                     if index==1 ,
%                         putativeGoodFileName = 'untitled.wsp' ;
%                     else
%                         putativeGoodFileName = sprintf('untitled-%d.wsp', index) ;
%                     end
%                     putativeGoodFilePath = fullfile(protocolFileFolderPath, putativeGoodFileName) ;
%                     didFindGoodFileName = ~exist(putativeGoodFilePath, 'file') ;
%                 end                
                self.AbsoluteProtocolFileName_ = fullfile(protocolFileFolderPath, 'untitled.wsp') ;
                self.HasUserSpecifiedProtocolFileName_ = false ;
                self.DoesProtocolNeedSave_ = false ;
            end            
        end  % function
        
        function delete(self)            
            %fprintf('WavesurferModel::delete()\n');
            if self.IsAwake_ ,
                % Signal to others that we are going away
                %self.Looper_.frontendIsBeingDeleted() ;
                %self.Refiller_.frontendIsBeingDeleted() ;
                %self.IPCPublisher_.send('frontendIsBeingDeleted') ;
                
                try
                    % Stop, delete the timer
                    if ~isempty(self.TheBigTimer_)
                        if isvalid(self.TheBigTimer_) ,
                            stop(self.TheBigTimer_) ;
                            delete(self.TheBigTimer_) ;
                        end
                        self.TheBigTimer_ = [] ;
                    end                        
                catch me
                    fprintf('Unable to delete the timer in WSM delete() method.  Error was: %s\n', me.message) ;
                end
                
                % Try to save the preferences
                try
                    self.savePreferences_(self.CurrentProfileName_) ;  
                        % if the above fails, we don't try to save self.CurrentProfileName to the
                        % LastProfileName as stored on disk
                    self.saveLastProfileName_() ;
                catch me
                    fprintf('Unable to save the preferences in WSM delete() method.  Error was: %s\n', me.message) ;
                end                    
                
                % Close the sockets
                self.Looper_ = [] ;
                self.Refiller_ = [] ;
                %self.LooperIPCSubscriber_ = [] ;
                %self.RefillerIPCSubscriber_ = [] ;
                %self.LooperIPCRequester_ = [] ;
                %self.RefillerIPCRequester_ = [] ;
                %self.IPCPublisher_ = [] ;                                
                
                % If yoked, tell SI that we're quitting
                try
                    self.notifyScanImageThatWavesurferIsQuittingIfYoked_() ;
                catch exception  
                    % If something goes wrong, just ignore it, since this is delete() method   
                    fprintf('notifyScanImageThatWavesurferIsQuittingIfYoked_() threw in the WSM delete() method, with the error message:\n') ;
                    fprintf('%s', exception.message) ;                    
                end                    
            end
            % Delete CommandServer_, client
            ws.deleteIfValidHandle(self.CommandServer_) ;
            ws.deleteIfValidHandle(self.CommandClient_) ;
            ws.deleteIfValidHandle(self.UserCodeManager_) ;  % Causes user object to be explicitly deleted, if there is one
            self.UserCodeManager_ = [] ;            
        end  % function
        
        function debug(self) %#ok<MANU>
            keyboard
        end  % function        
                
        function play(self)
            % Start a run without recording data to disk.
            self.Logging_.IsEnabled = false ;
            self.run_() ;
        end  % function
        
        function record(self)
            % Start a run, recording data to disk.
            self.Logging_.IsEnabled = true ;
            self.run_() ;
        end  % function

        function stop(self)
            % Called when you press the "Stop" button in the UI, for
            % instance.  Stops the current run, if any.

            %fprintf('WavesurferModel::stop()\n');
            if isequal(self.State,'idle') , 
                % do nothing except re-sync the view to the model
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;                
            else
                % Actually stop the ongoing sweep
                %self.abortSweepAndRun_('user');
                %self.WasRunStoppedByUser_ = true ;
                %fprintf('About to publish "frontendWantsToStopRun"\n') ;
                %self.IPCPublisher_.send('frontendWantsToStopRun');  
                self.Looper_.stoppingRun(self.IsPerformingSweep) ;
                didStopEpisode = self.Refiller_.stoppingRun() ;
                if self.IsPerformingSweep_ ,
                    self.stopTheOngoingSweep_() ;
                end
                if didStopEpisode ,
                    self.callUserMethod_('stoppingEpisode') ;
                end
                self.stopTheOngoingRun_();
                
                %fprintf('Just published "frontendWantsToStopRun"\n') ;
                  % the looper gets this message and stops the run, then
                  % publishes 'looperStoppedRun'.
                  % similarly, the refiller gets this message, stops the
                  % run, then publishes 'refillerStoppedRun'.
            end
        end  % function
    end
    
    methods  % These are all the methods that get called in response to ZMQ messages
%         function result = samplesAcquired(self, scanIndex, rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData)
%             %fprintf('got data.  scanIndex: %d\n',scanIndex) ;
%             
%             % If we are not performing sweeps, just ignore.  This can
%             % happen after stopping a run and then starting another, where some old messages from the last run are
%             % still in the queue
%             if self.IsPerformingSweep_ ,
%                 %fprintf('About to call samplesAcquired_()\n') ;
%                 self.samplesAcquired_(scanIndex, rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData) ;
%             end
%             result = [] ;
%         end  % function
        
%         function result = looperCompletedSweep(self)
%             % Call by the Looper, via ZMQ pub-sub, when it has completed a sweep
%             %fprintf('WavesurferModel::looperCompletedSweep()\n');
%             %self.DidLooperCompleteSweep_ = true ;
%             self.DidLooperCompleteSweep_ = true ;
%             result = [] ;
%         end
        
%         function result = refillerCompletedEpisodes(self)
%             % Call by the Refiller, via ZMQ pub-sub, when it has completed
%             % all the episodes in the run
%             %fprintf('WavesurferModel::refillerCompletedRun()\n');            
%             self.DidRefillerCompleteEpisodes_ = true ;
%             %fprintf('Just did self.DidRefillerCompleteEpisodes_ = true\n');
%             %self.DidMostRecentSweepComplete_ = self.DidLooperCompleteSweep_ ;
%             result = [] ;
%         end
        
%         function result = looperStoppedRun(self)
%             % Call by the Looper, via ZMQ pub-sub, when it has stopped the
%             % run (in response to a frontendWantsToStopRun message)
%             %fprintf('WavesurferModel::looperStoppedRun()\n');
%             self.WasRunStoppedInLooper_ = true ;
%             %self.WasRunStopped_ = self.WasRunStoppedInRefiller_ ;
%             if self.WasRunStoppedInRefiller_ ,
%                 if self.IsPerformingSweep_ ,
%                     self.stopTheOngoingSweep_() ;
%                 end
%                 self.stopTheOngoingRun_();
%             end
% 
%             result = [] ;
%         end
        
%         function result = looperReadyForRunOrPerhapsNot(self, err)  %#ok<INUSL>
%             % Call by the Looper, via ZMQ pub-sub, when it has finished its
%             % preparations for a run.  Currrently does nothing, we just
%             % need a message to tell us it's OK to proceed.
%             result = err ;
%         end
%         
%         function result = looperReadyForSweep(self) %#ok<MANU>
%             % Call by the Looper, via ZMQ pub-sub, when it has finished its
%             % preparations for a sweep.  Currrently does nothing, we just
%             % need a message to tell us it's OK to proceed.
%             result = [] ;
%         end        
        
%         function result = looperIsAlive(self)  %#ok<MANU>
%             % Doesn't need to do anything
%             %fprintf('WavesurferModel::looperIsAlive()\n');
%             result = [] ;
%         end
        
%         function result = refillerReadyForRunOrPerhapsNot(self, err) %#ok<INUSL>
%             % Call by the Refiller, via ZMQ pub-sub, when it has finished its
%             % preparations for a run.  Currrently does nothing, we just
%             % need a message to tell us it's OK to proceed.
%             result = err ;
%         end
%         
%         function result = refillerReadyForSweep(self) %#ok<MANU>
%             % Call by the Refiller, via ZMQ pub-sub, when it has finished its
%             % preparations for a sweep.  Currrently does nothing, we just
%             % need a message to tell us it's OK to proceed.
%             result = [] ;
%         end        
        
%         function result = refillerStoppedRun(self)
%             % Call by the Looper, via ZMQ pub-sub, when it has stopped the
%             % run (in response to a frontendWantsToStopRun message)
%             %fprintf('WavesurferModel::refillerStoppedRun()\n');
%             self.WasRunStoppedInRefiller_ = true ;
%             %self.WasRunStopped_ = self.WasRunStoppedInLooper_ ;
%             if self.WasRunStoppedInLooper_ ,
%                 if self.IsPerformingSweep_ ,
%                     self.stopTheOngoingSweep_() ;
%                 end
%                 self.stopTheOngoingRun_();
%             end
%             
%             result = [] ;
%         end
        
%         function result = refillerIsAlive(self)  %#ok<MANU>
%             % Doesn't need to do anything
%             %fprintf('WavesurferModel::refillerIsAlive()\n');
%             result = [] ;
%         end
        
%         function result = looperDidReleaseTimedHardwareResources(self) %#ok<MANU>
%             result = [] ;
%         end
%         
%         function result = refillerDidReleaseTimedHardwareResources(self) %#ok<MANU>
%             result = [] ;
%         end       
        
%         function result = gotMessageHeyRefillerIsDigitalOutputTimedWasSetInFrontend(self) %#ok<MANU>
%             result = [] ;
%         end       
    end  % ZMQ methods block
    
    methods
        function value = get.VersionString(self)  %#ok<MANU>
            value = ws.versionString() ;
        end  % function
        
        function value=get.State(self)
            value=self.State_;
        end  % function
        
%         function value=get.NextSweepIndex(self)
%             % This is a pass-through method to get the NextSweepIndex from
%             % the Logging subsystem, where it is actually stored.
%             if isempty(self.Logging) || ~isvalid(self.Logging),
%                 value=[];
%             else
%                 value=self.Logging_.NextSweepIndex;
%             end
%         end  % function
        
%         function out = get.Acquisition(self)
%             out = self.Acquisition_ ;
%         end
%         
%         function out = get.Stimulation(self)
%             out = self.Stimulation_ ;
%         end
%         
%         function out = get.Triggering(self)
%             out = self.Triggering_ ;
%         end
%         
%         function out = get.UserCodeManager(self)
%             out = self.UserCodeManager_ ;
%         end
% 
%         function out = get.Display(self)
%             out = self.Display_ ;
%         end
%         
%         function out = get.Logging(self)
%             out = self.Logging_ ;
%         end
        
%         function out = get.Ephys(self)
%             out = self.Ephys_ ;
%         end
        
        function out = get.NSweepsCompletedInThisRun(self)
            out = self.NSweepsCompletedInThisRun_ ;
        end
        
        function out = get.HasUserSpecifiedProtocolFileName(self)
            out = self.HasUserSpecifiedProtocolFileName_ ;
        end
        
        function out = get.AbsoluteProtocolFileName(self)
            out = self.AbsoluteProtocolFileName_ ;
        end
        
        function out = get.HasUserSpecifiedUserSettingsFileName(self)
            out = self.HasUserSpecifiedUserSettingsFileName_ ;
        end
        
        function out = get.AbsoluteUserSettingsFileName(self)
            out = self.AbsoluteUserSettingsFileName_ ;
        end
        
        function out = get.IndexOfSelectedFastProtocol(self)
            out = self.IndexOfSelectedFastProtocol_ ;
        end
        
        function set.IndexOfSelectedFastProtocol(self, newValue)
            if isnumeric(newValue) && isscalar(newValue) && isreal(newValue) && 1<=newValue && newValue<=self.NFastProtocols && (round(newValue)==newValue) ,
                self.IndexOfSelectedFastProtocol_ = double(newValue) ;
            else
                error('ws:invalidPropertyValue', ...
                      'IndexOfSelectedFastProtocol (scalar) positive integer between 1 and NFastProtocols, inclusive');
            end
        end
        
        function val = get.NSweepsPerRun(self)
            if self.AreSweepsContinuous ,
                val = 1;
            else
                %val = self.Triggering_.SweepTrigger.Source.RepeatCount;
                val = self.NSweepsPerRun_;
            end
        end  % function
        
        function set.NSweepsPerRun(self, newValue)
            if isscalar(newValue) && isnumeric(newValue) && isreal(newValue) && newValue>=1 && ...
                    (round(newValue)==newValue || isinf(newValue)) ,
                % If get here, value is a valid value for this prop
                if self.AreSweepsFiniteDuration ,
                    %self.Triggering_.willSetNSweepsPerRun();
                    self.NSweepsPerRun_ = double(newValue) ;
                    self.didSetNSweepsPerRun_(self.NSweepsPerRun_) ;
                    self.DoesProtocolNeedSave_ = true ;
                else
                    self.broadcast('UpdateMain');
                    self.broadcast('UpdateGeneral') ;
                    self.broadcast('UpdateChannels') ;
                    error('ws:invalidPropertyValue', ...
                          'NSweepsPerRun cannot be set when sweeps are continuous') ;

                end
            else
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;                
                error('ws:invalidPropertyValue', ...
                      'NSweepsPerRun must be a (scalar) positive integer, or inf') ;       
            end
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function out = get.SweepDurationIfFinite(self)
            out = self.SweepDurationIfFinite_ ;
        end  % function
        
        function set.SweepDurationIfFinite(self, value)
            %fprintf('Acquisition::set.Duration()\n');
            if isnumeric(value) && isscalar(value) && isfinite(value) && value>0 ,
                valueToSet = max(value,0.1);
                %self.willSetSweepDurationIfFinite();
                self.SweepDurationIfFinite_ = valueToSet;
                self.overrideOrReleaseStimulusMapDurationAsNeeded_();
                if self.IsXSpanSlavedToAcquistionDuration ,
                    self.syncTraces_() ;
                    self.broadcast('UpdateTraces') ;
                end
                self.broadcast('DidSetXSpan') ;
                self.DoesProtocolNeedSave_ = true ;
            else
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;                
                error('ws:invalidPropertyValue', ...
                      'SweepDurationIfFinite must be a (scalar) positive finite value');
            end
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function value = get.SweepDuration(self)
            if self.AreSweepsContinuous ,
                value=inf;
            else
                value=self.SweepDurationIfFinite_ ;
            end
        end  % function
        
        function set.SweepDuration(self, newValue)
            % Check value and set if valid
            if isnumeric(newValue) && isscalar(newValue) && ~isnan(newValue) && newValue>0 ,
                % If get here, newValue is a valid value for this prop
                if isfinite(newValue) ,
                    self.AreSweepsFiniteDuration = true ;
                    self.SweepDurationIfFinite = newValue ;
                else                        
                    self.AreSweepsContinuous = true ;
                end                        
                self.DoesProtocolNeedSave_ = true ;
            else
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;
                error('ws:invalidPropertyValue', ...
                      'SweepDuration must be a (scalar) positive value');
            end
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function value=get.AreSweepsFiniteDuration(self)
            value = self.AreSweepsFiniteDuration_ ;
        end
        
        function set.AreSweepsFiniteDuration(self,newValue)
            %fprintf('inside set.AreSweepsFiniteDuration.  self.AreSweepsFiniteDuration_: %d\n', self.AreSweepsFiniteDuration_);
            %newValue            
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                %fprintf('setting self.AreSweepsFiniteDuration_ to %d\n',logical(newValue));
                %self.willSetAreSweepsFiniteDuration();
                self.AreSweepsFiniteDuration_=logical(newValue);
                %self.AreSweepsContinuous=nan.The;
                %self.NSweepsPerRun=nan.The;
                %self.SweepDuration=nan.The;
                self.DoesProtocolNeedSave_ = true ;                
                self.overrideOrReleaseStimulusMapDurationAsNeeded_();
                self.didSetAreSweepsFiniteDuration_(self.AreSweepsFiniteDuration_, self.NSweepsPerRun_);
            end
            %self.broadcast('DidSetAreSweepsFiniteDurationOrContinuous');            
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end
        
%         function set.ReferenceClockSource(self, newValue)
%             if ws.isString(newValue) && ismember(newValue, self.AvailableReferenceClockSources) ,
%                 self.ReferenceClockSource_ = newValue ;
%                 %self.syncReferenceClockRateFromDeviceNameAndAvailableReferenceClockSourcesEtc_() ;
%                 isNewValueValid = true ;
%             else
%                 isNewValueValid = false ;
%             end                
%             self.broadcast('Update');
%             if ~isNewValueValid ,
%                 error('ws:invalidPropertyValue', ...
%                       'ReferenceClockSource must be a string, and equal to some element of AvailableReferenceClockSources');
%             end
%         end
        
        function value = get.AreSweepsContinuous(self)
            value = ~self.AreSweepsFiniteDuration_ ;
        end
        
        function set.AreSweepsContinuous(self,newValue)
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                self.AreSweepsFiniteDuration=~logical(newValue);
            end
        end        
    end
    
    methods (Access=protected)    
        function overrideOrReleaseStimulusMapDurationAsNeeded_(self)
            isSweepBased = self.AreSweepsFiniteDuration ;
            doesStimulusUseAcquisitionTriggerScheme = self.StimulationUsesAcquisitionTrigger ;
            %isStimulationTriggerIdenticalToAcquisitionTrigger = self.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
                % Is this the best design?  Seems maybe over-complicated.
                % Future Adam: If you change to using
                % isStimulationTriggerIdenticalToAcquisitionTrigger() to
                % determine whether map durations are overridden, leave a
                % note as to why you did that, preferably with a particular
                % "failure mode" in the current scheme.
            if isSweepBased && doesStimulusUseAcquisitionTriggerScheme ,
                self.Stimulation_.overrideStimulusLibraryMapDuration(self.SweepDuration) ;
            else
                self.Stimulation_.releaseStimulusLibraryMapDuration() ;
            end 
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateStimulusPreview') ;
        end  
        
        function notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_(self)
            % Called after setting an analog channel unit or scale (directly), to
            % notify other systems of the change.  Here "other" means "other than the
            % Acquisition or Stimulation subsystem where the change was made".
            %self.Display_.didSetAnalogChannelUnitsOrScales() ;
            self.Ephys_.didSetAnalogChannelUnitsOrScales() ;
            self.syncTraces_() ;
            self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;
            self.broadcast('UpdateTestPulser') ;
        end

    end  % protected methods block
    
    methods
%         function electrodesRemoved(self)
%             % Called by the Ephys to notify that one or more electrodes
%             % was removed
%             % Currently, tells Acquisition and Stimulation about the change.
%             self.didSetAnalogChannelUnitsOrScales() ;      
%         end       

%         function setSingleAIChannelTerminalID(self, i, newValue)
%             self.Acquisition_.setSingleAnalogTerminalID_(i, newValue) ;
%             self.didSetAnalogInputTerminalID();
%         end

        function setSingleAIChannelTerminalName(self, i, terminalName)
            % newValue should be a string like 'AI4' or 'AI11'
            terminalIDAsString = terminalName(3:end) ;
            terminalID = str2double(terminalIDAsString) ;
            self.setSingleAIChannelTerminalID(i, terminalID) ;
        end

        function setSingleAIChannelTerminalID(self, i, newValue)
            self.Acquisition_.setSingleAnalogTerminalID(i, newValue) ;
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsAIChannelTerminalOvercommitted_() ;
            %self.Display_.didSetAnalogInputTerminalID() ;
            %self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateChannels') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end

        function setSingleDIChannelTerminalName(self, iChannel, terminalName)
            terminalID = ws.dioTerminalIDFromName(terminalName) ;
            self.setSingleDIChannelTerminalID(iChannel, terminalID) ;
        end
        
        function setSingleDIChannelTerminalID(self, iChannel, terminalID)
            wasSet = self.Acquisition_.setSingleDigitalTerminalID_(iChannel, terminalID) ;            
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsDIOChannelTerminalOvercommitted_() ;
            %self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateChannels') ;
            self.broadcast('DidMaybeChangeProtocol') ;
            if wasSet ,
                %value = self.Acquisition_.DigitalTerminalIDs(iChannel) ;  % value is possibly normalized, terminalID is not
                self.Looper_.singleDigitalInputTerminalIDWasSetInFrontend(self.PrimaryDeviceName, ...
                                                                          self.IsPrimaryDeviceAPXIDevice, ...
                                                                          self.DOChannelTerminalIDs, ...
                                                                          self.IsDOChannelTimed, ...
                                                                          self.DOChannelStateIfUntimed, ...
                                                                          self.IsDOChannelTerminalOvercommitted) ;
                self.Refiller_.singleDigitalInputTerminalIDWasSetInFrontend(self.IsDOChannelTerminalOvercommitted) ;                
                %self.IPCPublisher_.send('singleDigitalInputTerminalIDWasSetInFrontend', ...
                %                        self.IsDOChannelTerminalOvercommitted ) ;
            end            
        end
        
%         function didSetAnalogInputChannelName(self, didSucceed, oldValue, newValue)
%             display=self.Display_;
%             if ~isempty(display)
%                 display.didSetAnalogInputChannelName(didSucceed, oldValue, newValue);
%             end            
%             ephys=self.Ephys_;
%             if ~isempty(ephys)
%                 ephys.didSetAnalogInputChannelName(didSucceed, oldValue, newValue);
%             end            
%             self.broadcast('UpdateChannels') ;
%         end
        
%         function didSetDigitalInputChannelName(self, didSucceed, oldValue, newValue)
%             self.Display_.didSetDigitalInputChannelName(didSucceed, oldValue, newValue);
%             self.broadcast('UpdateChannels') ;
%         end
        
%         function didSetAnalogOutputTerminalID(self)
%             self.syncIsAOChannelTerminalOvercommitted_() ;
%             self.broadcast('UpdateChannels') ;
%         end
        
        function setSingleAOChannelName(self, channelIndex, newValue)
            nAOChannels = self.Stimulation_.NAnalogChannels ;
            if ws.isIndex(channelIndex) && 1<=channelIndex && channelIndex<=nAOChannels ,
                oldValue = self.Stimulation_.AnalogChannelNames{channelIndex} ;
                if ws.isString(newValue) && ~isempty(newValue) && ~ismember(newValue,self.AllChannelNames) ,
                     self.Stimulation_.setSingleAnalogChannelName(channelIndex, newValue) ;
                     self.DoesProtocolNeedSave_ = true ;
                     didSucceed = true ;
                else
                    didSucceed = false;
                end
            else
                didSucceed = false ;
            end
            ephys = self.Ephys_ ;
            if ~isempty(ephys) ,
                ephys.didSetAnalogOutputChannelName(didSucceed, oldValue, newValue);
            end            
            self.broadcast('UpdateElectrodeManager') ;
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateChannels') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function 
        
%         function didSetAnalogOutputChannelName(self, didSucceed, oldValue, newValue)
%             ephys=self.Ephys;
%             if ~isempty(ephys)
%                 ephys.didSetAnalogOutputChannelName(didSucceed, oldValue, newValue);
%             end            
%             self.broadcast('UpdateChannels') ;
%         end
        
%         function didSetDigitalOutputChannelName(self, didSucceed, oldValue, newValue) %#ok<INUSD>
%             self.broadcast('UpdateChannels') ;
%         end

        function setSingleDOChannelName(self, channelIndex, newValue)
            nDOChannels = self.Stimulation_.NDigitalChannels ;
            if ws.isIndex(channelIndex) && 1<=channelIndex && channelIndex<=nDOChannels ,
                %oldValue = self.Stimuluation_.DigitalChannelNames{channelIndex} ;
                if ws.isString(newValue) && ~isempty(newValue) && ~ismember(newValue,self.AllChannelNames) ,
                     self.Stimulation_.setSingleDigitalChannelName(channelIndex, newValue) ;
                     self.DoesProtocolNeedSave_ = true ;
                     %didSucceed = true ;
                else
                    %didSucceed = false;
                end
            else
                %didSucceed = false ;
            end
%             ephys = self.Ephys ;
%             if ~isempty(ephys) ,
%                 ephys.didSetAnalogOutputChannelName(didSucceed, oldValue, newValue);
%             end            
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateChannels') ;          
            self.broadcast('DidMaybeChangeProtocol') ;            
        end  % function 
    end  % public methods block
       
    methods (Access=protected)
        function setIsYokedToScanImage_(self, newValue)
            err = [] ;
            
            % Set the value by enabling/disabling the command connector
            if islogical(newValue) && isscalar(newValue) ,
                try
                    self.CommandClient_.IsEnabled = (self.IsAwake_ && newValue) ;
                catch me
                    err = me ;
                end                
                %self.IsYokedToScanImage_ = self.CommandClient_.IsEnabled ;
            end            

            % Do an update
            self.broadcast('UpdateIsYokedToScanImage');

            % Throw any error that occured
            if ~isempty(err) ,
                throw(err) ;
            end
        end  % function
    end  % protected methods block
    
    methods
        function setIsYokedToScanImageForTesting_(self, newValue)
            self.setIsYokedToScanImage_(newValue) ;
        end        
        
        function value=get.IsYokedToScanImage(self)
            value = self.CommandClient_.IsEnabled ;
        end  % function        
        
        function value=get.IsITheOneTrueWavesurferModel(self)
            % This is deprecated
            value = self.IsAwake_ ;
        end  % function        
        
        function value=get.IsAwake(self)
            value = self.IsAwake_ ;
        end  % function        
        
%         function willPerformTestPulse(self)
%             % Called by the TestPulserModel to inform the WavesurferModel that
%             % it is about to start test pulsing.
%             
%             % I think the main thing we want to do here is to change the
%             % Wavesurfer mode to TestPulsing.
%             if isequal(self.State,'idle') ,
%                 self.setState_('test_pulsing');
%             end
%         end
        
%         function didPerformTestPulse(self)
%             % Called by the TestPulserModel to inform the WavesurferModel that
%             % it has just finished test pulsing.
%             
%             if isequal(self.State,'test_pulsing') ,
%                 self.setState_('idle');
%             end
%         end  % function
        
%         function didAbortTestPulse(self)
%             % Called by the TestPulserModel when a problem arises during test
%             % pulsing, that (hopefully) the TestPulseModel has been able to
%             % gracefully recover from.
%             
%             if isequal(self.State,'test_pulsing') ,
%                 self.setState_('idle');
%             end
%         end  % function
        
    end  % public methods block
    
    methods (Access=protected)
        function releaseTimedHardwareResources_(self)
            %self.Acquisition_.releaseHardwareResources();
            %self.Stimulation_.releaseHardwareResources();
            self.Triggering_.releaseTimedHardwareResources();
            %self.Ephys_.releaseTimedHardwareResources();
        end

        function didSetAreSweepsFiniteDuration_(self, areSweepsFiniteDuration, nSweepsPerRun) %#ok<INUSL>
            self.Triggering_.didSetAreSweepsFiniteDuration(nSweepsPerRun);
            self.syncTraces_() ;
            self.broadcast('UpdateTriggering') ;
            self.broadcast('UpdateTraces') ;
            self.broadcast('DidSetXSpan') ;
        end        

        function releaseTimedHardwareResourcesOfAllProcesses_(self)
            % Release our own hardware resources, and also tell the
            % satellites to do so.
            self.releaseTimedHardwareResources_() ;
            %self.IPCPublisher_.send('satellitesReleaseTimedHardwareResources') ;
            
            % Wait for the looper to respond
            %timeout = 10 ;  % s
            err = self.Looper_.releaseTimedHardwareResources() ;
            %self.LooperIPCRequester_.send('releaseTimedHardwareResources') ;
            %err = self.LooperIPCRequester_.waitForResponse(timeout, 'releaseTimedHardwareResources') ;
            if ~isempty(err) ,
                % Something went wrong
                throw(err);
            end
            
            % Wait for the refiller to respond
            err = self.Refiller_.releaseTimedHardwareResources() ;
            %self.RefillerIPCRequester_.send('releaseTimedHardwareResources') ;
            %err = self.RefillerIPCRequester_.waitForResponse(timeout, 'releaseTimedHardwareResources') ;
            if ~isempty(err) ,
                % Something went wrong
                throw(err);
            end            
        end  % function        
    end
    
    methods (Access = protected)
        function setState_(self,newValue)
            %self.broadcast('WillSetState');
            if ws.isAnApplicationState(newValue) ,
                if ~isequal(self.State_,newValue) ,
                    self.State_ = newValue ;
                end
            end
            self.broadcast('DidSetState');
        end  % function
        
        function run_(self)
            %fprintf('WavesurferModel::run_()\n');     
            
            % Can't run unless we are currently idle
            if ~isequal(self.State,'idle') ,
                return
            end

            % Change the readiness (this changes the pointer in the view)
            self.changeReadiness_(-1);
            
            % If yoked to scanimage, write to the command file, wait for a
            % response
            if self.AreSweepsFiniteDuration ,
                try
                    self.commandScanImageToStartLoopIfYoked_() ;
                catch excp
                    self.abortOngoingRun_() ;
                    self.changeReadiness_(+1) ;
                    rethrow(excp) ;
                end
            end
            
            if ~isempty(self.CommandServer_) ,
                self.CommandServer_.clearIsProcessingIncomingCommand_() ;  % awkward workaround because self.play() and self.record() are blocking calls
            end
            
            % Initialize the sweep counter, etc.
            self.NSweepsCompletedInThisRun_ = 0 ;
            self.AreAllSweepsCompleted_ = (self.NSweepsCompletedInThisRun_>=self.NSweepsPerRun) ;            
            %self.DidRefillerCompleteEpisodes_ = false ;
            %fprintf('Just did self.DidRefillerCompleteEpisodes_ = false\n');
            
            % Call the user method, if any
            self.callUserMethod_('startingRun');  
                        
            % Tell all the subsystems except the logging subsystem to prepare for the run
            % The logging subsystem has to wait until we obtain the analog
            % scaling coefficients from the Looper.
            try
                if self.DoTrodeUpdateBeforeRun ,
                    self.updateSmartElectrodeGainsAndModes() ;
                end
                self.Acquisition_.startingRun(self.AreSweepsContinuous, self.AreSweepsFiniteDuration, self.SweepDuration) ;
                if self.Stimulation_.IsEnabled ,
                    self.Stimulation_.startingRun() ;
                end
                if self.Display_.IsEnabled ,
                    self.XOffset = 0 ;
                    self.Display_.startingRun(self.XSpan, self.SweepDuration) ;
                    self.syncTraces_() ;
                    self.broadcast('UpdateTraces') ;
                end
                primaryDeviceName = self.PrimaryDeviceName ;
                isPrimaryDeviceAPXIDevice = self.IsPrimaryDeviceAPXIDevice ;
                self.Triggering_.startingRun(primaryDeviceName, isPrimaryDeviceAPXIDevice) ;
                self.UserCodeManager_.startingRun() ;
            catch me
                % Something went wrong
                self.abortOngoingRun_();
                self.changeReadiness_(+1);
                me.rethrow();
            end
            
            % Determine the keystone tasks for acq and stim
            [acquisitionKeystoneTaskType, acquisitionKeystoneTaskDeviceName, ...
             stimulationKeystoneTaskType, stimulationKeystoneTaskDeviceName] = ...
               self.determineKeystoneTasks_() ;
            
            % Tell the looper to start
            try
                coeffsAndClock = ...
                    self.Looper_.startingRun(acquisitionKeystoneTaskType, ...
                                             acquisitionKeystoneTaskDeviceName, ...
                                             self.IsAIChannelActive, ...
                                             self.IsDIChannelActive, ...
                                             self.AcquisitionSampleRate, ...
                                             self.SweepDuration, ...
                                             self.DataCacheDurationWhenContinuous, ...
                                             self.AIChannelDeviceNames, ...
                                             self.AIChannelTerminalIDs, ...
                                             self.PrimaryDeviceName, ...
                                             self.IsPrimaryDeviceAPXIDevice, ...
                                             self.acquisitionTriggerProperty('DeviceName'), ...
                                             self.acquisitionTriggerProperty('PFIID'), ...
                                             self.acquisitionTriggerProperty('Edge'), ...
                                             self.DIChannelTerminalIDs) ;
            catch err
                self.abortOngoingRun_();
                self.changeReadiness_(+1);
                rethrow(err) ;
            end            
            
            % Tell the refiller to start
            try
                stimulationTriggerClass = self.stimulationTriggerProperty('class') ;
                if isequal(stimulationTriggerClass, 'ws.CounterTrigger') ,                    
                    stimulationTriggerRepeatCount = self.stimulationTriggerProperty('RepeatCount') ;
                else
                    stimulationTriggerRepeatCount = nan ;
                end                    
                stimulationTriggerDeviceName = self.stimulationTriggerProperty('DeviceName') ;
                stimulationTriggerPFIID = self.stimulationTriggerProperty('PFIID') ;
                stimulationTriggerEdge = self.stimulationTriggerProperty('Edge') ;
                isStimulationTriggerIdenticalToAcquisitionTrigger = self.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
                self.Refiller_.startingRun(stimulationKeystoneTaskType, ...
                                           stimulationKeystoneTaskDeviceName, ...
                                           self.IsStimulationEnabled, ...
                                           stimulationTriggerClass, ...
                                           self.NSweepsPerRun, ...
                                           stimulationTriggerRepeatCount, ...
                                           isStimulationTriggerIdenticalToAcquisitionTrigger, ...
                                           self.AOChannelDeviceNames, ...
                                           self.IsAOChannelTerminalOvercommitted, ...
                                           self.AOChannelTerminalIDs, ...
                                           self.PrimaryDeviceName, ...
                                           self.IsPrimaryDeviceAPXIDevice, ...
                                           self.IsDOChannelTerminalOvercommitted, ...
                                           self.IsDOChannelTimed, ...
                                           self.DOChannelTerminalIDs, ...
                                           self.StimulationSampleRate, ...
                                           stimulationTriggerDeviceName, ...
                                           stimulationTriggerPFIID, ...
                                           stimulationTriggerEdge) ;
                % (Maybe) start an episode, which will wait for a trigger to *really*
                % start.
                % It's important to do this here, so that we get ready to
                % receive the first stim trigger *before* we tell the frontend
                % that we're ready for the run.
                if self.Refiller_.NEpisodesPerRun > 0 ,
                    if ~isStimulationTriggerIdenticalToAcquisitionTrigger ,
                        self.callUserMethod_('startingEpisode') ;
                        [aoData, doData] = self.getStimulationData(self.Refiller_.NEpisodesCompletedSoFarThisRun+1) ;                        
                        self.Refiller_.startEpisode(aoData, doData) ;
                    end
                end
            catch err
                self.abortOngoingRun_();
                self.changeReadiness_(+1);
                rethrow(err) ;                
            end            
            
%             % Create the timer
%             self.TheBigTimer_ = timer('ExecutionMode', 'fixedRate', ...
%                                       'BusyMode', 'drop', ...
%                                       'Period', 0.1, ...
%                                       'TimerFcn', @(timerObject, event)(self.handleTimerTick()), ...
%                                       'ErrorFcn', @(timerObject, event)(self.handleTimerError(event))) ;                                      
            
            % Stash the analog scaling coefficients (have to do this now,
            % instead of in Acquisiton.startingRun(), b/c we get them from
            % the looper
            analogScalingCoefficients = coeffsAndClock.ScalingCoefficients ;
            clockAtRunStartTic = coeffsAndClock.ClockAtRunStartTic ;
            self.Acquisition_.cacheAnalogScalingCoefficients(analogScalingCoefficients) ;
            self.ClockAtRunStart_ = clockAtRunStartTic ;  % store the value returned from the looper
            
            % Now tell the logging subsystem that a run is about to start,
            % since the analog scaling coeffs have been set
            try
                logging = self.Logging_ ;
                if logging.IsEnabled ,
                    headerStruct = ws.encodeForHeader(self) ;
                    logging.startingRun(self.NextRunAbsoluteFileName, self.IsAIChannelActive, self.ExpectedSweepScanCount, ...
                                        self.AcquisitionSampleRate, self.AreSweepsFiniteDuration, headerStruct);
                end
            catch me
                % Something went wrong
                self.abortOngoingRun_();
                self.changeReadiness_(+1);
                me.rethrow();
            end
            
            % Change our own state to running
            self.NTimesDataAvailableCalledSinceRunStart_=0;  % Have to do this now so that progress bar doesn't start in old position when continuous acq
            self.setState_('running') ;
            self.IsPerformingSweep_ = false ;  % the first sweep starts later, if at all
            self.IsPerformingRun_ = true ;
            
            % Handle timing stuff
            self.TimeOfLastWillPerformSweep_ = [] ;
            self.FromRunStartTicId_ = tic() ;
            rawUpdateDt = 1/self.Display_.UpdateRate ;  % s
            updateDt = min(rawUpdateDt,self.SweepDuration);  % s
            desiredNScansPerUpdate = max(1,round(updateDt*self.AcquisitionSampleRate)) ;  % don't want this to be zero!
            self.DesiredNScansPerUpdate_ = desiredNScansPerUpdate ;
            self.NScansPerUpdate_ = self.DesiredNScansPerUpdate_ ;  % at start, will be modified depending on how long stuff takes
            
            % Set up the samples buffer
            bufferSizeInScans = 30*self.DesiredNScansPerUpdate_ ;
            self.SamplesBuffer_ = ws.SamplesBuffer(self.Acquisition_.NActiveAnalogChannels, ...
                                                   self.Acquisition_.NActiveDigitalChannels, ...
                                                   bufferSizeInScans) ;
            
            self.changeReadiness_(+1);  % do this now to give user hint that they can press stop during run...

            %
            % Move on to the main within-run loop
            %
            %self.DidMostRecentSweepComplete_ = true ;  % Set this to true just so the first sweep gets started
            self.DrawnowTicId_ = tic() ;  % Need this for timing between drawnow()'s
            self.TimeOfLastDrawnow_ = toc(self.DrawnowTicId_) ;  % we don't really do a a drawnow() here, but need to init
            %didPerformFinalDrawnow = false ;
            %for iSweep = 1:self.NSweepsPerRun ,
            % Can't use a for loop b/c self.NSweepsPerRun can be Inf
            %while iSweep<=self.NSweepsPerRun && ~self.RefillerCompletedSweep ,            
            %self.WasRunStoppedInLooper_ = false ;
            %self.WasRunStoppedInRefiller_ = false ;
            %self.WasRunStopped_ = false ;
            %self.NSweepsPerformed_ = 0 ;  % in this run
            self.DidAnySweepFailToCompleteSoFar_ = false ;
            %didLastSweepComplete = true ;  % It is convenient to pretend the zeroth sweep completed successfully
            
            % Create the timer that runs during a run
            self.TheBigTimer_ = timer('ExecutionMode', 'fixedRate', ...
                                      'BusyMode', 'drop', ...
                                      'Period', 0.05, ...
                                      'TimerFcn', @(timerObject, event)(self.handleTimerTick())) ;
            
            % Finally, start the timer
            start(self.TheBigTimer_) ;
        end  % run_() function
        
        function openSweep_(self)
            % time between subsequent calls to this, etc
            t=toc(self.FromRunStartTicId_);
            self.TimeOfLastWillPerformSweep_=t;
            
            % clear timing stuff that is strictly within-sweep
            %self.TimeOfLastSamplesAcquired_=[];            
            
            % Reset the sample count for the sweep
            %self.NScansAcquiredSoFarThisSweep_ = 0;
            
            % update the current time
            self.t_=0;            
            
            % Pretend that we last polled at time 0
            self.TimeOfLastPollInSweep_ = 0 ;  % s 
            
            % As of the next line, the wavesurferModel is offically
            % performing a sweep
            % Actually, we'll wait until a bit later to set this true
            %self.IsPerformingSweep_ = true ;
            
            % Call user functions 
            self.callUserMethod_('startingSweep');            
            
            % Call startingSweep() on all the enabled subsystems
%             for idx = 1:numel(self.Subsystems_)
%                 if self.Subsystems_{idx}.IsEnabled ,
%                     self.Subsystems_{idx}.startingSweep();
%                 end
%             end
            %self.Ephys_.startingSweep() ;
            self.Acquisition_.startingSweep() ;
            %self.Stimulation_.startingSweep() ;
            if self.Display_.IsEnabled ,
                self.Display_.startingSweep() ;
            end
            %self.Triggering_.startingSweep() ;
            if self.Logging_.IsEnabled ,
                self.Logging_.startingSweep() ;
            end
            %self.UserCodeManager_.startingSweep() ;

            % Notify the refiller that we're starting a sweep, wait for the refiller to respond
            if self.Stimulation_.IsEnabled && (self.StimulationTriggerIndex==self.AcquisitionTriggerIndex) ,
                try
                    self.callUserMethod_('startingEpisode') ;
                    [aoData, doData] = self.getStimulationData(self.Refiller_.NEpisodesCompletedSoFarThisRun+1) ;
                    self.Refiller_.startEpisode(aoData, doData) ;
                    err = [] ;
                catch err
                end
                %self.RefillerIPCRequester_.send('startingSweep', self.NSweepsCompletedInThisRun_+1) ;
                %timeout = 11 ;  % s
                %err = self.RefillerIPCRequester_.waitForResponse(timeout, 'startingSweep') ;
                if ~isempty(err) ,
                    % Something went wrong
                    self.abortOngoingRun_();
                    self.changeReadiness_(+1);
                    throw(err);
                end
            end
            
            % Notify the looper that we're starting a sweep, wait for the looper to respond
            try
                self.Looper_.startingSweep(self.NSweepsCompletedInThisRun_+1) ;
                err = [] ;
            catch err
            end
            %self.LooperIPCRequester_.send('startingSweep', self.NSweepsCompletedInThisRun_+1) ;
            %timeout = 12 ;  % s
            %err = self.LooperIPCRequester_.waitForResponse(timeout, 'startingSweep') ;
            if ~isempty(err) 
                % Something went wrong
                self.abortOngoingRun_();
                self.changeReadiness_(+1);
                throw(err);
            end

            % Set the sweep timer
            self.FromSweepStartTicId_=tic();

            % Pulse the master trigger to start the sweep!
            %fprintf('About to pulse the master trigger!\n');
            %self.DidLooperCompleteSweep_ = false ;
            %self.DidAnySweepFailToCompleteSoFar_ = true ;
            self.IsPerformingSweep_ = true ;  
                % have to wait until now to set this true, so that old
                % samplesAcquired messages are ignored when we're
                % waiting for looperReadyForSweep and
                % refillerReadyForSweep messages above.
            self.Triggering_.pulseBuiltinTrigger();

            % Now poll
            %self.DidLooperCompleteSweep_ = false ;
            %self.DidRefillerCompleteSweep_ = ~self.Stimulation_.IsEnabled ;  % if stim subsystem is disabled, then the refiller is automatically done
        end

        function handleTimerTick(self)
            % Called every so often once data acquition has started            
            if ~self.AllowTimerCallback_ ,
                return
            end
            try
                if ~self.DidAnySweepFailToCompleteSoFar_ && ~(self.AreAllSweepsCompleted_ && self.Refiller_.DidCompleteEpisodes) ,
                    %fprintf('wasRunStopped: %d\n', self.WasRunStopped_) ;
                    if self.IsPerformingSweep_ ,
                        if self.Looper_.DidCompleteSweep ,
                            self.completeTheOngoingSweep_() ;
                        else
                            %fprintf('At top of within-sweep loop...\n') ;
                            timeSinceSweepStart = toc(self.FromSweepStartTicId_) ;
                            self.performOneLooperIterationDuringOngoingSweep_(timeSinceSweepStart, ...
                                                                              self.FromRunStartTicId_) ;
                            isStimulationTriggerIdenticalToAcquisitionTrigger = self.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
                            self.performOneRefillerIterationDuringOngoingRun_(isStimulationTriggerIdenticalToAcquisitionTrigger) ;
                            % do a drawnow() if it's been too long...
                            %timeSinceLastDrawNow = toc(self.DrawnowTicId_) - self.TimeOfLastDrawnow_ ;
                            %if timeSinceLastDrawNow > 0.1 ,  % 0.1 s, hence 10 Hz
                            drawnow('limitrate', 'nocallbacks') ;
                            %    self.TimeOfLastDrawnow_ = toc(self.DrawnowTicId_) ;
                            %end                    
                        end
                    else
                        % We are not currently performing a sweep, so check if we need to start one
                        if self.AreAllSweepsCompleted_ ,
                            % All sweeps are were performed, but the refiller must not be done yet if we got here
                            % Keep checking messages so we know when the
                            % refiller is done.  Also keep listening for looper
                            % messages, although I'm not sure we need to...
                            %fprintf('About to check for messages after completing all sweeps\n');
                            isStimulationTriggerIdenticalToAcquisitionTrigger = self.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
                            self.performOneRefillerIterationDuringOngoingRun_(isStimulationTriggerIdenticalToAcquisitionTrigger) ;                            
                            %fprintf('Check for messages after completing all sweeps\n');
                            % do a drawnow() if it's been too long...
                            %timeSinceLastDrawNow = toc(self.DrawnowTicId_) - self.TimeOfLastDrawnow_ ;
                            %if timeSinceLastDrawNow > 0.1 ,  % 0.1 s, hence 10 Hz
                            drawnow('limitrate', 'nocallbacks') ;
                            %    self.TimeOfLastDrawnow_ = toc(self.DrawnowTicId_) ;
                            %end
                        else                        
                            self.openSweep_() ;
                        end
                    end
                else
                    % This means we should end the run

                    % At this point, self.IsPerformingRun_ is true
                    % Do post-run clean up
                    if self.AreAllSweepsCompleted_ ,
                        % Wrap up run in which all sweeps completed
                        self.wrapUpRunInWhichAllSweepsCompleted_() ;
                    else
                        % Something went wrong
                        self.abortOngoingRun_();
                    end
                    % At this point, self.IsPerformingRun_ is false
                end
            catch exception
                if self.IsPerformingRun_ ,
                    if self.IsPerformingSweep_ ,
                        self.abortTheOngoingSweep_() ;
                    end
                    self.abortOngoingRun_() ;
                end
                if self.IsHeaded_ ,
                    self.broadcast('RaiseDialogOnException', exception);
                else
                    throw(exception) ;
                end
            end                
        end  % handleTimerTick()

%         function handleTimerError(self, event)
% %             %event  %#ok<NOPRT>
% %             %event.Data
% %             %fprintf('The timer had an error\n') ;
% %             if self.IsPerformingRun_ ,
% %                 if self.IsPerformingSweep_ ,
% %                     self.abortTheOngoingSweep_() ;
% %                 end
% %                 self.abortOngoingRun_() ;
% %             end
% %             exception = MException(event.Data.messageID, event.Data.message) ;
% %             if self.IsHeaded_ ,
% %                 self.broadcast('RaiseDialogOnException', exception);
% %             else
% %                 throw(exception) ;
% %             end
%         end  % handleTimerError()
        
%         function closeSweep_(self)
%             % End the sweep in the way appropriate, either a "complete",
%             % a stop, or an abort.
%             % Precondition: self.IsPerformingSweep_ (is true)
%             if self.DidLooperCompleteSweep_ ,
%                 self.completeTheOngoingSweep_() ;
%                 %didCompleteSweep = true ;
%             elseif self.WasRunStopped_ ,
%                 self.stopTheOngoingSweep_() ;
%                 %didCompleteSweep = false ;
%             else
%                 % Something must have gone wrong...
%                 self.abortTheOngoingSweep_() ;
%                 %didCompleteSweep = false ;
%             end
% 
%             % No errors thrown, so set return values accordingly
%             %didThrow = false ;  
%             %exception = [] ;
%         end  % function
        
%         function performSweep_(self)
%             % time between subsequent calls to this, etc
%             t=toc(self.FromRunStartTicId_);
%             self.TimeOfLastWillPerformSweep_=t;
%             
%             % clear timing stuff that is strictly within-sweep
%             %self.TimeOfLastSamplesAcquired_=[];            
%             
%             % Reset the sample count for the sweep
%             self.NScansAcquiredSoFarThisSweep_ = 0;
%             
%             % update the current time
%             self.t_=0;            
%             
%             % Pretend that we last polled at time 0
%             self.TimeOfLastPollInSweep_ = 0 ;  % s 
%             
%             % As of the next line, the wavesurferModel is offically
%             % performing a sweep
%             % Actually, we'll wait until a bit later to set this true
%             %self.IsPerformingSweep_ = true ;
%             
%             % Call user functions 
%             self.callUserMethod_('startingSweep');            
%             
%             % Call startingSweep() on all the enabled subsystems
%             for idx = 1:numel(self.Subsystems_)
%                 if self.Subsystems_{idx}.IsEnabled ,
%                     self.Subsystems_{idx}.startingSweep();
%                 end
%             end
% 
%             % Wait for the looper to respond
%             self.LooperIPCRequester_.send('startingSweep', self.NSweepsCompletedInThisRun_+1) ;
%             timeout = 12 ;  % s
%             err = self.LooperIPCRequester_.waitForResponse(timeout, 'startingSweep') ;
%             if ~isempty(err) ,
%                 % Something went wrong
%                 self.abortOngoingRun_();
%                 self.changeReadiness_(+1);
%                 throw(err);
%             end
% 
%             % Set the sweep timer
%             self.FromSweepStartTicId_=tic();
% 
%             % Pulse the master trigger to start the sweep!
%             %fprintf('About to pulse the master trigger!\n');
%             self.IsPerformingSweep_ = true ;  
%                 % have to wait until now to set this true, so that old
%                 % samplesAcquired messages are ignored when we're
%                 % waiting for looperReadyForSweep and
%                 % refillerReadyForSweep messages above.
%             self.Triggering_.pulseBuiltinTrigger();
% 
%             % Now poll
%             ~self.DidAnySweepFailToCompleteSoFar_ = false ;
%             %self.DidLooperCompleteSweep_ = false ;
%             %self.DidRefillerCompleteSweep_ = ~self.Stimulation_.IsEnabled ;  % if stim subsystem is disabled, then the refiller is automatically done
%             self.WasRunStoppedInLooper_ = false ;
%             self.WasRunStoppedInRefiller_ = false ;
%             self.WasRunStopped_ = false ;
%             %self.TimeOfLastDrawnow_ = toc(drawnowTicId) ;  % don't really do a drawnow() now, but that's OK            
%             %profile resume ;
%             while ~(~self.DidAnySweepFailToCompleteSoFar_ || self.WasRunStopped_) ,
%                 %fprintf('At top of within-sweep loop...\n') ;
%                 self.LooperIPCSubscriber_.processMessagesIfAvailable() ;  % process all available messages, to make sure we keep up
%                 self.RefillerIPCSubscriber_.processMessagesIfAvailable() ;  % process all available messages, to make sure we keep up
%                 % do a drawnow() if it's been too long...
%                 timeSinceLastDrawNow = toc(self.DrawnowTicId_) - self.TimeOfLastDrawnow_ ;
%                 if timeSinceLastDrawNow > 0.1 ,  % 0.1 s, hence 10 Hz
%                     drawnow() ;
%                     self.TimeOfLastDrawnow_ = toc(self.DrawnowTicId_) ;
%                 end
%             end    
%             %profile off ;
%             %isSweepComplete = ~self.DidAnySweepFailToCompleteSoFar_ ;  % don't want to rely on this state more than we have to
%             %didUserRequestStop = self.WasRunStopped_ ;  
%             %~self.DidAnySweepFailToCompleteSoFar_ = [] ;
%             %self.WasRunStoppedInLooper_ = [] ;
%             %self.WasRunStoppedInRefiller_ = [] ;
%             %self.WasRunStopped_ = [] ;  % want this to stay as true/false
%             %end  % function                
% 
%             % End the sweep in the way appropriate, either a "complete",
%             % a stop, or an abort
%             if ~self.DidAnySweepFailToCompleteSoFar_ ,
%                 self.completeTheOngoingSweep_() ;
%                 %didCompleteSweep = true ;
%             elseif self.WasRunStopped_ ,
%                 self.stopTheOngoingSweep_();
%                 %didCompleteSweep = false ;
%             else
%                 % Something must have gone wrong...
%                 self.abortTheOngoingSweep_();
%                 %didCompleteSweep = false ;
%             end
% 
%             % No errors thrown, so set return values accordingly
%             %didThrow = false ;  
%             %exception = [] ;
%         end  % function

        % A run can end for one of three reasons: it is manually stopped by
        % the user in the middle, something goes wrong, and it completes
        % normally.  We call the first a "stop", the second an "abort", and
        % the third a "complete".  This terminology is (hopefully) used
        % consistently for methods like cleanUpAfterCompletedSweep_,
        % abortTheOngoingSweep_, etc.  Also note that this sense of
        % "abort" is not the same as it is used for DAQmx tasks.  For DAQmx
        % tasks, it means, roughly: stop the task, uncommit it, and
        % unreserve it.  Sometimes a Waversurfer "abort" leads to a DAQmx "abort"
        % for one or more DAQmx tasks, but this is not necessarily the
        % case.
        
%         function cleanUpAfterCompletedSweep_(self)
%             % Clean up after a sweep completes successfully.
%             
%             %fprintf('WavesurferModel::cleanUpAfterSweep_()\n');
%             %dbstack
%             
%             % Notify all the subsystems that the sweep is done
%             for idx = 1: numel(self.Subsystems_)
%                 if self.Subsystems_{idx}.IsEnabled
%                     self.Subsystems_{idx}.completingSweep();
%                 end
%             end
%             
%             % Bump the number of completed sweeps
%             self.NSweepsCompletedInThisRun_ = self.NSweepsCompletedInThisRun_ + 1;
%         
%             % Broadcast event
%             self.broadcast('DidCompleteSweep');
%             
%             % Call user method
%             self.callUserMethod_('didCompleteSweep');
%         end  % function
                
        function completeTheOngoingSweep_(self)
            %fprintf('WavesurferModel::completeTheOngoingSweep_()\n') ;
            
            self.dataAvailable_() ;  % Process any remaining data in the samples buffer

            % Bump the number of completed sweeps
            self.NSweepsCompletedInThisRun_ = self.NSweepsCompletedInThisRun_ + 1;

            % Notify all the necessary subsystems that the sweep is done
            %for idx = 1: numel(self.Subsystems_)
            %    if self.Subsystems_{idx}.IsEnabled
            %        self.Subsystems_{idx}.completingSweep();
            %    end
            %end
            if self.Logging_.IsEnabled ,
                self.Logging_.completingSweep() ;
            end                           

            % Call user method
            self.callUserMethod_('completingSweep');

            % As of next line, sweep is officially over
            %self.NSweepsPerformed_ = self.NSweepsPerformed_ + 1 ;
            %self.IsMostRecentSweepComplete_ = true ;
            self.AreAllSweepsCompleted_ = (self.NSweepsCompletedInThisRun_>=self.NSweepsPerRun) ;
            self.IsPerformingSweep_ = false ;

            % Broadcast event
            self.broadcast('DidCompleteSweep');
        end
        
        function stopTheOngoingSweep_(self)
            %fprintf('WavesurferModel::stopTheOngoingSweep_()\n') ;

            % Stop the ongoing sweep following user request            
            self.dataAvailable_() ;  % Process any remaining data in the samples buffer

            % Notify all the needed subsystems that the sweep was stopped
            % (in reverse order if multiple)
            if self.Logging_.IsEnabled , 
                self.Logging_.stoppingSweep() ;
            end
            
            % Call user method
            self.callUserMethod_('stoppingSweep');                    

            % As of next line, sweep is officially over
            %self.NSweepsPerformed_ = self.NSweepsPerformed_ + 1 ;
            self.DidAnySweepFailToCompleteSoFar_ = true ;
            self.IsPerformingSweep_ = false ;                    
        end  % function

        function abortTheOngoingSweep_(self)
            % Clean up after a sweep shits the bed.
            %fprintf('WavesurferModel::abortTheOngoingSweep_()\n') ;
            
            % Notify all needed subsystems that the sweep aborted, in
            % reverse order
%             for i = numel(self.Subsystems_):-1:1 ,
%                 if self.Subsystems_{i}.IsEnabled ,
%                     try 
%                         self.Subsystems_{i}.abortingSweep();
%                     catch me
%                         % In theory, Subsystem::abortingSweep() never
%                         % throws an exception
%                         % But just in case, we catch it here and ignore it
%                         disp(me.getReport());
%                     end
%                 end
%             end
            if self.Logging_.IsEnabled ,
                try
                    self.Logging_.abortingSweep();
                catch me
                    % In theory, ws.Logging::abortingSweep() never
                    % throws an exception
                    % But just in case, we catch it here and ignore it
                    disp(me.getReport());
                end
            end
            
            % Call user method
            self.callUserMethod_('abortingSweep');
            
            % declare the sweep over
            %self.NSweepsPerformed_ = self.NSweepsPerformed_ + 1 ;
            self.DidAnySweepFailToCompleteSoFar_ = true ;
            self.IsPerformingSweep_ = false ;            
        end  % function
                
        function wrapUpRunInWhichAllSweepsCompleted_(self)
            % Clean up after all sweeps complete successfully.
            
            % Stop, delete the timer
            if ~isempty(self.TheBigTimer_)
                if isvalid(self.TheBigTimer_) ,
                    stop(self.TheBigTimer_) ;
                    delete(self.TheBigTimer_) ;
                end
                self.TheBigTimer_ = [] ;
            end                
                                    
            % Notify other processes
            %self.IPCPublisher_.send('completingRun') ;
            self.Looper_.completingRun() ;
            didStopEpisode = self.Refiller_.completingRun() ;
            if didStopEpisode ,
                self.callUserMethod_('stoppingEpisode');
            end            
            
            % Notify subsystems that need to know, in forward order (why
            % forward?)
%             for idx = 1: numel(self.Subsystems_) ,
%                 if self.Subsystems_{idx}.IsEnabled ,
%                     self.Subsystems_{idx}.completingRun() ;
%                 end
%             end
            if self.Display_.IsEnabled ,
                self.Display_.completingRun() ;
            end
            self.Triggering_.completingRun() ;
            if self.Logging_.IsEnabled ,
                self.Logging_.completingRun() ;
            end
            
            % Call user method
            self.callUserMethod_('completingRun') ;

            % Finalize
            self.IsPerformingRun_ =  false ;
            self.setState_('idle') ;            
            self.NRunsCompleted_ = self.NRunsCompleted_ + 1 ;
            
            % Notify SI, if called for.  If there's a problem, log a
            % warning and then proceed.
            try
                %self.notifyScanImageThatRunCompletedNormallyIfYoked_() ;
            catch exception
                self.logWarning('ws:errorNotifyingScanImageThatRunFinishedNormally', ...
                                'The run finished normally, but there was a problem notifying ScanImage of this', ...
                                exception) ;
            end
            
        end  % function
        
        function stopTheOngoingRun_(self)
            % Called when the user stops the run in the middle, typically
            % by pressing the stop button.            
            
            % Stop, delete the timer
            if ~isempty(self.TheBigTimer_)
                if isvalid(self.TheBigTimer_) ,
                    stop(self.TheBigTimer_) ;
                    delete(self.TheBigTimer_) ;
                end
                self.TheBigTimer_ = [] ;
            end                
                                    
            %fprintf('WavesurferModel::stopTheOngoingRun_()\n') ;
            
            % Notify other processes --- or not, we don't currently need
            % this
            %self.IPCPublisher_.send('didStopRun') ;

            % No need to notify other processes, already did this by
            % sending 'frontendWantsToStopRun' message
            % % Notify other processes
            % self.IPCPublisher_.send('stoppingRun') ;

            % Notify subsystems, in reverse of starting order
%             for idx = numel(self.Subsystems_):-1:1 ,
%                 if self.Subsystems_{idx}.IsEnabled ,
%                     self.Subsystems_{idx}.stoppingRun() ;
%                 end
%             end
            if self.Logging_.IsEnabled ,
                self.Logging_.stoppingRun(self) ;
            end
            self.Triggering_.stoppingRun() ;
            if self.Display_.IsEnabled ,
                self.Display_.stoppingRun() ;
            end
            
            % Call user method
            self.callUserMethod_('stoppingRun');                
            
            % Set state back to idle
            self.IsPerformingRun_ =  false;
            self.setState_('idle');
        end  % function

        function abortOngoingRun_(self)
            % Called when a run fails for some undesirable reason.
            
            %fprintf('WavesurferModel::abortOngoingRun_()\n') ;

            % Stop, delete the timer
            if ~isempty(self.TheBigTimer_)
                if isvalid(self.TheBigTimer_) ,
                    stop(self.TheBigTimer_) ;
                    delete(self.TheBigTimer_) ;
                end
                self.TheBigTimer_ = [] ;
            end                
                        
            % Notify other processes
            %self.IPCPublisher_.send('abortingRun') ;
            self.Looper_.abortingRun() ;
            didAbortEpisode = self.Refiller_.abortingRun() ;
            if didAbortEpisode ,
                self.callUserMethod_('abortingEpisode');
            end                                
            
            % Notify subsystems, in reverse of starting order
%             for idx = numel(self.Subsystems_):-1:1 ,
%                 if self.Subsystems_{idx}.IsEnabled ,
%                     self.Subsystems_{idx}.abortingRun() ;
%                 end
%             end
            if self.Logging_.IsEnabled ,
                self.Logging_.abortingRun(self) ;
            end
            self.Triggering_.abortingRun() ;
            if self.Display_.IsEnabled ,
                self.Display_.abortingRun() ;
            end
            
            % Call user method
            self.callUserMethod_('abortingRun');
            
            % Set state back to idle
            self.IsPerformingRun_ = false ; 
            self.setState_('idle');
            
            % Notify SI if yoked
            try
                self.notifyScanImageThatRunAbortedIfYoked_() ;
            catch exception
                self.logWarning('ws:errorNotifyingScanImageThatRunAborted', ...
                                'The run aborted, and then there was a problem notifying ScanImage of this', ...
                                exception) ;
            end
            
        end  % function
        
        function samplesAcquired_(self, rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData)
            % Record the time
            %self.TimeInSweep_ = toc(self.FromSweepStartTicId_) ;
            
            % % Get the number of scans
            % nScans = size(rawAnalogData,1) ;

            % % Update the number of scans acquired
            % self.NScansAcquiredSoFarThisSweep_ = self.NScansAcquiredSoFarThisSweep_ + nScans ;

            % Add the new data to the storage buffer
            self.SamplesBuffer_.store(rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData) ;
            
            nScansInBuffer = self.SamplesBuffer_.nScansInBuffer() ;
            nScansPerUpdate = self.NScansPerUpdate_ ;
            %fprintf('nScansInBuffer, nScansPerUpdate: %10d, %10d\n',nScansInBuffer,nScansPerUpdate);
            if nScansInBuffer >= nScansPerUpdate ,
                %profile resume
                ticId = tic() ;
                self.dataAvailable_() ;
                durationOfDataAvailableCall = toc(ticId) ;
                % Update NScansPerUpdate to make sure we don't fall behind,
                % but must update no slower than 10x slower than desired.
                % (The buffer size is set to 100x
                % self.DesiredNScansPerUpdate_, and it's important that the
                % buffer be larger than the largest possible
                % nScansPerUpdate)
                fs = self.AcquisitionSampleRate ;
                nScansPerUpdateNew = min(10*self.DesiredNScansPerUpdate_ , ...
                                         max(2*round(durationOfDataAvailableCall*fs), ...
                                             self.DesiredNScansPerUpdate_ ) ) ;                                      
                self.NScansPerUpdate_= nScansPerUpdateNew ;
                %profile off
            end
        end
        
        function dataAvailable_(self)
            % The central method for handling incoming data.
            % Calls the dataAvailable() method on all the relevant subsystems, which handle display, logging, etc.            
            
            %fprintf('At top of WavesurferModel::dataAvailable_()\n') ;
            %tHere=tic();
            self.NTimesDataAvailableCalledSinceRunStart_ = self.NTimesDataAvailableCalledSinceRunStart_ + 1 ;
            [rawAnalogData,rawDigitalData,timeSinceRunStartAtStartOfData] = self.SamplesBuffer_.empty() ;            
            nScans = size(rawAnalogData,1) ;

            % Scale the new data, notify subsystems that we have new data
            if (nScans>0)
                % update the current time
                dt=1/self.AcquisitionSampleRate;
                self.t_=self.t_+nScans*dt;  % Note that this is the time stamp of the sample just past the most-recent sample

                % Scale the analog data
                channelScales=self.AIChannelScales(self.IsAIChannelActive);
                scalingCoefficients = self.Acquisition_.AnalogScalingCoefficients ;
                scaledAnalogData = ws.scaledDoubleAnalogDataFromRawMex(rawAnalogData, channelScales, scalingCoefficients) ;                
                %scaledAnalogData = ws.scaledDoubleAnalogDataFromRaw(rawAnalogData, channelScales) ;                
%                 inverseChannelScales=1./channelScales;  % if some channel scales are zero, this will lead to nans and/or infs
%                 if isempty(rawAnalogData) ,
%                     scaledAnalogData=zeros(size(rawAnalogData));
%                 else
%                     data = double(rawAnalogData);
%                     combinedScaleFactors = 3.0517578125e-4 * inverseChannelScales;  % counts-> volts at AI, 3.0517578125e-4 == 10/2^(16-1)
%                     scaledAnalogData=bsxfun(@times,data,combinedScaleFactors); 
%                 end
                
                % Store the data in the user cache
                self.Acquisition_.addDataToCache(rawAnalogData, rawDigitalData, self.AreSweepsFiniteDuration_) ;
                
                % 

                % Notify each relevant subsystem that data has just been acquired
                isSweepBased = self.AreSweepsFiniteDuration_ ;
                t = self.t_;
                if self.Logging_.IsEnabled ,
                    nActiveAnalogChannels = self.Acquisition_.NActiveAnalogChannels ;
                    nActiveDigitalChannels = self.Acquisition_.NActiveDigitalChannels ;                    
                    expectedSweepScanCount = self.ExpectedSweepScanCount ;
                    self.Logging_.dataAvailable(isSweepBased, ...
                                               t, ...
                                               scaledAnalogData, ...
                                               rawAnalogData, ...
                                               rawDigitalData, ...
                                               timeSinceRunStartAtStartOfData, ...
                                               nActiveAnalogChannels, ...
                                               nActiveDigitalChannels, ...
                                               expectedSweepScanCount);
                end
                if self.Display_.IsEnabled ,
                    [doesNeedClear, doesNeedDidSetXOffset] = ...
                        self.Display_.dataAvailable(isSweepBased, ...
                                                    t, ...
                                                    scaledAnalogData, ...
                                                    rawAnalogData, ...
                                                    rawDigitalData, ...
                                                    timeSinceRunStartAtStartOfData, ...
                                                    self.XSpan);
                    if doesNeedDidSetXOffset ,
                        self.broadcast('DidSetXOffset') ;
                    end
                    if doesNeedClear ,
                        self.syncTraces_() ;
                        self.broadcast('UpdateTraces') ;  
                    else
                        %self.broadcast('UpdateAfterDataAdded', t, scaledAnalogData, rawDigitalData) ;                        
                        self.addData_(t, scaledAnalogData, rawDigitalData) ;  % This will broadcast
                    end
                end
                self.callUserMethod_('dataAvailable');

                % Fire an event to cause views to sync
                self.broadcast('UpdateForNewData');
                
                % Do a drawnow(), to make sure user sees the changes, and
                % to process any button presses, etc.
                %drawnow();
            end
            %toc(tHere)
        end  % function
        
    end % protected methods block
        
    methods 
        % Allows access to protected and protected variables from ws.Encodable.
        function out = getPropertyValue_(self, name)
            out = self.(name) ;
        end  % function
        
        % Allows access to protected and protected variables from ws.Encodable.
        function setPropertyValue_(self, name, value)
            %if isequal(name,'IsDIChannelTerminalOvercommitted_')
            %    dbstack
            %    dbstop
            %end
            self.(name) = value ;
        end  % function
    end  % methods ( Access = protected )
    
    methods
        function addStarterChannelsAndStimulusLibrary(self)
            % Adds an AI channel, an AO channel, creates a stimulus and a
            % map in the stim library, and sets the current outputable to
            % the newly-created map.  This is intended to be run on a
            % "virgin" wavesurferModel.
        
            if self.DoesProtocolNeedSave || self.NAIChannels>0 || self.NAOChannels>0 || ~self.isStimulusLibraryEmpty() ,
                % do nothing
            else
                self.disableBroadcasts() ;
                self.addAIChannel() ;
                aoChannelIndex = self.addAOChannel() ;
                if ~isempty(aoChannelIndex) ,
                    aoChannelName = self.AOChannelNames{aoChannelIndex} ;
                    self.Stimulation_.setStimulusLibraryToSimpleLibraryWithUnitPulse({aoChannelName}) ;
                end
                self.IsDisplayEnabled = true ;
                self.DoesProtocolNeedSave_ = false ;  % this is a special case
                self.enableBroadcastsMaybe() ;            
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
            end
        end
    end  % methods block
    
    methods
        function executeIncomingCommand(self, command)
            % Executes a single received command.  There may be several
            % commands in a single command file.
            commandName = command.name ;
            parameters = command.parameters ;
            switch commandName ,
                case 'connect'
                    self.setIsYokedToScanImage_(true) ;
                case 'disconnect'
                    self.setIsYokedToScanImage_(false) ;
                case 'set-index-of-first-sweep-in-run'
                    self.Logging_.NextSweepIndex = str2double(parameters{1}) ;
                case 'set-number-of-sweeps-in-run'
                    self.NSweepsPerRun = str2double(parameters{1}) ;
                case 'set-data-file-folder-path'
                    self.DataFileLocation = parameters{1} ;
                case 'set-data-file-base-name'
                    self.DataFileBaseName = parameters{1} ;
                case 'set-data-file-session-index'
                    self.Logging_.SessionIndex = str2double(parameters{1}) ;
                case 'set-is-session-number-included-in-data-file-name'
                    self.Logging_.DoIncludeSessionIndex = logical(str2double(parameters{1})) ;
                case 'set-is-date-included-in-data-file-name'
                    self.Logging_.DoIncludeDate = logical(str2double(parameters{1})) ;                    
                case 'saving-configuration-file-at-full-path'
                    configurationFileName = parameters{1} ;
                    protocolFileName = ws.replaceFileExtension(configurationFileName, '.wsp') ;
                    self.saveProtocolFileGivenFileName(protocolFileName) ;
                case 'loading-configuration-file-at-full-path'
                    configurationFileName = parameters{1} ;
                    protocolFileName = ws.replaceFileExtension(configurationFileName, '.wsp') ;
                    self.openProtocolFileGivenFileName(protocolFileName) ;
                case 'saving-user-file-at-full-path'
                    siUserFileName = parameters{1} ;
                    [~, putativeProfileName, ~] = fileparts(siUserFileName) ;
                    if isequal(putativeProfileName, self.CurrentProfileName) ,
                        % Write the preferences to disk, just for the hell of it
                        self.savePreferences_(self.CurrentProfileName) ;
                    else
                        if ismember(putativeProfileName, self.ProfileNames) ,
                            newProfileName = putativeProfileName ;
                            self.CurrentProfileName = newProfileName ;  % this will update view                            
                        else
                            newProfileName = putativeProfileName ;
                            temporaryProfileName = self.createNewProfile() ;  % this will update view
                            self.CurrentProfileName = temporaryProfileName ;  % this will update view
                            self.renameCurrentProfile(newProfileName) ;  % this will update view          
                            self.savePreferences_(self.CurrentProfileName) ;  % what the hell                                                        
                        end
                    end                        
                case 'loading-user-file-at-full-path'
                    siUserFileName = parameters{1} ;
                    [~, putativeProfileName, ~] = fileparts(siUserFileName) ;
                    if isequal(putativeProfileName, self.CurrentProfileName) ,
                        % Do nothing
                    else
                        if ismember(putativeProfileName, self.ProfileNames) ,
                            newProfileName = putativeProfileName ;
                            self.CurrentProfileName = newProfileName ;  % this will update view                            
                        else
                            error('Can''t find profile named "%s"', putativeProfileName) ;
                        end
                    end                        
                case 'record'
                    self.record()
                case 'play'
                    self.play()
                case 'stop'
                    % This is similar to the user pressing the stop button
                    % in the WS UI
                    self.stop();
                case 'did-complete-acquisition-mode-normally'
                    % SI sends this when "acquisiition mode" (aka a loop or
                    % grab) finishes without problems on its end.
                    % Currently, we don't do anything in response to this.
                case 'error'
                    % SI send this when something goes wrong during
                    % acquisition mode.  This is equivalent to an "abort"
                    % in WS terminology.
                    if isempty(parameters) ,
                        message = '' ;  %#ok<NASGU>
                    else
                        message = parameters{1} ;                     %#ok<NASGU>
                    end
                    self.stop() ;
                case 'ping'
                    % do nothing (exception send an "OK" back to SI
                otherwise
                    error('WavesurferModel:UnknownScanImageCommand', ...
                          'Received unknown command ''%s'' from ScanImage', commandName) ;
            end
        end % function                
    end  % public methods block
    
    methods (Access=protected)            
        function commandScanImageToStartLoopIfYoked_(self)
            % Sends acquisition parameters to ScanImage and start
            % acquisition
            
            nAcqsInSet=self.NSweepsPerRun;
            iFirstAcqInSet=self.Logging_.NextSweepIndex;
            
            %fprintf(fid,'Internally generated: %d\n',);
            %fprintf(fid,'InputPFI| %d\n',pfiLineId);
            %fprintf(fid,'Edge type| %s\n',edgeTypeString);
            commandFileAsString = ...
                sprintf('%s\n%s|%d\n%s|%d\n%s|%d\n%s|%s\n%s|%s\n%s\n', ...
                        '6', ...
                        'set-log-file-counter', iFirstAcqInSet, ...
                        'set-acq-count-in-loop', nAcqsInSet, ...
                        'set-log-enabled', self.Logging_.IsEnabled, ...
                        'set-log-file-folder-path', fileparts(self.NextRunAbsoluteFileName), ...
                        'set-log-file-base-name', self.Logging_.AugmentedBaseName, ...
                        'loop') ;
            
            self.CommandClient_.sendCommandFileAsString(commandFileAsString);
        end  % function        
        
        function notifyScanImageThatRunCompletedNormallyIfYoked_(self)
            commandFileAsString = sprintf('1\ndid-complete-run-normally\n') ;  % sprintf converts '\n' to a newline
            self.CommandClient_.sendCommandFileAsString(commandFileAsString);
        end
        
        function notifyScanImageThatRunAbortedIfYoked_(self)
            commandFileAsString = sprintf('1\nerror\n') ;  % sprintf converts '\n' to a newline
            self.CommandClient_.sendCommandFileAsString(commandFileAsString);
        end
        
        function notifyScanImageThatWavesurferIsQuittingIfYoked_(self)
            commandFileAsString = sprintf('1\nwavesurfer-is-quitting\n') ;  % sprintf converts '\n' to a newline
            commandClient = self.CommandClient_ ;
            if ~isempty(commandClient) && isvalid(commandClient) ,
                commandClient.sendCommandFileAsString(commandFileAsString) ;
            end
        end
        
        function notifyScanImageThatSavingProtocolFileIfYoked_(self, absoluteProtocolFileName)
            commandFileAsString = sprintf('1\nsaving-protocol-file-at-full-path| %s\n', absoluteProtocolFileName) ;
            self.CommandClient_.sendCommandFileAsString(commandFileAsString) ;
        end  % function
        
        function notifyScanImageThatOpeningProtocolFileIfYoked_(self, absoluteProtocolFileName)
            commandFileAsString = sprintf('1\nopening-protocol-file-at-full-path| %s\n', absoluteProtocolFileName) ;
            self.CommandClient_.sendCommandFileAsString(commandFileAsString) ;
        end  % function
        
        function notifyScanImageThatSavingPreferencesIfYoked_(self, profileName)
            absolutePreferencesFilePath = ws.preferencesFileNameFromProfileName(profileName) ;
            commandFileAsString = sprintf('1\nsaving-user-file-at-full-path| %s\n', absolutePreferencesFilePath) ;
            self.CommandClient_.sendCommandFileAsString(commandFileAsString) ;
        end  % function

        function notifyScanImageThatLoadingPreferencesIfYoked_(self, profileName)
            absolutePreferencesFilePath = ws.preferencesFileNameFromProfileName(profileName) ;
            commandFileAsString = sprintf('1\nopening-user-file-at-full-path| %s\n', absolutePreferencesFilePath) ;
            self.CommandClient_.sendCommandFileAsString(commandFileAsString) ;
        end  % function
    end % methods
    
%     properties (Hidden, SetAccess=protected)
%         mdlPropAttributes = struct();   % ws.WavesurferModel.propertyAttributes();        
%         mdlHeaderExcludeProps = {};
%     end
    
%     methods (Static)
%         function s = propertyAttributes()
%             s = struct();
%             
%             %s.NSweepsPerRun = struct('Attributes',{{'positive' 'integer' 'finite' 'scalar' '>=' 1}});
%             %s.SweepDuration = struct('Attributes',{{'positive' 'finite' 'scalar'}});
%             s.AreSweepsFiniteDuration = struct('Classes','binarylogical');  % dependency on AreSweepsContinuous handled in the setter
%             s.AreSweepsContinuous = struct('Classes','binarylogical');  % dependency on IsTrailBased handled in the setter
%         end  % function
%     end  % class methods block
    
%     methods
%         function nScopesMayHaveChanged(self)
%             self.broadcast('NScopesMayHaveChanged');
%         end
%     end
    
    methods
        function value = get.NTimesDataAvailableCalledSinceRunStart(self)
            value=self.NTimesDataAvailableCalledSinceRunStart_;
        end
    end

    methods
        function value = get.ClockAtRunStart(self)
            value = self.ClockAtRunStart_ ;
        end
    end
    
    methods
        function openProtocolFileGivenFileName(self, fileName)
            % Actually loads the named protocol file.  fileName should be a
            % file name referring to a file that is known to be
            % present, at least as of a few milliseconds ago.
            if ~self.isIdleSensuLato() || self.DoesProtocolNeedSave ,
                return
            end
            self.changeReadiness_(-1);
            if ws.isFileNameAbsolute(fileName) ,
                absoluteFileName = fileName ;
            else
                absoluteFileName = fullfile(pwd(),fileName) ;
            end
            saveStruct = load('-mat',absoluteFileName) ;
            wavesurferModelSettings = saveStruct.ws_WavesurferModel ;
            newModel = ws.decodeEncodingContainer(wavesurferModelSettings, self) ;
            self.broadcast('RequestLayoutForAllWindows');  % Have to prompt the figure/controller to tell us this
              % We do this in case something in the protocol file is weird, so that we can
              % fallback on the current positions and visibilities.
            self.mimicProtocolThatWasJustLoaded_(newModel) ;
            %if isfield(saveStruct, 'layoutForAllWindows') ,
            %    self.LayoutForAllWindows_ = saveStruct.layoutForAllWindows ;
            %end
            %self.Acquisition_.invalidateDataCache() ;
            self.AbsoluteProtocolFileName_ = absoluteFileName ;
            self.HasUserSpecifiedProtocolFileName_ = true ; 
            self.DoesProtocolNeedSave_ = false ;
            self.LastProtocolFilePath_ = absoluteFileName ;
            
            self.broadcast('Update') ;
%             self.broadcast('UpdateLogging') ;
%             self.broadcast('UpdateUserCodeManager') ;
%             self.broadcast('UpdateElectrodeManager') ;
%             self.broadcast('UpdateTestPulser') ;            
%             self.broadcast('UpdateTraces') ;
%             self.broadcast('UpdateDisplay') ;
%             self.broadcast('UpdateStimulusLibrary') ;
%             self.broadcast('UpdateStimulusPreview') ;
%             self.broadcast('UpdateTriggering') ;
%             self.broadcast('UpdateMain');
%             self.broadcast('UpdateGeneral') ;
%             self.broadcast('UpdateChannels') ;
            self.broadcast('LayoutAllWindows') ;
            
            self.callUserMethod_('wake');  % wake the user object
            
            %siConfigFilePath = ws.replaceFileExtension(absoluteFileName, '.cfg') ;
            self.notifyScanImageThatOpeningProtocolFileIfYoked_(absoluteFileName);
            self.changeReadiness_(+1);
        end  % function
    end
    
%     methods
%         function setLayoutForAllWindows_(self, layoutForAllWindows)
%             self.LayoutForAllWindows_ = layoutForAllWindows ;
%         end
%     end
    
    methods
%         function saveProtocolFileGivenAbsoluteFileNameAndWindowsLayout(self, absoluteFileName, layoutForAllWindows)
%             self.LayoutForAllWindows_ = layoutForAllWindows ;
%             self.saveProtocolFileGivenAbsoluteFileName(absoluteFileName) ;
%         end

        function saveProtocolFileGivenFileName(self, fileName)
            %wavesurferModelSettings=self.encodeConfigurablePropertiesForFileType('cfg');
            self.changeReadiness_(-1);       
            if ws.isFileNameAbsolute(fileName) ,
                absoluteFileName = fileName ;
            else
                absoluteFileName = fullfile(pwd(),fileName) ;
            end
            self.broadcast('RequestLayoutForAllWindows');  % Have to prompt the figure/controller to tell us this
              % If headless, window position and visibility fields will not change
            self.callUserMethod_('willSaveToProtocolFile');  % notify the user object we're about to save  
            wavesurferModelSettings = ws.encodeForPersistence(self) ;
            %wavesurferModelSettingsVariableName=self.getEncodedVariableName();
            wavesurferModelSettingsVariableName = 'ws_WavesurferModel' ;
            versionString = ws.versionString() ;
            %layoutForAllWindows = self.LayoutForAllWindows_ ;
%             saveStruct=struct(wavesurferModelSettingsVariableName,wavesurferModelSettings, ...
%                               'layoutForAllWindows',layoutForAllWindows, ...
%                               'versionString',versionString);  %#ok<NASGU>
            saveStruct=struct(wavesurferModelSettingsVariableName,wavesurferModelSettings, ...
                              'versionString',versionString);  %#ok<NASGU>
            save('-mat','-v7.3',absoluteFileName,'-struct','saveStruct');     
            self.AbsoluteProtocolFileName_ = absoluteFileName ;
            %self.broadcast('DidSetAbsoluteProtocolFileName');            
            self.HasUserSpecifiedProtocolFileName_ = true ;
            self.DoesProtocolNeedSave_ = false ;
            self.LastProtocolFilePath_ = absoluteFileName ;
            %siConfigFilePath = ws.replaceFileExtension(absoluteFileName, '.cfg') ;
            self.notifyScanImageThatSavingProtocolFileIfYoked_(absoluteFileName) ;
            self.changeReadiness_(+1);            
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end

        function pretendThatProtocolWasSaved_(self) 
            % Don't use this unless you know what you're doing --- it doesn't even do an
            % update
            self.DoesProtocolNeedSave_ = false ;
        end
        
%         function saveProtocolFileGivenAbsoluteFileName(self, absoluteFileName)
%             % This is here for backwards-compatibility
%             self.saveProtocolFileGivenFileName(absoluteFileName) ;
%         end
        
    end        
    
%     methods
%         function openUserFileGivenFileName(self, fileName)
%             % Actually opens the named user file.  fileName should be an
%             % file name referring to a file that is known to be
%             % present, at least as of a few milliseconds ago.
%             self.changeReadiness_(-1) ;
%             if ws.isFileNameAbsolute(fileName) ,
%                 absoluteFileName = fileName ;
%             else
%                 absoluteFileName = fullfile(pwd(),fileName) ;
%             end                        
%             saveStruct=load('-mat',absoluteFileName) ;
%             wavesurferModelSettingsVariableName = 'ws_WavesurferModel' ;            
%             wavesurferModelSettings=saveStruct.(wavesurferModelSettingsVariableName) ;
%             newModel = ws.decodeEncodingContainer(wavesurferModelSettings, self) ;
%             self.mimicUserSettings_(newModel) ;
%             self.AbsoluteUserSettingsFileName_ = absoluteFileName ;
%             self.HasUserSpecifiedUserSettingsFileName_ = true ;            
%             if self.DoUsePreferences ,
%                 ws.setPreference('LastUserFilePath', absoluteFileName) ;
%             end
%             self.notifyScanImageThatLoadingPreferencesIfYoked_(absoluteFileName) ;
%             self.changeReadiness_(+1) ;            
%             self.broadcast('UpdateFastProtocols') ;
%             self.broadcast('Update') ;
%         end
%     end        

%     methods
%         function saveUserFileGivenFileName(self, fileName)
%             self.changeReadiness_(-1) ;
%             if ws.isFileNameAbsolute(fileName) ,
%                 absoluteFileName = fileName ;
%             else
%                 absoluteFileName = fullfile(pwd(),fileName) ;
%             end                        
%             userSettings = ws.encodeForPersistence(self) ;
%             wavesurferModelSettingsVariableName = 'ws_WavesurferModel' ;
%             versionString = ws.versionString() ;
%             saveStruct=struct(wavesurferModelSettingsVariableName,userSettings, ...
%                               'versionString',versionString) ;  %#ok<NASGU>
%             save('-mat','-v7.3',absoluteFileName,'-struct','saveStruct') ;     
%             self.AbsoluteUserSettingsFileName_ = absoluteFileName ;
%             self.HasUserSpecifiedUserSettingsFileName_ = true ;            
%             if self.DoUsePreferences ,
%                 ws.setPreference('LastUserFilePath', absoluteFileName) ;
%             end
%             self.notifyScanImageThatSavingUserFileIfYoked_(absoluteFileName) ;
%             self.changeReadiness_(+1) ;            
%             self.broadcast('Update') ;            
%         end  % function
%     end
    
%     methods
%         function set.AbsoluteProtocolFileName(self,newValue)
%             % SetAccess is protected, no need for checks here
%             self.AbsoluteProtocolFileName=newValue;
%             self.broadcast('DidSetAbsoluteProtocolFileName');
%         end
%     end

%     methods
%         function set.AbsoluteUserSettingsFileName(self,newValue)
%             % SetAccess is protected, no need for checks here
%             self.AbsoluteUserSettingsFileName=newValue;
%             self.broadcast('DidSetAbsoluteUserSettingsFileName');
%         end
%     end

%     methods
%         function updateFastProtocol(self)
%             % Called by one of the child FastProtocol's when it is changed
%             self.broadcast('UpdateFastProtocols');
%             self.broadcast('Update');  % need to update main window also
%         end
%     end

%     methods
%         function result = allDigitalTerminalIDs(self)
%             nDigitalTerminalIDsInHardware = self.NDIOTerminals ;
%             result = 0:(nDigitalTerminalIDsInHardware-1) ;              
%         end
        
%         function result = digitalTerminalIDsInUse(self)
%             inputDigitalTerminalIDs = self.Acquisition_.DigitalTerminalIDs ;
%             outputDigitalTerminalIDs = self.Stimulation_.DigitalTerminalIDs ;
%             result = sort([inputDigitalTerminalIDs outputDigitalTerminalIDs]) ;
%         end
        
%         function result = freeDigitalTerminalIDs(self)
%             allIDs = self.allDigitalTerminalIDs() ;  
%             inUseIDs = self.digitalTerminalIDsInUse() ;
%             result = setdiff(allIDs, inUseIDs) ;
%         end
        
%         function result = isDigitalTerminalIDInUse(self, DigitalTerminalID)
%             inUseDigitalTerminalIDs = self.digitalTerminalIDsInUse() ;
%             result = ismember(DigitalTerminalID, inUseDigitalTerminalIDs) ;
%         end
%     end
    
    methods
        function setSingleDOChannelTerminalID(self, iChannel, terminalID)
            wasSet = self.Stimulation_.setSingleDigitalTerminalID_(iChannel, terminalID) ;            
            %self.didSetDigitalOutputTerminalID() ;
            self.syncIsDIOChannelTerminalOvercommitted_() ;
            self.broadcast('UpdateChannels') ;
            if wasSet ,
                %self.Parent.singleDigitalOutputTerminalIDWasSetInStimulationSubsystem(i) ;
                value = self.Stimulation_.DigitalTerminalIDs(iChannel) ;  % value is possibly normalized, terminalID is not
%                 self.IPCPublisher_.send('singleDigitalOutputTerminalIDWasSetInFrontend', ...
%                                         iChannel, value, self.IsDOChannelTerminalOvercommitted ) ;
                self.Looper_.singleDigitalOutputTerminalIDWasSetInFrontend(self.PrimaryDeviceName, ...
                                                                           self.IsPrimaryDeviceAPXIDevice, ...
                                                                           self.DOChannelTerminalIDs, ...
                                                                           self.IsDOChannelTimed, ...
                                                                           self.DOChannelStateIfUntimed, ...
                                                                           self.IsDOChannelTerminalOvercommitted) ;
                self.Refiller_.singleDigitalOutputTerminalIDWasSetInFrontend(iChannel, ...
                                                                             value, ...
                                                                             self.IsDOChannelTerminalOvercommitted ) ;
            end
        end
        
%         function singleDigitalOutputTerminalIDWasSetInStimulationSubsystem(self, i)
%             % This only gets called if the value was actually set.
%             value = self.Stimulation_.DigitalTerminalIDs(i) ;
%             keyboard
%             self.IPCPublisher_.send('singleDigitalOutputTerminalIDWasSetInFrontend', ...
%                                     i, value, self.IsDOChannelTerminalOvercommitted ) ;
%         end

%         function digitalOutputStateIfUntimedWasSetInStimulationSubsystem(self)
%             value = self.DOChannelStateIfUntimed ;
%             self.Looper_.digitalOutputStateIfUntimedWasSetInFrontend(value, self.IsDOChannelTimed) ;
%             self.Refiller_.digitalOutputStateIfUntimedWasSetInFrontend(value) ;
%         end
        
        function isDigitalChannelTimedWasSetInStimulationSubsystem(self)
            value = self.IsDOChannelTimed ;
            % Notify the refiller first, so that it can release all the DO
            % channels
            %self.IPCPublisher_.send('isDigitalOutputTimedWasSetInFrontend',value) ;
            self.Looper_.isDigitalOutputTimedWasSetInFrontend(self.PrimaryDeviceName, ...
                                                              self.IsPrimaryDeviceAPXIDevice, ...
                                                              self.DOChannelTerminalIDs, ...
                                                              self.IsDOChannelTimed, ...
                                                              self.DOChannelStateIfUntimed, ...
                                                              self.IsDOChannelTerminalOvercommitted) ;
            self.Refiller_.isDigitalOutputTimedWasSetInFrontend(value) ;            
        end
        
        function newChannelIndex = addAIChannel(self)
            nextFreeDeviceNameAndTerminalIDMaybe = self.nextFreeAITerminal() ;
            if isempty(nextFreeDeviceNameAndTerminalIDMaybe) ,
                % No free AI terminals
                newChannelIndex = [] ;
            else
                nextFreeDeviceNameAndTerminalID = nextFreeDeviceNameAndTerminalIDMaybe(1) ;
                newChannelIndex = self.Acquisition_.addAnalogChannel(nextFreeDeviceNameAndTerminalID.deviceName, ...
                                                                     nextFreeDeviceNameAndTerminalID.terminalID) ;
                %self.Acquisition_.invalidateDataCache() ;                                                  
                self.DoesProtocolNeedSave_ = true ;                                                  
                self.syncIsAIChannelTerminalOvercommitted_() ;
                self.Display_.didAddAnalogInputChannel() ;
                self.syncTraces_() ;
                self.broadcast('UpdateTraces') ;
                self.broadcast('UpdateMain') ;
                %self.Ephys_.didChangeNumberOfInputChannels();
                self.broadcast('EMDidChangeNumberOfInputChannels');
                self.broadcast('UpdateTestPulser');
                self.broadcast('UpdateChannels');  % causes channels figure to update
                self.broadcast('DidChangeNumberOfInputChannels');
            end
        end  % function
        
        function newChannelIndex = addAOChannel(self)
            nextFreeDeviceNameAndTerminalIDMaybe = self.nextFreeAOTerminal() ;
            if isempty(nextFreeDeviceNameAndTerminalIDMaybe) ,
                newChannelIndex = [] ;
            else            
                nextFreeDeviceNameAndTerminalID = nextFreeDeviceNameAndTerminalIDMaybe(1) ;
                newChannelIndex = self.Stimulation_.addAnalogChannel(nextFreeDeviceNameAndTerminalID.deviceName, ...
                                                                     nextFreeDeviceNameAndTerminalID.terminalID) ;
                self.DoesProtocolNeedSave_ = true ;                                                  
                self.syncIsAOChannelTerminalOvercommitted_() ;
                %self.Ephys_.didChangeNumberOfOutputChannels() ;                
                self.broadcast('EMDidChangeNumberOfOutputChannels');
                self.broadcast('UpdateTestPulser') ;                
                self.broadcast('UpdateChannels') ;  % causes channels figure to update
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
            end
        end
        
        function channelIndex = addDIChannel(self)
            nextFreeDeviceNameAndTerminalIDMaybe = self.nextFreeDIOTerminal() ;
            if isempty(nextFreeDeviceNameAndTerminalIDMaybe) ,
                channelIndex = [] ;
            else
                nextFreeDeviceNameAndTerminalID = nextFreeDeviceNameAndTerminalIDMaybe(1) ;                
                channelIndex = self.Acquisition_.addDigitalChannel(nextFreeDeviceNameAndTerminalID.deviceName, ...
                                                                   nextFreeDeviceNameAndTerminalID.terminalID) ;
                %self.Acquisition_.invalidateDataCache() ;
                self.DoesProtocolNeedSave_ = true ;                                                  
                self.syncIsDIOChannelTerminalOvercommitted_() ;
                self.Display_.didAddDigitalInputChannel() ;
                self.syncTraces_() ;
                self.broadcast('UpdateTraces') ;
                self.broadcast('UpdateMain') ;
                %self.Ephys_.didChangeNumberOfInputChannels() ;                
                self.broadcast('EMDidChangeNumberOfInputChannels');
                self.broadcast('UpdateTestPulser');
                self.broadcast('UpdateChannels') ;  % causes channels figure to update
                self.broadcast('DidChangeNumberOfInputChannels');  % causes scope controllers to be synched with scope models
                self.Looper_.didAddDigitalInputChannelInFrontend(self.PrimaryDeviceName, ...
                                                                 self.IsPrimaryDeviceAPXIDevice, ...
                                                                 self.DOChannelTerminalIDs, ...
                                                                 self.IsDOChannelTimed, ...
                                                                 self.DOChannelStateIfUntimed, ...
                                                                 self.IsDOChannelTerminalOvercommitted) ;
                self.Refiller_.didAddDigitalInputChannelInFrontend(self.IsDOChannelTerminalOvercommitted) ;
            end
        end
        
        function newChannelIndex = addDOChannel(self)
            nextFreeDeviceNameAndTerminalIDMaybe = self.nextFreeDIOTerminal() ;
            if isempty(nextFreeDeviceNameAndTerminalIDMaybe) ,
                newChannelIndex = [] ;
            else
                nextFreeDeviceNameAndTerminalID = nextFreeDeviceNameAndTerminalIDMaybe(1) ;                
                %freeTerminalIDs = self.freeDigitalTerminalIDs() ;
                %allDeviceNames = self.AllDeviceNames ;
                newChannelIndex = self.Stimulation_.addDigitalChannel(nextFreeDeviceNameAndTerminalID.deviceName, ...
                                                                      nextFreeDeviceNameAndTerminalID.terminalID) ;
                %self.Display_.didAddDigitalOutputChannel() ;
                self.DoesProtocolNeedSave_ = true ;                                                  
                self.syncIsDIOChannelTerminalOvercommitted_() ;
                %self.Stimulation_.notifyLibraryThatDidChangeNumberOfOutputChannels_() ;
                self.broadcast('UpdateStimulusLibrary');
                self.broadcast('UpdateStimulusPreview') ;
                %self.Ephys_.didChangeNumberOfOutputChannels();
                self.broadcast('EMDidChangeNumberOfOutputChannels');
                self.broadcast('UpdateTestPulser') ;                                
                self.broadcast('UpdateChannels');  % causes channels figure to update
                %self.broadcast('DidChangeNumberOfOutputChannels');  % causes scope controllers to be synched with scope models
%                 channelNameForEachDOChannel = self.Stimulation_.DigitalChannelNames ;
%                 terminalIDForEachDOChannel = self.Stimulation_.DigitalTerminalIDs ;
%                 isTimedForEachDOChannel = self.IsDOChannelTimed ;
%                 onDemandOutputForEachDOChannel = self.DOChannelStateIfUntimed ;
%                 isTerminalOvercommittedForEachDOChannel = self.IsDOChannelTerminalOvercommitted ;
%                 self.IPCPublisher_.send('didAddDigitalOutputChannelInFrontend', ...
%                                         channelNameForEachDOChannel, ...
%                                         terminalIDForEachDOChannel, ...
%                                         isTimedForEachDOChannel, ...
%                                         onDemandOutputForEachDOChannel, ...
%                                         isTerminalOvercommittedForEachDOChannel) ;
                self.Looper_.didAddDigitalOutputChannelInFrontend(self.PrimaryDeviceName, ...
                                                                  self.IsPrimaryDeviceAPXIDevice, ...
                                                                  self.DOChannelTerminalIDs, ...
                                                                  self.IsDOChannelTimed, ...
                                                                  self.DOChannelStateIfUntimed, ...
                                                                  self.IsDOChannelTerminalOvercommitted) ;
%                 self.Refiller_.didAddDigitalOutputChannelInFrontend(channelNameForEachDOChannel, ...
%                                                                     terminalIDForEachDOChannel, ...
%                                                                     isTimedForEachDOChannel, ...
%                                                                     onDemandOutputForEachDOChannel, ...
%                                                                     isTerminalOvercommittedForEachDOChannel) ;
            end
        end        
        
        function deleteMarkedAIChannels(self)
            wasDeleted = self.Acquisition_.deleteMarkedAnalogChannels() ;
            %self.Acquisition_.invalidateDataCache() ;
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsAIChannelTerminalOvercommitted_() ;            
            self.Display_.didDeleteAnalogInputChannels(wasDeleted) ;
            self.syncTraces_() ;
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateTraces') ;
            self.broadcast('EMDidChangeNumberOfInputChannels');
            self.broadcast('UpdateTestPulser');
            self.broadcast('UpdateChannels');  % causes channels figure to update
            self.broadcast('DidChangeNumberOfInputChannels');  
        end
        
        function deleteMarkedDIChannels(self)
            wasDeleted = self.Acquisition_.deleteMarkedDigitalChannels() ;
            %self.Acquisition_.invalidateDataCache() ;
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsDIOChannelTerminalOvercommitted_() ;
            self.Display_.didDeleteDigitalInputChannels(wasDeleted) ;
            self.syncTraces_() ;
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateTraces') ;
            %self.Ephys_.didChangeNumberOfInputChannels() ;
            self.broadcast('EMDidChangeNumberOfInputChannels');
            self.broadcast('UpdateTestPulser');
            self.broadcast('UpdateChannels') ;  % causes channels figure to update
            self.broadcast('DidChangeNumberOfInputChannels') ;  
            %self.broadcast('DidMaybeChangeProtocol') ;
%             self.IPCPublisher_.send('didDeleteDigitalInputChannelsInFrontend', ...
%                                     self.IsDOChannelTerminalOvercommitted) ;
            self.Looper_.didDeleteDigitalInputChannelsInFrontend(self.PrimaryDeviceName, ...
                                                                 self.IsPrimaryDeviceAPXIDevice, ...
                                                                 self.DOChannelTerminalIDs, ...
                                                                 self.IsDOChannelTimed, ...
                                                                 self.DOChannelStateIfUntimed, ...
                                                                 self.IsDOChannelTerminalOvercommitted) ;
            self.Refiller_.didDeleteDigitalInputChannelsInFrontend(self.IsDOChannelTerminalOvercommitted) ;
        end
        
        function deleteMarkedAOChannels(self)
            self.Stimulation_.deleteMarkedAnalogChannels() ;
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsAOChannelTerminalOvercommitted_() ;            
            %self.Display_.didRemoveAnalogOutputChannel(nameOfRemovedChannel) ;
            %self.Ephys_.didChangeNumberOfOutputChannels();
            self.broadcast('EMDidChangeNumberOfOutputChannels');
            self.broadcast('UpdateTestPulser') ;            
%             self.Stimulation_.notifyLibraryThatDidChangeNumberOfOutputChannels_();  
%               % we might be able to call this from within
%               % self.Stimulation_.deleteMarkedAnalogChannels, and that would
%               % generally be better, but I'm afraid of introducing new
%               % bugs...
            self.broadcast('UpdateStimulusLibrary');
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateChannels');  % causes channels figure to update
            %self.broadcast('DidChangeNumberOfOutputChannels');  % causes scope controllers to be synched with scope models
        end
        
        function deleteMarkedDOChannels(self)
            % Determine which to delete, which to keep
            isToBeDeleted = self.IsDOChannelMarkedForDeletion ;
            isKeeper = ~isToBeDeleted ;
            
            % Turn off any untimed DOs that are about to be deleted
            digitalOutputStateIfUntimed = self.DOChannelStateIfUntimed ;
            self.DOChannelStateIfUntimed = digitalOutputStateIfUntimed & isKeeper ;                        
            
            % Make the needed changed to the Stimulation subsystem
            self.Stimulation_.deleteMarkedDigitalChannels_(isToBeDeleted) ;
            
            % Do all the things that need doing after that
            self.DoesProtocolNeedSave_ = true ;
            self.syncIsDIOChannelTerminalOvercommitted_() ;
            self.broadcast('UpdateStimulusLibrary');
            self.broadcast('UpdateStimulusPreview') ;
            %self.Ephys_.didChangeNumberOfOutputChannels();
            self.broadcast('EMDidChangeNumberOfOutputChannels');            
            self.broadcast('UpdateTestPulser') ;            
            self.broadcast('UpdateChannels');  % causes channels figure to update
%             channelNameForEachDOChannel = self.Stimulation_.DigitalChannelNames ;
%             terminalIDForEachDOChannel = self.Stimulation_.DigitalTerminalIDs ;
%             isTimedForEachDOChannel = self.IsDOChannelTimed ;
%             onDemandOutputForEachDOChannel = self.DOChannelStateIfUntimed ;
%             isTerminalOvercommittedForEachDOChannel = self.IsDOChannelTerminalOvercommitted ;
%             self.IPCPublisher_.send('didRemoveDigitalOutputChannelsInFrontend', ...
%                                     channelNameForEachDOChannel, ...
%                                     terminalIDForEachDOChannel, ...
%                                     isTimedForEachDOChannel, ...
%                                     onDemandOutputForEachDOChannel, ...
%                                     isTerminalOvercommittedForEachDOChannel) ;
            self.Looper_.didRemoveDigitalOutputChannelsInFrontend(self.PrimaryDeviceName, ...
                                                                  self.IsPrimaryDeviceAPXIDevice, ...
                                                                  self.DOChannelTerminalIDs, ...
                                                                  self.IsDOChannelTimed, ...
                                                                  self.DOChannelStateIfUntimed, ...
                                                                  self.IsDOChannelTerminalOvercommitted) ;
%             self.Refiller_.didRemoveDigitalOutputChannelsInFrontend(channelNameForEachDOChannel, ...
%                                                                     terminalIDForEachDOChannel, ...
%                                                                     isTimedForEachDOChannel, ...
%                                                                     onDemandOutputForEachDOChannel, ...
%                                                                     isTerminalOvercommittedForEachDOChannel) ;
        end        
    end  % public methods block
    
%     methods
%         function triggeringSubsystemJustStartedFirstSweepInRun(self)
%             % Called by the triggering subsystem just after first sweep
%             % started.
%             
%             % This means we need to run the main polling loop.
%             self.runWithinSweepPollingLoop_();
%         end  % function
%     end
    
    methods
        function mimic(self, other)
            % Cause self to resemble other.
            
            % Disable broadcasts for speed
            self.disableBroadcasts();
            
            % Get the list of property names for this file type
            propertyNames = ws.listPropertiesForPersistence(self);
            
            % Set each property to the corresponding one
            for i = 1:length(propertyNames) ,
                thisPropertyName=propertyNames{i};
                if any(strcmp(thisPropertyName,{'Triggering_', 'Acquisition_', 'Stimulation_', 'Display_', 'Ephys_', 'UserCodeManager_'})) ,
                    %self.(thisPropertyName).mimic(other.(thisPropertyName)) ;
                    self.(thisPropertyName).mimic(other.getPropertyValue_(thisPropertyName)) ;                    
                else
                    if isprop(other,thisPropertyName) ,
                        source = other.getPropertyValue_(thisPropertyName) ;
                        self.setPropertyValue_(thisPropertyName, source) ;
                    end
                end
            end
            
            % Re-enable broadcasts
            self.enableBroadcastsMaybe();
            
            % Broadcast update
            self.broadcast('Update');
        end  % function        
    end  % public methods block

    methods (Access=protected)
        function [acquisitionKeystoneTaskType, acquisitionKeystoneTaskDeviceName, stimulationKeystoneTaskType, stimulationKeystoneTaskDeviceName] = ...
                determineKeystoneTasks_(self)
            % The acq and stim subsystems each have a "keystone" task.  This task is
            % identified by its type (one of "ai", "di", "ao", and "do") and its device
            % name (the name of one of the devices in use).  In some cases, the
            % keystone task for the acq subsystem is the same as that for the stim
            % subsystem.  All the tasks in the subsystem that are not the keystone task
            % have their start trigger set to <keystone task name>/<keystone task
            % type>/StartTrigger.  If a task is a keystone task, it is started after
            % all non-keystone tasks are started.
            %
            % If you're not careful about this stuff, the acquisition tasks can end up
            % getting triggered (e.g. by an external trigger) before the stimulation
            % tasks have been started, even though the user has configured both acq and
            % stim subsystems to use the same trigger.  This business with the keystone
            % tasks is designed to eliminate this in the common case of acq and stim
            % subsystems using the same trigger, and ameliorate it in cases where the
            % acq and stim subsystems use different triggers.
            
            % First figure out the acq keystone task
            isAIChannelActive = self.IsAIChannelActive ;
            deviceNamePerActiveChannel = self.AIChannelDeviceNames(isAIChannelActive) ;
            terminalIDPerActiveChannel = self.AIChannelTerminalIDs(isAIChannelActive) ;                
            % ws.collectTerminalsByDevice only returns device names that appear in deviceNamePerActiveChannel
            deviceNamePerActiveAIDevice = ...
                ws.collectTerminalsByDevice(deviceNamePerActiveChannel, terminalIDPerActiveChannel, self.PrimaryDeviceName) ;
            if isempty(deviceNamePerActiveAIDevice) ,
                nDIChannels = self.Acquisition_.NActiveDigitalChannels ;            
                if nDIChannels==0 ,
                    % If get here, no active AI channels, no active DI channels. WS will throw
                    % an error at run start in this case, so doesn't really matter what we
                    % return.
                    acquisitionKeystoneTaskType = '' ;
                    acquisitionKeystoneTaskDeviceName = '' ;
                else
                    acquisitionKeystoneTaskType = 'di' ;
                    acquisitionKeystoneTaskDeviceName = self.PrimaryDeviceName ;
                end
            else
                % There's at least one active AI channel
                acquisitionKeystoneTaskType = 'ai' ;
                acquisitionKeystoneTaskDeviceName = deviceNamePerActiveAIDevice{1} ;  
                  % this will be the primary device, if there are any active AI channels on it
            end
            
            % Now figure out the stim keystone task
            if self.AcquisitionTriggerIndex==self.StimulationTriggerIndex ,
                % Acq and stim subsystems are using the same trigger, so
                % acq and stim subsystems will have the same keystone task+device.
                stimulationKeystoneTaskType = acquisitionKeystoneTaskType ;
                stimulationKeystoneTaskDeviceName = acquisitionKeystoneTaskDeviceName ;
            else
                % Acq and stim subsystems are using different triggers, so
                % acq and stim subsystems will have distinct keystone task+device.
                
                % So now we have to determine the stim keystone tasks.                
                deviceNamePerChannel = self.AOChannelDeviceNames ;
                terminalIDPerChannel = self.AOChannelTerminalIDs ;
                % ws.collectTerminalsByDevice only returns device names that appear in deviceNamePerActiveChannel
                deviceNamePerAODevice = ...
                    ws.collectTerminalsByDevice(deviceNamePerChannel, terminalIDPerChannel, self.PrimaryDeviceName) ;
                if isempty(deviceNamePerAODevice) ,
                    nDOChannels = self.Stimulation_.NTimedDigitalChannels ;  
                    if nDOChannels == 0 ,
                        % No AO channel, no timed DO channels.
                        % In this case, no timed output tasks will be created, so having these
                        % empty shouldn't be a problem.
                        stimulationKeystoneTaskType = '' ;
                        stimulationKeystoneTaskDeviceName = '' ;
                    else
                        stimulationKeystoneTaskType = 'do' ;
                        stimulationKeystoneTaskDeviceName = self.PrimaryDeviceName ;
                    end
                else
                    stimulationKeystoneTaskType = 'ao' ;
                    stimulationKeystoneTaskDeviceName = deviceNamePerAODevice{1} ;
                      % this will be the primary device, if there are any active AO channels on it
                end
            end
            
%             fprintf('In WavesurferModel:determineKeystoneTasks():\n') ;
%             fprintf('  acquisitionKeystoneTask: %s\n', acquisitionKeystoneTaskType) ;
%             fprintf('  acquisitionKeystoneDevice: %s\n', acquisitionKeystoneTaskDeviceName) ;
%             fprintf('  stimulationKeystoneTask: %s\n', stimulationKeystoneTaskType) ;            
%             fprintf('  stimulationKeystoneDevice: %s\n', stimulationKeystoneTaskDeviceName) ;            
        end
        
        function mimicProtocolThatWasJustLoaded_(self, other)
            % Cause self to resemble other, but only w.r.t. the protocol.

            % Do this before replacing properties in place, or bad things
            % will happen
            self.releaseTimedHardwareResources_() ;

            % Get the list of property names for this file type
            propertyNames = ws.listPropertiesForPersistence(self);
            
            % Don't want to do broadcasts while we're in a
            % possibly-inconsistent state
            %self.disableBroadcasts() ;
            self.disableBroadcasts() ;
            
            % Set each property to the corresponding one
            for i = 1:length(propertyNames) ,
                thisPropertyName=propertyNames{i};
                if any(strcmp(thisPropertyName,{'Triggering_', 'Acquisition_', 'Stimulation_', 'Display_', 'Ephys_'})) ,
                    %self.(thisPropertyName).mimic(other.(thisPropertyName)) ;
                    self.(thisPropertyName).mimic(other.getPropertyValue_(thisPropertyName)) ;
                elseif any(strcmp(thisPropertyName,{'UserCodeManager_'})) ,
                    %self.(thisPropertyName).mimic(other.(thisPropertyName)) ;
                    self.(thisPropertyName).mimic(other.getPropertyValue_(thisPropertyName)) ;  % needs root model arg (not anymore...)
                elseif any(strcmp(thisPropertyName,{'FastProtocols_', 'Logging_'})) ,
                    % do nothing                
                elseif any(strcmp(thisPropertyName, ...
                                  {'MainFigurePosition_', 'GeneralSettingsFigurePosition_', 'ChannelsFigurePosition_', ...
                                   'StimulusLibraryFigurePosition_', 'StimulusPreviewFigurePosition_', 'TriggersFigurePosition_', ...
                                   'UserCodeManagerFigurePosition_', 'ElectrodeManagerFigurePosition_', 'TestPulserFigurePosition_'})) ,
                    if isprop(other, thisPropertyName) ,
                        source = other.getPropertyValue_(thisPropertyName) ;
                        if isempty(source) ,
                            % Do nothing, no reason to replace a possibly-nonempty position with an empty
                            % one.  Empty ones sometimes happen for old protocol files...
                        else
                            self.setPropertyValue_(thisPropertyName, source) ;
                        end
                    end
                else
                    if isprop(other,thisPropertyName) ,
                        source = other.getPropertyValue_(thisPropertyName) ;
                        self.setPropertyValue_(thisPropertyName, source) ;
                    end
                end
            end
            
            % Do sanity-checking on persisted state
            self.sanitizePersistedState_() ;
            
            % Make sure the transient state is consistent with
            % the non-transient state
            self.synchronizeTransientStateToPersistedState_() ;            
            
            % Do some things to help with old protocol files.  Have to do
            % it here b/c have to do it after
            % synchronizeTransientStateToPersistedState_(), since that's
            % what probes the device to get number of counters, PFI lines,
            % etc, and we need that info to be right before we call the
            % methods below.
            primaryDeviceName = self.PrimaryDeviceName_ ;
            isPrimaryDeviceAPXIDevice = ws.isDeviceAPXIDevice(primaryDeviceName) ;
            self.IsPrimaryDeviceAPXIDevice_ = isPrimaryDeviceAPXIDevice ;
            deviceIndex = self.getDeviceIndexFromName(primaryDeviceName) ;
            if isempty(deviceIndex) ,
                % this means the device name does not specify a currently-valid device name
                nCounters = 0 ;
                nPFITerminals = 0 ;
                %nAITerminals = 0 ;
                %nAOTerminals = 0 ;                                    
                %nDIOTerminals = 0 ;
            else
                nCounters = self.NCountersPerDevice_(deviceIndex) ;                
                nPFITerminals = self.NPFITerminalsPerDevice_(deviceIndex) ;                
                %nAITerminals = self.NAITerminalsPerDevice_(deviceIndex) ;
                %nAOTerminals = self.NAOTerminalsPerDevice_(deviceIndex) ;          
                %nDIOTerminals = self.NDIOTerminalsPerDevice_(deviceIndex) ;
            end
            self.informSubsystemsThatWeAreSettingPrimaryDeviceName_(primaryDeviceName, nCounters, nPFITerminals) ;  
                % Old protocol files don't store the 
                % device name in subobjects, so we call
                % this to set the PrimaryDeviceName
                % throughout to the one set in self
            self.didSetAreSweepsFiniteDuration_(self.AreSweepsFiniteDuration_, self.NSweepsPerRun_) ;  % Ditto
            self.didSetNSweepsPerRun_(self.NSweepsPerRun_) ;  % Ditto

            % Safe to do broadcasts again
            %self.enableBroadcastsMaybe() ;
            self.enableBroadcastsMaybe() ;
            
            % Make sure the looper knows which output channels are timed vs
            % on-demand
            %keyboard
            if self.IsAwake_ ,
                %isTerminalOvercommittedForEachDOChannel = self.IsDOChannelTerminalOvercommitted ;  % this is transient, so isn't in the wavesurferModelSettings
%                 self.IPCPublisher_.send('didSetPrimaryDeviceInFrontend', ...
%                                         primaryDeviceName, ...
%                                         isPrimaryDeviceAPXIDevice, ...
%                                         isTerminalOvercommittedForEachDOChannel) ;
                self.Looper_.didSetPrimaryDeviceInFrontend(self.PrimaryDeviceName, ...
                                                           self.IsPrimaryDeviceAPXIDevice, ...
                                                           self.DOChannelTerminalIDs, ...
                                                           self.IsDOChannelTimed, ...
                                                           self.DOChannelStateIfUntimed, ...
                                                           self.IsDOChannelTerminalOvercommitted) ;
                self.Refiller_.didSetPrimaryDeviceInFrontend() ;
                %looperProtocol = self.getLooperProtocol_() ;
                %self.IPCPublisher_.send('frontendJustLoadedProtocol', looperProtocol, isTerminalOvercommittedForEachDOChannel) ;
                self.Looper_.loadingProtocol(self.PrimaryDeviceName, ...
                                             self.IsPrimaryDeviceAPXIDevice, ...
                                             self.DOChannelTerminalIDs, ...
                                             self.IsDOChannelTimed, ...
                                             self.DOChannelStateIfUntimed, ...
                                             self.IsDOChannelTerminalOvercommitted) ;
                %self.Refiller_.frontendJustLoadedProtocol(isTerminalOvercommittedForEachDOChannel) ;                
            end
        end  % function
    end  % protected methods block
    
    methods
        function sanitizePersistedState_(self)
            % This method should perform any sanity-checking that might be
            % advisable after loading the persistent state from disk.
            % This is often useful to provide backwards compatibility
            
            % Set the override state for the stimulus map durations
            self.overrideOrReleaseStimulusMapDurationAsNeeded_() ;
            
            self.Display_.sanitizePersistedStateGivenChannelCounts_(self.NAIChannels, self.NDIChannels) ;            
        end
    end  % protected methods block
    
    methods (Access=protected) 
        function mimicUserSettings_(self, other)
            % Cause self to resemble other, but only w.r.t. the user settings            
            source = other.getPropertyValue_('FastProtocols_') ;
            self.FastProtocols_ = ws.copyCellArrayOfHandles(source) ;
        end  % function        
    end  % protected methods block
    
    methods (Access=protected)    
%         function disableAllBroadcastsDammit_(self)
%             self.disableBroadcasts() ;
%             %self.Triggering_.disableBroadcasts() ;
%             %self.Acquisition_.disableBroadcasts() ;
%             %self.Stimulation_.disableBroadcasts() ;
%             %self.Display_.disableBroadcasts() ;
%             %self.Ephys_.TestPulser.disableBroadcasts() ;
%             %self.Ephys_.ElectrodeManager.disableBroadcasts() ;
%             %self.Ephys_.disableAllBroadcastsDammit_() ;
%             %self.UserCodeManager_.disableBroadcasts() ;
%             %self.Logging_.disableBroadcasts() ;            
%         end
%         
%         function enableBroadcastsMaybeDammit_(self)
%             self.Logging_.enableBroadcastsMaybe() ;                        
%             self.UserCodeManager_.enableBroadcastsMaybe() ;
%             self.Ephys_.enableBroadcastsMaybeDammit_() ;
% %             self.Ephys_.ElectrodeManager.enableBroadcastsMaybe() ;
% %             self.Ephys_.TestPulser.enableBroadcastsMaybe() ;
%             self.Display_.enableBroadcastsMaybe() ;
%             self.Stimulation_.enableBroadcastsMaybe() ;
%             self.Acquisition_.enableBroadcastsMaybe() ;            
%             self.Triggering_.enableBroadcastsMaybe() ;
%             self.enableBroadcastsMaybe() ;
%         end
        
%         function updateEverythingAfterProtocolFileOpen_(self)
%             self.broadcast('UpdateLogging') ;
%             self.broadcast('UpdateUserCodeManager') ;
%             self.broadcast('UpdateElectrodeManager') ;
%             self.broadcast('UpdateTestPulser') ;            
%             self.broadcast('UpdateTraces') ;
%             self.broadcast('UpdateDisplay') ;
%             self.broadcast('UpdateStimulusLibrary') ;
%             self.broadcast('UpdateTriggering') ;
%             self.broadcast('Update') ;            
%             self.broadcast('LayoutAllWindows') ;
%         end
    end  % protected methods block
    
    methods (Static)
        function pathToRepoRoot = pathNamesThatNeedToBeOnSearchPath()
            % Allow user to invoke Wavesurfer from the Matlab command line, for
            % this Matlab session only.

            pathToWavesurferModel = mfilename('fullpath') ;
            pathToWsModulerFolder = fileparts(pathToWavesurferModel) ;  % should be +ws folder
            pathToRepoRoot = fileparts(pathToWsModulerFolder) ;  % should be repo root
            %pathToMatlabZmqLib = fullfile(pathToRepoRoot,'matlab-zmq','lib') ;
            
            %result = { pathToRepoRoot , pathToMatlabZmqLib } ;
        end
        
%         function portNumbers = getFreeEphemeralPortNumbers(nPorts)
%             % Determine which three free ports to use:
%             % First bind to three free ports to get their addresses
%             freePorts = struct('context',{},'socket',{},'endpoint',{},'portNumber',{});
%             for i=1:nPorts ,
%                 %freePorts(i).context = zmq.core.ctx_new();
%                 freePorts(i).context = zmq.Context();
%                 %freePorts(i).socket  = zmq.core.socket(freePorts(i).context, 'ZMQ_PUSH');
%                 freePorts(i).socket  = freePorts(i).context.socket('ZMQ_PUSH');
%                 address = 'tcp://127.0.0.1:*';
%                 %zmq.core.bind(freePorts(i).socket, address);
%                 freePorts(i).socket.bind(address);
%                 %freePorts(i).endpoint = zmq.core.getsockopt(freePorts(i).socket, 'ZMQ_LAST_ENDPOINT');
%                 %freePorts(i).endpoint = freePorts(i).socket.getsockopt('ZMQ_LAST_ENDPOINT');
%                 %freePorts(i).endpoint = freePorts(i).socket.get('ZMQ_LAST_ENDPOINT');
%                 freePorts(i).endpoint = freePorts(i).socket.bindings{end} ;
%                 splitString = strsplit(freePorts(i).endpoint,'tcp://127.0.0.1:');
%                 freePorts(i).portNumber = str2double(splitString{2});
%             end
% 
%             % Unbind the ports to free them up for the actual
%             % processes. Doing it in this way (rather than
%             % binding/unbinding each port sequentially) will minimize amount of
%             % time between a port being unbound and bound by a process.
%             for i=1:nPorts ,
%                 %zmq.core.disconnect(freePorts(i).socket, freePorts(i).endpoint);
%                 %zmq.core.close(freePorts(i).socket);
%                 %zmq.core.ctx_shutdown(freePorts(i).context);
%                 %zmq.core.ctx_term(freePorts(i).context);
%                 freePorts(i).socket = [] ;
%                 freePorts(i).context = [] ;
%             end
% 
%             portNumbers = [ freePorts(:).portNumber ] ;            
%         end
    end  % static methods block
    
    methods
%         function varargout = execute(methodName, varargin)
%             % A method call, but with logging of warnings
%             self.clearWarningLog_() ;
%             self.startWarningLogging_() ;
%             varargout = self.(methodName)(varargin{:}) ;
%             self.stopWarningLogging_() ;            
%         end  % method
%         
%         function result = didWarningsOccur(self) 
%             result = ~isempty(self.WarningLog_) ;
%         end
%         
%         function result = getWarningLog(self)
%             result = self.WarningLog_ ;
%         end  % method
        
        function do(self, methodName, varargin)
            % This is intended to be the usual way of calling model
            % methods.  For instance, a call to a ws.Controller
            % controlActuated() method should generally result in a single
            % call to .do() on it's model object, and zero direct calls to
            % model methods.  This gives us a
            % good way to implement functionality that is common to all
            % model method calls, when they are called as the main "thing"
            % the user wanted to accomplish.  For instance, we start
            % warning logging near the beginning of the .do() method, and turn
            % it off near the end.  That way we don't have to do it for
            % each model method, and we only do it once per user command.            
            self.AllowTimerCallback_ = false ;
            self.startLoggingWarnings() ;
            try
                self.(methodName)(varargin{:}) ;
            catch exception
                % If there's a real exception, the warnings no longer
                % matter.  But we want to restore the model to the
                % non-logging state.
                self.stopLoggingWarnings() ;  % discard the result, which might contain warnings
                self.resetReadiness_() ;  % Need to do this to make sure we don't stay unready for the rest of the WSM lifetime
                self.AllowTimerCallback_ = true ;
                rethrow(exception) ;
            end
            warningExceptionMaybe = self.stopLoggingWarnings() ;
            if ~isempty(warningExceptionMaybe) ,
                warningException = warningExceptionMaybe{1} ;
                self.AllowTimerCallback_ = true ;
                throw(warningException) ;
            end
            self.AllowTimerCallback_ = true ;
        end

        function logWarning(self, identifier, message, causeOrEmpty)
            % This is public b/c subsystem need to call it, but it should
            % only be called by subsystems.
            if nargin<4 ,
                causeOrEmpty = [] ;
            end
            warningException = MException(identifier, message) ;
            if ~isempty(causeOrEmpty) ,
                warningException = warningException.addCause(causeOrEmpty) ;
            end
            if self.UnmatchedLogWarningStartCount_>0 ,
                self.WarningCount_ = self.WarningCount_ + 1 ;
                if self.WarningCount_ < 10 ,
                    self.WarningLog_ = vertcat(self.WarningLog_, ...
                                               warningException) ;            
                else
                    % Don't want to log a bazillion warnings, so do nothing
                end
            else
                % Just issue a normal warning
                warning(identifier, message) ;
                if ~isempty(causeOrEmpty) ,
                    fprintf('Cause of warning:\n');
                    display(causeOrEmpty.getReport());
                end
            end
        end  % method
    end  % public methods block
    
    methods 
        function startLoggingWarnings(self)
            % fprintf('\n\n\n\n');
            % dbstack
            % fprintf('At entry to startLogginWarnings: self.UnmatchedLogWarningStartCount_ = %d\n', self.UnmatchedLogWarningStartCount_) ;
            self.UnmatchedLogWarningStartCount_ = self.UnmatchedLogWarningStartCount_ + 1 ;
            if self.UnmatchedLogWarningStartCount_ <= 1 ,
                % Unless things have gotten weird,
                % self.UnmatchedLogWarningStartCount_ should *equal* one
                % here.
                % This means this is the first unmatched start (or the
                % first one after a matched set of starts and stops), so we
                % reset the warning count and the warning log.
                self.UnmatchedLogWarningStartCount_ = 1 ;  % If things have gotten weird, fix them.
                self.WarningCount_ = 0 ;
                self.WarningLog_ = MException.empty(0, 1) ;
            end
        end        
        
        function exceptionMaybe = stopLoggingWarnings(self)
            % fprintf('\n\n\n\n');
            % dbstack
            % fprintf('At entry to stopLogginWarnings: self.UnmatchedLogWarningStartCount_ = %d\n', self.UnmatchedLogWarningStartCount_) ;
            % Get the warnings, if any
            self.UnmatchedLogWarningStartCount_ = self.UnmatchedLogWarningStartCount_ - 1 ;
            if self.UnmatchedLogWarningStartCount_ <= 0 , 
                % Technically, this should only happen when
                % self.UnmatchedLogWarningStartCount_==0, unless there was a
                % programming error.  But in any case, we
                % produce a summary of the warnings, if any, and set the
                % (unmatched) start counter to zero, even though it should
                % be zero already.
                self.UnmatchedLogWarningStartCount_ = 0 ;
                loggedWarnings = self.WarningLog_ ;
                % Process them, summarizing them in a maybe (a list of length
                % zero or one) of exceptions.  The individual warnings are
                % stored in the causes of the exception, if it exists.
                nWarnings = self.WarningCount_ ;
                if nWarnings==0 ,
                    exceptionMaybe = {} ;                
                else
                    if nWarnings==1 ,
                        exceptionMessage = loggedWarnings(1).message ;
                    else
                        exceptionMessage = sprintf('%d warnings occurred.\nThe first one was: %s', ...
                                                   nWarnings, ...
                                                   loggedWarnings(1).message) ;
                    end                                           
                    exception = MException('ws:warningsOccurred', exceptionMessage) ;
                    for i = 1:length(loggedWarnings) ,
                        exception = exception.addCause(loggedWarnings(i)) ;
                    end
                    exceptionMaybe = {exception} ;
                end
                % Clear the warning log before returning
                self.WarningCount_ = 0 ;
                self.WarningLog_ = MException.empty(0, 1) ;
            else
                % self.UnmatchedLogWarningStartCount_ > 1, so decrement the number of
                % unmatched starts, but don't do much else.
                exceptionMaybe = {} ;
            end
        end  % method
        
        function clearSelectedFastProtocol(self)
            selectedIndex = self.IndexOfSelectedFastProtocol ;  % this is never empty
            fastProtocol = self.FastProtocols_{selectedIndex} ;
            fastProtocol.ProtocolFileName = '' ;
            fastProtocol.AutoStartType = 'do_nothing' ;
            self.broadcast('UpdateFastProtocols');
            self.broadcast('UpdateMain');
        end  % method
        
%         function result = selectedFastProtocolFileName(self)
%             % Returns a maybe of strings (i.e. a cell array of string, the
%             % cell array of length zero or one)
%             selectedIndex = self.IndexOfSelectedFastProtocol ;  % this is never empty
%             fastProtocol = self.FastProtocols{selectedIndex} ;
%             result = fastProtocol.ProtocolFileName ;            
%         end  % method
        
%         function setSelectedFastProtocolFileName(self, newFileName)
%             % newFileName should be an absolute file path
%             selectedIndex = self.IndexOfSelectedFastProtocol ;  % this is never empty
%             self.setFastProtocolFileName(selectedIndex, newFileName) ;
%         end  % method
        
%         function setFastProtocolFileName(self, index, newFileName)
%             % newFileName should be an absolute file path
%             if isscalar(index) && isnumeric(index) && isreal(index) && round(index)==index && 1<=index && index<=self.NFastProtocols ,                
%                 if ws.isString(newFileName) && ~isempty(newFileName) && ws.isFileNameAbsolute(newFileName) && exist(newFileName,'file') ,
%                     fastProtocol = self.FastProtocols_{index} ;
%                     fastProtocol.ProtocolFileName = newFileName ;
%                 else
%                     self.updateFastProtocol() ;
%                     error('ws:invalidPropertyValue', ...
%                           'Fast protocol file name must be an absolute path');              
%                 end
%             else
%                 self.updateFastProtocol() ;
%                 error('ws:invalidPropertyValue', ...
%                       'Fast protocol index must a real numeric scalar integer between 1 and %d', self.NFastProtocols);
%             end                
%             self.updateFastProtocol() ;
%         end  % method

        function result = getSelectedFastProtocolProperty(self, propertyName)
            selectedIndex = self.IndexOfSelectedFastProtocol ;  % this is never empty
            result = self.getFastProtocolProperty(selectedIndex, propertyName) ;
        end
        
        function setSelectedFastProtocolProperty(self, propertyName, newValue)
            selectedIndex = self.IndexOfSelectedFastProtocol ;  % this is never empty
            self.setFastProtocolProperty(selectedIndex, propertyName, newValue) ;
        end
        
        function result = getFastProtocolProperty(self, index, propertyName)
            if isscalar(index) && isnumeric(index) && isreal(index) && round(index)==index && 1<=index && index<=self.NFastProtocols ,               
                fastProtocol = self.FastProtocols_{index} ;
                result = fastProtocol.(propertyName) ;
            else
                error('ws:invalidPropertyValue', ...
                      'Fast protocol index must a real numeric scalar integer between 1 and %d', self.NFastProtocols);
            end                
        end        
        
        function setFastProtocolProperty(self, index, propertyName, newValue)
            if isscalar(index) && isnumeric(index) && isreal(index) && round(index)==index && 1<=index && index<=self.NFastProtocols ,               
                fastProtocol = self.FastProtocols_{index} ;
                try 
                    fastProtocol.(propertyName) = newValue ;
                    %if isequal(propertyName, 'ProtocolFileName') && self.DoUsePreferences ,
                    %    ws.setPreference('LastProtocolFilePath', newValue) ;
                    %end
                catch exception
                    self.broadcast('UpdateFastProtocols');
                    self.broadcast('UpdateMain');
                    rethrow(exception) ;
                end
            else
                self.broadcast('UpdateFastProtocols');
                self.broadcast('UpdateMain');  % need to update main window also                    
                error('ws:invalidPropertyValue', ...
                      'Fast protocol index must a real numeric scalar integer between 1 and %d', self.NFastProtocols);
            end                
            self.broadcast('UpdateFastProtocols');
            self.broadcast('UpdateMain');  % need to update main window also            
        end  % method        
                
        function incrementSessionIndex(self)
            self.Logging_.incrementSessionIndex() ;
        end
        
        function setSelectedOutputableByIndex(self, index)            
            self.Stimulation_.setSelectedOutputableByIndex(index) ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateGeneral');
        end  % method

        function setSelectedOutputableByClassNameAndIndex(self, className, indexWithinClass)
            self.Stimulation_.setSelectedOutputableByClassNameAndIndex(className, indexWithinClass) ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateGeneral');
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % method        
        
        function openFastProtocolByIndex(self, index)
            if ws.isIndex(index) && 1<=index && index<=self.NFastProtocols ,
                fastProtocol = self.FastProtocols_{index} ;
                fileName = fastProtocol.ProtocolFileName ;
                if ~isempty(fileName) , 
                    if exist(fileName, 'file') ,
                        self.openProtocolFileGivenFileName(fileName) ;
                    else
                        error('ws:fastProtocolFileMissing', ...
                              'The protocol file %s is missing.',fileName);
                    end
                end
            end
        end
        
        function performAutoStartForFastProtocolByIndex(self, index) 
            if ws.isIndex(index) && 1<=index && index<=self.NFastProtocols ,
                fastProtocol = self.FastProtocols_{index} ;            
                if isequal(fastProtocol.AutoStartType,'play') ,
                    self.play();
                elseif isequal(fastProtocol.AutoStartType,'record') ,
                    self.record();
                end
            end
        end  % method        
        
%         function result = get.LayoutForAllWindows(self)
%             result = self.LayoutForAllWindows_ ;
%         end
        
        function value = get.AllDeviceNames(self)
            value = self.AllDeviceNames_ ;
        end  % function

        function value = get.PrimaryDeviceName(self)
            value = self.PrimaryDeviceName_ ;
        end  % function

        function value = get.IsPrimaryDeviceAPXIDevice(self)
            value = self.IsPrimaryDeviceAPXIDevice_ ;
        end  % function
        
        function set.PrimaryDeviceName(self, newValue)
            if ws.isString(newValue) && ~isempty(newValue) ,
                allDeviceNames = self.AllDeviceNames ;
                isAMatch = strcmpi(newValue,allDeviceNames) ;  % DAQmx device names are not case-sensitive
                if any(isAMatch) ,
                    iMatch = find(isAMatch,1) ;
                    primaryDeviceName = allDeviceNames{iMatch} ;
                    self.PrimaryDeviceName_ = primaryDeviceName ;
                    isPrimaryDeviceAPXIDevice = ws.isDeviceAPXIDevice(primaryDeviceName) ;
                    self.IsPrimaryDeviceAPXIDevice_ = isPrimaryDeviceAPXIDevice ;                    
                    self.DoesProtocolNeedSave_ = true ;
                    
                    % Tell the subsystems that we've changed the device
                    % name
                    nCounters = self.NCountersPerDevice_(iMatch) ;
                    nPFITerminals = self.NPFITerminalsPerDevice_(iMatch) ;
                    self.informSubsystemsThatWeAreSettingPrimaryDeviceName_(primaryDeviceName, nCounters, nPFITerminals) ;

                    % Recalculate which digital terminals are now
                    % overcommitted, since that also updates which are
                    % out-of-range for the device
                    self.syncIsAIChannelTerminalOvercommitted_() ;
                    self.syncIsAOChannelTerminalOvercommitted_() ;
                    self.syncIsDIOChannelTerminalOvercommitted_() ;
                    
                    % Change our state to reflect the presence of the
                    % device
                    self.setState_('idle') ;

                    if self.IsAwake_ ,
                        self.Looper_.didSetPrimaryDeviceInFrontend(self.PrimaryDeviceName, ...
                                                                   self.IsPrimaryDeviceAPXIDevice, ...
                                                                   self.DOChannelTerminalIDs, ...
                                                                   self.IsDOChannelTimed, ...
                                                                   self.DOChannelStateIfUntimed, ...
                                                                   self.IsDOChannelTerminalOvercommitted) ;
                        self.Refiller_.didSetPrimaryDeviceInFrontend(primaryDeviceName, ...
                                                                     isPrimaryDeviceAPXIDevice, ...
                                                                     self.IsDOChannelTerminalOvercommitted) ;
                    end                        
                else
                    self.broadcast('UpdateMain');
                    self.broadcast('UpdateGeneral') ;
                    self.broadcast('UpdateChannels') ;
                    self.broadcast('UpdateTriggering') ;
                    error('ws:invalidPropertyValue', ...
                          'PrimaryDeviceName must be the name of an NI DAQmx device');       
                end                        
            else
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;
                self.broadcast('UpdateTriggering') ;
                error('ws:invalidPropertyValue', ...
                      'PrimaryDeviceName must be a nonempty string');       
            end
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
            self.broadcast('UpdateTriggering') ;
        end  % function
    end  % public methods block
    
    methods
        function result = get.NAITerminalsPerDevice(self)
            % The number of AI channels available, if you used them all in
            % differential mode, which is what we do.
            result = self.NAITerminalsPerDevice_ ;
        end
        
        function result = get.AITerminalIDsOnEachDevice(self)
            % A list of the available AI terminal IDs on the current
            % device, if you used all AIs in differential mode, which is
            % what we do.
            result = self.AITerminalIDsOnEachDevice_ ;
        end
        
        function result = get.NAOTerminalsPerDevice(self)
            % The number of AO channels available.
            result = self.NAOTerminalsPerDevice_ ;
        end

        function result = get.NDIOTerminalsPerDevice(self)
            % The number of DIO channels available.  We only count the DIO
            % channels capable of timed operation, i.e. the P0.x channels.
            % This is a conscious design choice.  We treat the PFIn/Pm.x
            % channels as being only PFIn channels.
            result = self.NDIOTerminalsPerDevice_ ;
        end  % function
        
        function result = get.NPFITerminalsPerDevice(self)
            result = self.NPFITerminalsPerDevice_ ;
        end  % function
        
        function result = get.NCountersPerDevice(self)
            % The number of counters (CTRs) on the board.
            result = self.NCountersPerDevice_ ;
        end  % function        
        
        function result = getAllAITerminalNames(self, deviceName)             
            deviceIndex = self.getDeviceIndexFromName(deviceName) ;
            if isempty(deviceIndex) ,
                allAITerminalIDs = zeros(1,0) ;
            else
                nAIsInHardware = self.NAITerminalsPerDevice(deviceIndex) ;  % this is the number of terminals if all are differential, which they are
                allAITerminalIDs = ws.differentialAITerminalIDsGivenCount(nAIsInHardware) ;
            end
            result = arrayfun(@(id)(sprintf('AI%d',id)), allAITerminalIDs, 'UniformOutput', false ) ;
        end        
        
        function result = getAllAOTerminalNames(self, deviceName)             
            deviceIndex = self.getDeviceIndexFromName(deviceName) ;
            if isempty(deviceIndex) ,
                nAOsInHardware = 0 ;
            else
                nAOsInHardware = self.NAOTerminalsPerDevice(deviceIndex) ;
            end
            result = arrayfun(@(id)(sprintf('AO%d',id)), 0:(nAOsInHardware-1), 'UniformOutput', false ) ;
        end
        
        function result = getAllDIOTerminalNames(self, deviceName)             
            deviceIndex = self.getDeviceIndexFromName(deviceName) ;
            if isempty(deviceIndex) ,
                nChannelsInHardware = 0 ;
            else
                nChannelsInHardware = self.NDIOTerminalsPerDevice(deviceIndex) ;
            end
            result = arrayfun(@(id)(sprintf('P0.%d',id)), 0:(nChannelsInHardware-1), 'UniformOutput', false ) ;
        end        
        
%         function result = get.NDigitalChannels(self)
%             nDIs = self.Acquisition_.NDigitalChannels ;
%             nDOs = self.Stimulation_.NDigitalChannels ;
%             result =  nDIs + nDOs ;
%         end
        
        function result = get.AllChannelNames(self)
            aiNames = self.Acquisition_.AnalogChannelNames ;
            diNames = self.Acquisition_.DigitalChannelNames ;
            aoNames = self.Stimulation_.AnalogChannelNames ;
            doNames = self.Stimulation_.DigitalChannelNames ;
            result = [ aiNames diNames aoNames doNames ] ;
        end

        function result = get.AIChannelNames(self)
            result = self.Acquisition_.AnalogChannelNames ;
        end

        function result = get.DIChannelNames(self)
            result = self.Acquisition_.DigitalChannelNames ;
        end

        function result = get.AOChannelNames(self)
            result = self.Stimulation_.AnalogChannelNames ;
        end

        function result = get.DOChannelNames(self)
            result = self.Stimulation_.DigitalChannelNames ;
        end

        function result = get.IsAIChannelTerminalOvercommitted(self)
            result = self.IsAIChannelTerminalOvercommitted_ ;
        end
        
        function result = get.IsAOChannelTerminalOvercommitted(self)
            result = self.IsAOChannelTerminalOvercommitted_ ;
        end
        
        function result = get.IsDIChannelTerminalOvercommitted(self)
            result = self.IsDIChannelTerminalOvercommitted_ ;
        end
        
        function result = get.IsDOChannelTerminalOvercommitted(self)
            result = self.IsDOChannelTerminalOvercommitted_ ;
        end        
    end  % public methods block
        
    methods (Access=protected)
        function syncDeviceResourceCountsFromDeviceNames_(self)
            % Probe the devices to find out their capabilities
            allDeviceNames = self.AllDeviceNames ;
            nDevices = length(allDeviceNames) ;
            nDIOTerminalsPerDevice = zeros(1, nDevices) ;
            nPFITerminalsPerDevice = zeros(1, nDevices) ;
            nCountersPerDevice = zeros(1, nDevices) ;
            nAITerminalsPerDevice = zeros(1, nDevices) ;
            nAOTerminalsPerDevice = zeros(1, nDevices) ;
            aiTerminalIDsOnEachDevice = cell(1, nDevices) ;
            for i = 1:nDevices ,            
                deviceName = allDeviceNames{i} ;
                [nDIOTerminals, nPFITerminals] = ws.getNumberOfDIOAndPFITerminalsFromDevice(deviceName) ;
                nCounters = ws.getNumberOfCountersFromDevice(deviceName) ;
                nAITerminals = ws.getNumberOfDifferentialAITerminalsFromDevice(deviceName) ;
                aiTerminalsOnDevice = ws.differentialAITerminalIDsGivenCount(nAITerminals) ;
                nAOTerminals = ws.getNumberOfAOTerminalsFromDevice(deviceName) ;            
                % Put in arrays
                nDIOTerminalsPerDevice(i) = nDIOTerminals ;                
                nPFITerminalsPerDevice(i) = nPFITerminals ;
                nCountersPerDevice(i) = nCounters ;
                nAITerminalsPerDevice(i) = nAITerminals ;
                nAOTerminalsPerDevice(i) = nAOTerminals ;                
                aiTerminalIDsOnEachDevice{i} = aiTerminalsOnDevice ;
            end
            self.NDIOTerminalsPerDevice_ = nDIOTerminalsPerDevice ;
            self.NPFITerminalsPerDevice_ = nPFITerminalsPerDevice ;
            self.NCountersPerDevice_ = nCountersPerDevice ;
            self.NAITerminalsPerDevice_ = nAITerminalsPerDevice ;
            self.NAOTerminalsPerDevice_ = nAOTerminalsPerDevice ;
            self.AITerminalIDsOnEachDevice_ = aiTerminalIDsOnEachDevice ;
        end  % function     
        
        function syncIsAIChannelTerminalOvercommitted_(self)            
            % For each channel, determines if the terminal ID for that
            % channel is "overcommited".  I.e. if two channels specify the
            % same terminal ID, that terminal ID is overcommitted.  Also,
            % if that specified terminal ID is not a legal terminal ID for
            % the current device, then we say that that terminal ID is
            % overcommitted.  (Because there's one in use, and zero
            % available.)
            
            % For AI terminals
            deviceNameForEachAIChannel = self.AIChannelDeviceNames ;
            aiTerminalIDForEachAIChannel = self.Acquisition_.AnalogTerminalIDs ;
            allDeviceNames = self.AllDeviceNames ;
            aiTerminalIDsOnEachDevice = self.AITerminalIDsOnEachDevice ;
                       
            nOccurencesOfDeviceAndTerminal = ws.nOccurencesOfDeviceNameAndID(deviceNameForEachAIChannel, aiTerminalIDForEachAIChannel) ;            
            isDeviceNameAITerminalIDPairValid = ...
                ws.isDeviceNameTerminalIDPairValid(deviceNameForEachAIChannel, aiTerminalIDForEachAIChannel, allDeviceNames, aiTerminalIDsOnEachDevice) ;
            
            self.IsAIChannelTerminalOvercommitted_ = ...
                (nOccurencesOfDeviceAndTerminal>1) | ~isDeviceNameAITerminalIDPairValid ;
        end
        
        function syncIsAOChannelTerminalOvercommitted_(self)            
            % For each channel, determines if the terminal ID for that
            % channel is "overcommited".  I.e. if two channels specify the
            % same terminal ID, that terminal ID is overcommitted.  Also,
            % if that specified terminal ID is not a legal terminal ID for
            % the current device, then we say that that terminal ID is
            % overcommitted.
            
            % For AO terminals
            deviceNameForEachAOChannel = self.AOChannelDeviceNames ;
            aoTerminalIDForEachAOChannel = self.Stimulation_.AnalogTerminalIDs ;
            allDeviceNames = self.AllDeviceNames ;
            nAOTerminalsPerDevice = self.NAOTerminalsPerDevice ;
            
            nOccurencesOfDeviceAndTerminal = ws.nOccurencesOfDeviceNameAndID(deviceNameForEachAOChannel, aoTerminalIDForEachAOChannel) ;            
            isDeviceNameAOTerminalIDPairValid = ...
                ws.isDeviceNameTerminalIDPairValid(deviceNameForEachAOChannel, aoTerminalIDForEachAOChannel, allDeviceNames, nAOTerminalsPerDevice) ;            
            
            self.IsAOChannelTerminalOvercommitted_ = (nOccurencesOfDeviceAndTerminal>1) | ~isDeviceNameAOTerminalIDPairValid ;            
        end
        
        function syncIsDIOChannelTerminalOvercommitted_(self)
            [nOccurencesOfDeviceAndTerminalForEachDIChannel,nOccurencesOfDeviceAndTerminalForEachDOChannel] = self.computeDIOTerminalCommitments() ;
            deviceNameForEachDIChannel = self.DIChannelDeviceNames ;
            terminalIDForEachDIChannel = self.Acquisition_.DigitalTerminalIDs ;
            deviceNameForEachDOChannel = self.DOChannelDeviceNames ;
            terminalIDForEachDOChannel = self.Stimulation_.DigitalTerminalIDs ;
            deviceNameForEachDIOChannel = horzcat(deviceNameForEachDIChannel, deviceNameForEachDOChannel) ;
            terminalIDForEachDIOChannel = horzcat(terminalIDForEachDIChannel, terminalIDForEachDOChannel) ;
            allDeviceNames = self.AllDeviceNames ;
            nDIOTerminalsPerDevice = self.NDIOTerminalsPerDevice ;
            
            isDeviceNameTerminalIDPairValidForEachDIOChannel = ...
                ws.isDeviceNameTerminalIDPairValid(deviceNameForEachDIOChannel, terminalIDForEachDIOChannel, allDeviceNames, nDIOTerminalsPerDevice) ;            
            
            nDIChannels = length(terminalIDForEachDIChannel) ;
            isDeviceNameTerminalIDPairValidForEachDIChannel = isDeviceNameTerminalIDPairValidForEachDIOChannel(1:nDIChannels) ;
            isDeviceNameTerminalIDPairValidForEachDOChannel = isDeviceNameTerminalIDPairValidForEachDIOChannel(nDIChannels+1:end) ;
            
            self.IsDIChannelTerminalOvercommitted_ = (nOccurencesOfDeviceAndTerminalForEachDIChannel>1) | ~isDeviceNameTerminalIDPairValidForEachDIChannel ;
            self.IsDOChannelTerminalOvercommitted_ = (nOccurencesOfDeviceAndTerminalForEachDOChannel>1) | ~isDeviceNameTerminalIDPairValidForEachDOChannel ;
        end  % function
    end  % protected methods block
    
    methods
        function [nOccurencesOfDeviceAndTerminalForEachDIChannel, nOccurencesOfDeviceAndTerminalForEachDOChannel] = computeDIOTerminalCommitments(self) 
            % Determine how many channels are "claiming" the terminal ID of
            % each digital channel.  On return,
            % nOccurencesOfAcquisitionTerminal is 1 x (the number of DI
            % channels) and is the number of channels that currently have
            % their terminal ID set to the same terminal ID as that channel.
            % nOccurencesOfStimulationTerminal is similar, but for DO
            % channels.
            deviceNameForEachDIChannel = self.DIChannelDeviceNames ;
            terminalIDForEachDIChannel = self.Acquisition_.DigitalTerminalIDs ;
            deviceNameForEachDOChannel = self.DOChannelDeviceNames ;
            terminalIDForEachDOChannel = self.Stimulation_.DigitalTerminalIDs ;
            deviceNameForEachDigitalChannel = horzcat(deviceNameForEachDIChannel, deviceNameForEachDOChannel) ;
            terminalIDForEachDigitalChannel = horzcat(terminalIDForEachDIChannel, terminalIDForEachDOChannel) ;
            nOccurencesOfDeviceAndTerminal = ws.nOccurencesOfDeviceNameAndID(deviceNameForEachDigitalChannel, terminalIDForEachDigitalChannel) ;
            % Sort them into the acq, stim ones
            nDIChannels = length(terminalIDForEachDIChannel) ;
            nOccurencesOfDeviceAndTerminalForEachDIChannel = nOccurencesOfDeviceAndTerminal(1:nDIChannels) ;
            nOccurencesOfDeviceAndTerminalForEachDOChannel = nOccurencesOfDeviceAndTerminal(nDIChannels+1:end) ;
        end
        
        function result = get.IsProcessingIncomingCommand(self)
            if ~isempty(self.CommandServer_) && isvalid(self.CommandServer_) ,
                result = self.CommandServer_.IsProcessingIncomingCommand ;
            else
                result = false ;
            end
        end                
    end  % public methods block

    methods (Access=protected)
        function synchronizeTransientStateToPersistedState_(self)            
            % This method should set any transient state variables to
            % ensure that the object invariants are met, given the values
            % of the persisted state variables.  The default implementation
            % does nothing, but subclasses can override it to make sure the
            % object invariants are satisfied after an object is decoded
            % from persistant storage.  This is called by
            % ws.Encodable.decodeEncodingContainerGivenParent() after
            % a new object is instantiated, and after its persistent state
            % variables have been set to the encoded values.
            
            self.Acquisition_.synchronizeTransientStateToPersistedStateHelper() ;
            self.Stimulation_.synchronizeTransientStateToPersistedStateHelper() ;
            
            %self.syncDeviceResourceCountsFromDeviceName_() ;
            %self.syncAvailableReferenceClockSourcesFromDeviceName_() ;
            self.syncIsAIChannelTerminalOvercommitted_() ;
            self.syncIsAOChannelTerminalOvercommitted_() ;
            self.syncIsDIOChannelTerminalOvercommitted_() ;
            % Need something here for yoking...
            %self.CommandClient_.IsEnabled = self.IsYokedToScanImage_ && self.IsAwake_ ;            
            %self.CommandServer_.IsEnabled = self.IsYokedToScanImage_ && self.IsAwake_ ;            
            
            self.Display_.synchronizeTransientStateToPersistedStateHelper() ;
        end  % method
        
        function informSubsystemsThatWeAreSettingPrimaryDeviceName_(self, primaryDeviceName, nCounters, nPFITerminals)
            %self.IsPrimaryDeviceAPXIDevice_ = ws.isDeviceAPXIDevice(primaryDeviceName) ;
            self.Acquisition_.settingPrimaryDeviceName(primaryDeviceName) ;
            self.Stimulation_.settingPrimaryDeviceName(primaryDeviceName) ;
            self.Triggering_.settingPrimaryDeviceName(primaryDeviceName, nCounters, nPFITerminals) ;
            %self.Ephys_.settingPrimaryDeviceName(primaryDeviceName) ;
        end  % method
        
        function didSetNSweepsPerRun_(self, nSweepsPerRun)
            self.Triggering_.didSetNSweepsPerRun(nSweepsPerRun) ;
            self.broadcast('UpdateTriggering') ;
        end  % function        
    end  % protected methods block 
    
    methods  % Present the user-relevant triggering methods at the WSM level
        function result = get.AcquisitionTriggerIndex(self)            
            result = self.Triggering_.AcquisitionTriggerSchemeIndex ;
        end  % function
        
        function set.AcquisitionTriggerIndex(self, newValue)
            try
                self.Triggering_.setAcquisitionTriggerIndex(newValue, self.NSweepsPerRun) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering');         
            self.broadcast('DidMaybeChangeProtocol') ;
        end
                
        function result = get.StimulationTriggerIndex(self)            
            result = self.Triggering_.StimulationTriggerSchemeIndex ;
        end  % function
        
        function set.StimulationTriggerIndex(self, newValue)
            try
                self.Triggering_.setStimulationTriggerIndex(newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering');                        
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function value = get.StimulationUsesAcquisitionTrigger(self)
            value = self.Triggering_.StimulationUsesAcquisitionTriggerScheme ;
        end  % function        
        
        function set.StimulationUsesAcquisitionTrigger(self, newValue)
            try
                self.Triggering_.StimulationUsesAcquisitionTriggerScheme = newValue ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.overrideOrReleaseStimulusMapDurationAsNeeded_() ;
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function        

        function result = isStimulationTriggerIdenticalToAcquisitionTrigger(self)
            result = self.Triggering_.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
        end  % function
        
        function addCounterTrigger(self)    
            try
                self.Triggering_.addCounterTrigger(self.PrimaryDeviceName) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function deleteMarkedCounterTriggers(self)
            try
                self.Triggering_.deleteMarkedCounterTriggers() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function addExternalTrigger(self)    
            try
                self.Triggering_.addExternalTrigger(self.PrimaryDeviceName) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function deleteMarkedExternalTriggers(self)
            try
                self.Triggering_.deleteMarkedExternalTriggers() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function setTriggerProperty(self, triggerType, triggerIndexWithinType, propertyName, newValue)
            try
                self.Triggering_.setTriggerProperty(triggerType, triggerIndexWithinType, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || ~isequal(propertyName, 'IsMarkedForDeletion');
            catch exception
                self.broadcast('UpdateTriggering');
                rethrow(exception) ;
            end
            self.broadcast('UpdateTriggering') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function result = get.TriggerCount(self)
            result = self.Triggering_.TriggerCount ;
        end  % function
        
        function result = get.CounterTriggerCount(self)
            result = self.Triggering_.CounterTriggerCount ;
        end  % function
        
        function result = get.ExternalTriggerCount(self)
            result = self.Triggering_.ExternalTriggerCount ;
        end  % function
        
        function result = freePFIIDs(self)
            result = self.Triggering_.freePFIIDs() ;
        end  % function
        
        function result = freeCounterIDs(self)
            result = self.Triggering_.freeCounterIDs() ;
        end  % function
        
        function result = isCounterTriggerMarkedForDeletion(self)
            result = self.Triggering_.isCounterTriggerMarkedForDeletion() ;
        end  % function
        
        function result = isExternalTriggerMarkedForDeletion(self)
            result = self.Triggering_.isExternalTriggerMarkedForDeletion() ;
        end  % function
        
        function result = triggerNames(self)
            result = self.Triggering_.triggerNames() ; 
        end  % function        
        
        function result = acquisitionTriggerProperty(self, propertyName)
            result = self.Triggering_.acquisitionTriggerProperty(propertyName) ;
        end  % function        
            
        function result = stimulationTriggerProperty(self, propertyName)
            result = self.Triggering_.stimulationTriggerProperty(propertyName) ;
        end  % function        
            
        function result = counterTriggerProperty(self, index, propertyName)
            result = self.Triggering_.counterTriggerProperty(index, propertyName) ;
        end  % function
        
        function result = externalTriggerProperty(self, index, propertyName)
            result = self.Triggering_.externalTriggerProperty(index, propertyName) ;
        end  % function
    end  % public methods block
    
    methods  % expose stim library methods at WSM level
        function clearStimulusLibrary(self)
            try
                self.Stimulation_.clearStimulusLibrary() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                self.broadcast('UpdateMain');
                self.broadcast('UpdateGeneral') ;
                self.broadcast('UpdateChannels') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setSelectedStimulusLibraryItemByClassNameAndIndex(self, className, index)
            try
                self.Stimulation_.setSelectedStimulusLibraryItemByClassNameAndIndex(className, index) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function index = addNewStimulusSequence(self)
            try
                index = self.Stimulation_.addNewStimulusSequence() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateMain');
            self.broadcast('UpdateGeneral') ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function duplicateSelectedStimulusLibraryItem(self)
            try
                self.Stimulation_.duplicateSelectedStimulusLibraryItem() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function bindingIndex = addBindingToSelectedStimulusLibraryItem(self)
            try
                bindingIndex = self.Stimulation_.addBindingToSelectedStimulusLibraryItem() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
                
        function bindingIndex = addBindingToStimulusLibraryItem(self, className, itemIndex)
            try
                bindingIndex = self.Stimulation_.addBindingToStimulusLibraryItem(className, itemIndex) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
                
        function deleteMarkedBindingsFromSequence(self)
            try
                self.Stimulation_.deleteMarkedBindingsFromSequence() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function index = addNewStimulusMap(self)
            try
                index = self.Stimulation_.addNewStimulusMap() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateMain');
        end  % function        
        
%         function addChannelToSelectedStimulusLibraryItem(self)
%             try
%                 self.Stimulation_.addChannelToSelectedStimulusLibraryItem() ;
%             catch exception
%                 self.broadcast('UpdateStimulusLibrary') ;
%                 rethrow(exception) ;
%             end
%             self.broadcast('UpdateStimulusLibrary') ;                
%         end  % function
                
        function deleteMarkedChannelsFromSelectedStimulusLibraryItem(self)
            try
                self.Stimulation_.deleteMarkedChannelsFromSelectedStimulusLibraryItem() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function stimulusIndex = addNewStimulus(self)
            try
                stimulusIndex = self.Stimulation_.addNewStimulus() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function        
        
        function result = isSelectedStimulusLibraryItemInUse(self)
            result = self.Stimulation_.isSelectedStimulusLibraryItemInUse() ;
        end  % function        

        function deleteSelectedStimulusLibraryItem(self)
            try
                self.Stimulation_.deleteSelectedStimulusLibraryItem() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                self.broadcast('UpdateGeneral') ;  % Need to update the list of outputables, maybe
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateGeneral') ;  % Need to update the list of outputables, maybe
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function        
        
        function result = selectedStimulusLibraryItemClassName(self)
            result = self.Stimulation_.selectedStimulusLibraryItemClassName() ;
        end  % function        

        function result = selectedStimulusLibraryItemIndexWithinClass(self)
            result = self.Stimulation_.selectedStimulusLibraryItemIndexWithinClass() ;
        end  % function        

        function setSelectedStimulusLibraryItemProperty(self, propertyName, newValue)
            try
                didSetOutputableName = self.Stimulation_.setSelectedStimulusLibraryItemProperty(propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            if didSetOutputableName ,
                self.broadcast('UpdateGeneral');
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function       
        
        function setSelectedStimulusAdditionalParameter(self, iParameter, newString)
            try
                self.Stimulation_.setSelectedStimulusAdditionalParameter(iParameter, newString) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function                
        
        function setBindingOfSelectedSequenceToNamedMap(self, indexOfElementWithinSequence, newMapName)
            try
                self.Stimulation_.setBindingOfSelectedSequenceToNamedMap(indexOfElementWithinSequence, newMapName) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function                
            
%         function setIsMarkedForDeletionForElementOfSelectedSequence(self, indexOfElementWithinSequence, newValue)
%             try
%                 self.Stimulation_.setIsMarkedForDeletionForElementOfSelectedSequence(indexOfElementWithinSequence, newValue) ;
%             catch exception
%                 self.broadcast('UpdateStimulusLibrary') ;
%                 rethrow(exception) ;
%             end
%             self.broadcast('UpdateStimulusLibrary') ;                
%         end  % function        

        function setBindingOfSelectedMapToNamedStimulus(self, indexOfBindingWithinMap, newStimulusName)
            try
                self.Stimulation_.setBindingOfSelectedMapToNamedStimulus(indexOfBindingWithinMap, newStimulusName) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function                
        
%         function result = propertyForElementOfSelectedStimulusLibraryItem(self, indexOfElementWithinItem, propertyName)
%             result = self.Stimulation_.propertyForElementOfSelectedStimulusLibraryItem(indexOfElementWithinItem, propertyName) ;
%         end  % function        
        
        function setPropertyForElementOfSelectedMap(self, indexOfElementWithinMap, propertyName, newValue)
            try
                self.Stimulation_.setPropertyForElementOfSelectedMap(indexOfElementWithinMap, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function        
        
%         function plotSelectedStimulusLibraryItem(self, figureGH)
%             sampleRate = self.StimulationSampleRate ;  % Hz 
%             channelNames = [self.AOChannelNames self.DOChannelNames] ;
%             isChannelAnalog = [true(size(self.AOChannelNames)) false(size(self.DOChannelNames))] ;
%             self.Stimulation_.plotSelectedStimulusLibraryItem(figureGH, sampleRate, channelNames, isChannelAnalog) ;
%         end  % function
        
        function [y, t] = previewStimulus(self, stimulusIndex)
            sampleRate = self.StimulationSampleRate ;  % Hz 
            [y, t] = self.Stimulation_.previewStimulus(stimulusIndex, sampleRate) ;
        end
        
        function [y, t] = previewStimulusMap(self, mapIndex)
            sampleRate = self.StimulationSampleRate ;  % Hz 
            channelNames = [self.AOChannelNames self.DOChannelNames] ;
            isChannelAnalog = [true(size(self.AOChannelNames)) false(size(self.DOChannelNames))] ;
            [y, t] = self.Stimulation_.previewStimulusMap(mapIndex, sampleRate, channelNames, isChannelAnalog) ;
        end
        
        function result = selectedStimulusLibraryItemProperty(self, propertyName)
            result = self.Stimulation_.selectedStimulusLibraryItemProperty(propertyName) ;
        end  % method        
        
        function result = propertyFromEachStimulusLibraryItemInClass(self, className, propertyName) 
            % Result is a cell array, even it seems like it could/should be another kind of array
            result = self.Stimulation_.propertyFromEachStimulusLibraryItemInClass(className, propertyName) ;
        end  % function
        
        function result = indexOfStimulusLibraryClassSelection(self, className)
            result = self.Stimulation_.indexOfStimulusLibraryClassSelection(className) ;
        end  % method                    

        function result = stimulusLibraryClassSelectionProperty(self, className, propertyName)
            result = self.Stimulation_.stimulusLibraryClassSelectionProperty(className, propertyName) ;
        end  % method        
        
        function result = stimulusLibrarySelectedItemProperty(self, propertyName)
            result = self.Stimulation_.stimulusLibrarySelectedItemProperty(propertyName) ;
        end

        function result = stimulusLibrarySelectedItemBindingProperty(self, bindingIndex, propertyName)
            result = self.Stimulation_.stimulusLibrarySelectedItemBindingProperty(bindingIndex, propertyName) ;
        end
        
        function result = stimulusLibrarySelectedItemBindingTargetProperty(self, bindingIndex, propertyName)
            result = self.Stimulation_.stimulusLibrarySelectedItemBindingTargetProperty(bindingIndex, propertyName) ;
        end
        
        function result = stimulusLibraryItemProperty(self, className, index, propertyName)
            result = self.Stimulation_.stimulusLibraryItemProperty(className, index, propertyName) ;
        end  % function                        
        
        function result = stimulusLibraryItemBindingProperty(self, className, itemIndex, bindingIndex, propertyName)
            result = self.Stimulation_.stimulusLibraryItemBindingProperty(className, itemIndex, bindingIndex, propertyName) ;
        end  % function                        
        
        function result = stimulusLibraryItemBindingTargetProperty(self, className, itemIndex, bindingIndex, propertyName)
            result = self.Stimulation_.stimulusLibraryItemBindingTargetProperty(className, itemIndex, bindingIndex, propertyName) ;
        end  % function                                
        
        function result = isStimulusLibraryItemBindingTargetEmpty(self, className, itemIndex, bindingIndex)
            result = self.Stimulation_.isStimulusLibraryItemBindingTargetEmpty(className, itemIndex, bindingIndex) ;
        end  % function
        
        function result = isStimulusLibrarySelectedItemBindingTargetEmpty(self, bindingIndex)
            result = self.Stimulation_.isStimulusLibrarySelectedItemBindingTargetEmpty(bindingIndex) ;
        end  % function        
        
        function result = isStimulusLibraryEmpty(self)
            result = self.Stimulation_.isStimulusLibraryEmpty() ;
        end  % function        
        
        function result = isAStimulusLibraryItemSelected(self)
            result = self.Stimulation_.isAStimulusLibraryItemSelected() ;
        end  % function        

        function result = isAnyBindingMarkedForDeletionForStimulusLibrarySelectedItem(self)
            result = self.Stimulation_.isAnyBindingMarkedForDeletionForStimulusLibrarySelectedItem() ;            
        end  % function        

        function result = stimulusLibraryOutputableNames(self)
            result = self.Stimulation_.stimulusLibraryOutputableNames() ;            
        end
        
        function result = stimulusLibrarySelectedOutputableProperty(self, propertyName)
            result = self.Stimulation_.stimulusLibrarySelectedOutputableProperty(propertyName) ;            
        end
        
        function result = areStimulusLibraryMapDurationsOverridden(self)
            result = self.Stimulation_.areStimulusLibraryMapDurationsOverridden() ;            
        end
        
        function setStimulusLibraryItemProperty(self, className, index, propertyName, newValue)
            try
                didSetOutputableName = self.Stimulation_.setStimulusLibraryItemProperty(className, index, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            % Special case to deal with renaming an outputable
            if didSetOutputableName ,
                self.broadcast('UpdateGeneral') ;
            else
                self.broadcast('DidMaybeChangeProtocol') ;                
            end
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('UpdateStimulusPreview') ;
        end        
        
        function setStimulusLibraryItemBindingProperty(self, className, itemIndex, bindingIndex, propertyName, newValue)
            try
                self.Stimulation_.setStimulusLibraryItemBindingProperty(className, itemIndex, bindingIndex, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function                        
        
        function result = isStimulusLibraryItemInUse(self, className, itemIndex)
            result = self.Stimulation_.isStimulusLibraryItemInUse(className, itemIndex) ;
        end
        
        function deleteStimulusLibraryItem(self, className, itemIndex)
            try
                self.Stimulation_.deleteStimulusLibraryItem(className, itemIndex) ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                self.broadcast('UpdateGeneral') ;  % Need to update the list of outputables, maybe
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateGeneral') ;  % Need to update the list of outputables, maybe
            self.broadcast('DidMaybeChangeProtocol') ;
            
        end
        
        function setSelectedStimulusLibraryItemWithinClassBindingProperty(self, className, bindingIndex, propertyName, newValue)
            try
                self.Stimulation_.setSelectedStimulusLibraryItemWithinClassBindingProperty(className, bindingIndex, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || ~isequal(propertyName, 'IsMarkedForDeletion') ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function        
        
        function populateStimulusLibraryForTesting(self)
            try
                self.Stimulation_.populateStimulusLibraryForTesting() ;
                self.DoesProtocolNeedSave_ = true ;
            catch exception
                self.broadcast('UpdateStimulusLibrary') ;
                self.broadcast('UpdateStimulusPreview') ;
                self.broadcast('UpdateMain') ;
                rethrow(exception) ;
            end
            self.broadcast('UpdateStimulusLibrary') ;                
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('UpdateMain') ;
        end  % function        
        
        function result = getCommandServer_(self)
            % This is intended to be used for testing and debugging only.
            result = self.CommandServer_ ;
        end
        
        function result = getCommandClient_(self)
            % This is intended to be used for testing and debugging only.
            result = self.CommandClient_ ;
        end        
    end  % public methods block
    
    methods
        function value = get.AIChannelScales(self)
            ephys = self.Ephys_ ;
            %electrodeManager = ephys.ElectrodeManager ;
            aiChannelNames = self.Acquisition_.AnalogChannelNames ;
            [channelScalesFromElectrodes, isChannelScaleEnslaved] = ephys.getMonitorScalingsByName(aiChannelNames);
            value = ws.fif(isChannelScaleEnslaved, channelScalesFromElectrodes, self.Acquisition_.getAnalogChannelScales_());
        end
        
        function value = getNumberOfElectrodesClaimingAIChannel(self)
            ephys = self.Ephys_ ;
            %electrodeManager = ephys.ElectrodeManager ;
            channelNames = self.Acquisition_.AnalogChannelNames ;
            value = ephys.getNumberOfElectrodesClaimingMonitorChannel(channelNames) ;
        end
        
        function value = get.AIChannelUnits(self)            
            ephys=self.Ephys_;
            %electrodeManager=ephys.ElectrodeManager;
            channelNames=self.Acquisition_.AnalogChannelNames;
            [channelUnitsFromElectrodes, isChannelScaleEnslaved] = ephys.getMonitorUnitsByName(channelNames);
            value = ws.fif(isChannelScaleEnslaved, channelUnitsFromElectrodes, self.Acquisition_.getAnalogChannelUnits_());
        end
        
        function set.AIChannelUnits(self,newValue)
            isChangeable= ~(self.getNumberOfElectrodesClaimingAIChannel()==1);
            self.Acquisition_.setAnalogChannelUnits_(ws.fif(isChangeable, newValue, self.Acquisition_.getAnalogChannelUnits_())) ;
            self.DoesProtocolNeedSave_ = true ;
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function set.AIChannelScales(self,newValue)
            isChangeable= ~(self.getNumberOfElectrodesClaimingAIChannel()==1);
            self.Acquisition_.setAnalogChannelScales_(ws.fif(isChangeable, newValue, self.Acquisition_.getAnalogChannelScales_())) ;
            self.DoesProtocolNeedSave_ = true ;
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setAIChannelUnitsAndScales(self,newUnitsRaw,newScales)
            isChangeable= ~(self.getNumberOfElectrodesClaimingAIChannel()==1);
            newUnits = cellfun(@strtrim,newUnitsRaw,'UniformOutput',false) ;
            self.Acquisition_.setAnalogChannelUnits_( ws.fif(isChangeable, newUnits, self.Acquisition_.getAnalogChannelUnits_()) ) ;
            self.Acquisition_.setAnalogChannelScales_( ws.fif(isChangeable, newScales, self.Acquisition_.getAnalogChannelScales_()) ) ;
            self.DoesProtocolNeedSave_ = true ;
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setSingleAIChannelUnits(self,i,newValueRaw)
            isChangeableFull = ~(self.getNumberOfElectrodesClaimingAIChannel()==1) ;
            isChangeable = isChangeableFull(i) ;
            if isChangeable , 
                newValue = strtrim(newValueRaw) ;
                self.Acquisition_.setSingleAnalogChannelUnits_(i, newValue) ;                
                self.DoesProtocolNeedSave_ = true ;
            end
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setSingleAIChannelScale(self,i,newValue)
            isChangeableFull = ~(self.getNumberOfElectrodesClaimingAIChannel()==1);
            isChangeable = isChangeableFull(i);
            if isChangeable ,
                self.Acquisition_.setSingleAnalogChannelScale_(i, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            end
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function result = get.IsAIChannelActive(self)
            result = self.Acquisition_.IsAnalogChannelActive ;
        end
        
        function setSingleIsAIChannelActive(self, aiChannelIndex, newValue)
            % Boolean array indicating which of the AI channels is
            % active.
            self.Acquisition_.setSingleIsAnalogChannelActive(aiChannelIndex, newValue) ;
            self.syncTraces_() ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('EMDidSetIsInputChannelActive') ;
            self.broadcast('TPDidSetIsInputChannelActive') ;
            self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateChannels') ;
        end    

        function result = getTimestampsForDataInCache(self)
            dt = 1/self.AcquisitionSampleRate ;
            tPastLast = self.t_ ;  % the timestamp of the scan that would follow the last scan in the cache
            n = self.Acquisition_.getNScansInCache() ;
            result = dt*(0:(n-1))' - dt*n + tPastLast ;
        end
        
        function scaledAnalogData = getLatestAIData(self)
            % Get the data from the most-recent data available callback, as
            % doubles.
            rawAnalogData = self.Acquisition_.getLatestRawAnalogData() ;
            channelScales = self.AIChannelScales(self.IsAIChannelActive) ;
            scalingCoefficients = self.Acquisition_.AnalogScalingCoefficients ;
            scaledAnalogData = ws.scaledDoubleAnalogDataFromRawMex(rawAnalogData, channelScales, scalingCoefficients) ;
        end  % function
        
        function scaledAnalogData = getAIDataFromCache(self)
            % Get the data from the main-memory cache, as double-precision floats.  This
            % call unwraps the circular buffer for you.
            rawAnalogData = self.Acquisition_.getRawAnalogDataFromCache();
            if isempty(rawAnalogData) ,
                scaledAnalogData = double(rawAnalogData) ;
            else
                %channelScales=self.AIChannelScales(self.IsAIChannelActive);
                channelScales = self.AIChannelScales(self.Acquisition_.IsInCacheFromAnalogChannelIndex);
                scalingCoefficients = self.Acquisition_.AnalogScalingCoefficients ;
                scaledAnalogData = ws.scaledDoubleAnalogDataFromRawMex(rawAnalogData, channelScales, scalingCoefficients) ;
            end
        end  % function

%         function scaledData = getSinglePrecisionAIDataFromCache(self)
%             % Get the data from the main-memory cache, as single-precision floats.  This
%             % call unwraps the circular buffer for you.
%             rawAnalogData = self.Acquisition_.getRawAnalogDataFromCache();
%             channelScales=self.AIChannelScales(self.IsAIChannelActive);
%             scalingCoefficients = self.Acquisition_.AnalogScalingCoefficients ;
%             scaledData = ws.scaledSingleAnalogDataFromRaw(rawAnalogData, channelScales, scalingCoefficients) ;
%         end  % function

        function [result, signalCount] = getDIDataFromCache(self)
            % Get the data from the main-memory cache, as double-precision floats.  This
            % call unwraps the circular buffer for you.
            [result, signalCount] = self.Acquisition_.getRawDigitalDataFromCache();
        end  % function
        
        function result = get.IsDIChannelActive(self)
            result = self.Acquisition_.IsDigitalChannelActive ;
        end

        function setSingleIsDIChannelActive(self, diChannelIndex, newValue)
            % Boolean array indicating which of the AI channels is
            % active.
            self.Acquisition_.setSingleIsDigitalChannelActive(diChannelIndex, newValue) ;
            self.syncTraces_() ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('EMDidSetIsInputChannelActive') ;
            self.broadcast('TPDidSetIsInputChannelActive') ;
            self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;                        
            self.broadcast('UpdateChannels') ;
        end    
        
        function result = getLatestDIData(self)
            result = self.Acquisition_.getLatestRawDigitalData() ;
        end  % function
        
        function result = aiChannelUnitsFromName(self, channelName)
            if isempty(channelName) ,
                result='';
            else
                iChannel=self.Acquisition_.aiChannelIndexFromName(channelName) ;
                if isempty(iChannel) || isnan(iChannel) ,
                    result='';
                else
                    result=self.AIChannelUnits{iChannel} ;
                end
            end
        end
        
        function result = aiChannelScaleFromName(self, channelName)
            if isempty(channelName) ,
                result='';
            else
                iChannel=self.Acquisition_.aiChannelIndexFromName(channelName);
                if isempty(iChannel) || isnan(iChannel) ,
                    result='';
                else
                    result=self.AIChannelScales(iChannel);
                end
            end
        end  % function
        
        function result=get.IsAIChannelMarkedForDeletion(self)
            % Boolean array indicating which of the available AI channels is
            % active.
            result = self.Acquisition_.getIsAnalogChannelMarkedForDeletion_() ;
        end
        
        function set.IsAIChannelMarkedForDeletion(self,newValue)
            % Boolean array indicating which of the AI channels is
            % active.
            self.Acquisition_.setIsAnalogChannelMarkedForDeletion_(newValue) ;
            self.broadcast('UpdateChannels') ;
        end
        
        function result=get.IsDIChannelMarkedForDeletion(self)
            % Boolean array indicating which of the available AI channels is
            % active.
            result = self.Acquisition_.getIsDigitalChannelMarkedForDeletion_() ;
        end
        
        function set.IsDIChannelMarkedForDeletion(self,newValue)
            % Boolean array indicating which of the AI channels is
            % active.
            self.Acquisition_.setIsDigitalChannelMarkedForDeletion_(newValue) ;
            self.broadcast('UpdateChannels') ;
        end
        
        function setSingleAIChannelName(self, i, newValue)
            allChannelNames = self.AllChannelNames ;
            [didSucceed, oldValue] = self.Acquisition_.setSingleAnalogChannelName_(i, newValue, allChannelNames) ;
            if didSucceed, 
                self.DoesProtocolNeedSave_ = true ;
            end
            display=self.Display_;
            if ~isempty(display)
                %display.didSetAnalogInputChannelName(didSucceed, oldValue, newValue);
                %self.broadcast('UpdateTraces') ;
                self.broadcast('UpdateMain') ;
            end            
            ephys=self.Ephys_;
            if ~isempty(ephys)
                ephys.didSetAnalogInputChannelName(didSucceed, oldValue, newValue);
                self.broadcast('UpdateElectrodeManager') ;
            end            
            self.broadcast('UpdateChannels') ;
        end
        
        function setSingleDIChannelName(self, i, newValue)
            allChannelNames = self.AllChannelNames ;
            didSucceed = self.Acquisition_.setSingleDigitalChannelName_(i, newValue, allChannelNames) ;
            if didSucceed, 
                self.DoesProtocolNeedSave_ = true ;
            end
            %self.Display_.didSetDigitalInputChannelName(didSucceed, oldValue, newValue);
            %self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateChannels') ;
        end
        
        function set.AcquisitionSampleRate(self, newValue)
            if isscalar(newValue) && isnumeric(newValue) && isfinite(newValue) && newValue>0 ,                
                % Constrain value appropriately
                isValueValid = true ;
                newValue = double(newValue) ;
                sampleRate = ws.WavesurferModel.coerceSampleFrequencyToAllowedValue(self.PrimaryDeviceName, self.IsPrimaryDeviceAPXIDevice, newValue) ;
                self.Acquisition_.setSampleRate_(sampleRate) ;
                self.Ephys_.didSetAcquisitionSampleRate(newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            else
                isValueValid = false ;
            end
            self.broadcast('DidSetAcquisitionSampleRate');
            self.broadcast('DidMaybeChangeProtocol');
            if ~isValueValid ,
                error('ws:invalidPropertyValue', ...
                      'AcquisitionSampleRate must be a positive finite numeric scalar');
            end                
        end  % function
        
        function out = get.AcquisitionSampleRate(self)
            out = self.Acquisition_.getSampleRate_() ;
        end  % function
        
        function out = get.ExpectedSweepScanCount(self)            
            out = ws.nScansFromScanRateAndDesiredDuration(self.AcquisitionSampleRate, self.SweepDuration) ;
        end  % function
        
        function value = get.IsXSpanSlavedToAcquistionDuration(self)
            if self.AreSweepsContinuous ,
                value = false ;
            else
                value = self.Display_.IsXSpanSlavedToAcquistionDuration ;
            end
        end  % function
        
        function set.IsXSpanSlavedToAcquistionDuration(self, newValue)
            if self.IsXSpanSlavedToAcquistionDurationSettable ,
                if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && isfinite(newValue))) ,
                    isNewValueAllowed = true ;
                    self.Display_.IsXSpanSlavedToAcquistionDuration = logical(newValue) ;
                    self.syncTraces_() ;
                    self.broadcast('UpdateTraces');
                else
                    isNewValueAllowed = false ;
                end
            else
                isNewValueAllowed = true ;  % sort of in a trivial sense...
            end
            self.broadcast('UpdateMain');            
            self.broadcast('UpdateGeneral');            
            if ~isNewValueAllowed ,
                error('ws:invalidPropertyValue', ...
                      'IsXSpanSlavedToAcquistionDuration must be a logical scalar, or convertible to one') ;
            end                            
            
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol');
        end
        
        function value = get.IsXSpanSlavedToAcquistionDurationSettable(self)
            value = self.AreSweepsFiniteDuration ;
        end  % function             
        
        function value = get.XSpan(self)
            if self.IsXSpanSlavedToAcquistionDuration ,
                sweepDuration = self.SweepDuration ;
                value = ws.fif(isfinite(sweepDuration), sweepDuration, 1) ;
            else
                value = self.Display_.getXSpan() ;
            end
        end
        
        function set.XSpan(self, newValue)            
            if self.IsXSpanSlavedToAcquistionDuration ,
                % don't set anything
                didSucceed = true ;  % this is by convention
            else
                if isnumeric(newValue) && isscalar(newValue) && isfinite(newValue) && newValue>0 ,
                    self.Display_.setXSpan(double(newValue)) ;
                    self.syncTraces_() ;
                    self.broadcast('UpdateTraces') ;
                    didSucceed = true ;
                else
                    didSucceed = false ;
                end
            end
            self.broadcast('DidSetXSpan');
            if ~didSucceed ,
                error('ws:invalidPropertyValue', ...
                      'XSpan must be a scalar finite positive number') ;
            end                                      
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function toggleIsAIChannelDisplayed(self, aiChannelIndex) 
            nAIChannels = self.NAIChannels ;
            try
                self.Display_.toggleIsAnalogChannelDisplayed(aiChannelIndex, nAIChannels) ;
            catch err
                self.broadcast('UpdateMain') ;
                rethrow(err) ;
            end
            self.broadcast('UpdateMain') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function toggleIsDIChannelDisplayed(self, diChannelIndex) 
            nDIChannels = self.NDIChannels ;
            try
                self.Display_.toggleIsDigitalChannelDisplayed(diChannelIndex, nDIChannels) ;
            catch err
                self.broadcast('UpdateMain') ;
                rethrow(err) ;
            end
            self.broadcast('UpdateMain') ;            
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function result = get.IsAIChannelDisplayed(self)
            result = self.Display_.IsAnalogChannelDisplayed ;
        end
        
        function result = get.IsDIChannelDisplayed(self)
            result = self.Display_.IsDigitalChannelDisplayed ;
        end
        
        function result = get.NAIChannels(self)
            result = self.Acquisition_.NAnalogChannels ;
        end

        function result = get.NDIChannels(self)
            result = self.Acquisition_.NDigitalChannels ;
        end
        
        function result = get.NAOChannels(self)
            result = self.Stimulation_.NAnalogChannels ;
        end

        function result = get.NDOChannels(self)
            result = self.Stimulation_.NDigitalChannels ;
        end
        
        function result=getNumberOfElectrodesClaimingAOChannel(self)
            ephys=self.Ephys_;
            %electrodeManager=ephys.ElectrodeManager;
            channelNames=self.Stimulation_.AnalogChannelNames;
            result = ephys.getNumberOfElectrodesClaimingCommandChannel(channelNames);
        end  % function       
        
        function result = get.AOChannelScales(self)
            ephys=self.Ephys_;
            %electrodeManager=ephys.ElectrodeManager;
            channelNames = self.Stimulation_.AnalogChannelNames ;
            [analogChannelScalesFromElectrodes, isChannelScaleEnslaved] = ephys.getCommandScalingsByName(channelNames) ;
            result = ws.fif(isChannelScaleEnslaved,analogChannelScalesFromElectrodes,self.Stimulation_.getAnalogChannelScales_()) ;
        end  % function
        
%         function set.AOChannelScales(self, newValue)
%             oldValue = self.Stimulation_.getAnalogChannelScales_() ;
%             isChangeable = ~(self.getNumberOfElectrodesClaimingAOChannel()==1) ;
%             editedNewValue = ws.fif(isChangeable,newValue,oldValue) ;
%             self.Stimulation_.setAnalogChannelScales_(editedNewValue) ;
%             self.didSetAnalogChannelUnitsOrScales() ;            
%         end  % function        
        
        function result = get.AOChannelUnits(self)
            ephys = self.Ephys_ ;
            %electrodeManager = ephys.ElectrodeManager ;
            channelNames = self.Stimulation_.AnalogChannelNames ;
            [channelUnitsFromElectrodes, isChannelScaleEnslaved] = ephys.getCommandUnitsByName(channelNames) ;
            result = ws.fif(isChannelScaleEnslaved, channelUnitsFromElectrodes, self.Stimulation_.getAnalogChannelUnits_()) ;
        end  % function

        function setSingleAOChannelUnits(self, i, newValue)
            isChangeableFull=~(self.getNumberOfElectrodesClaimingAOChannel()==1);
            isChangeable= isChangeableFull(i);
            if isChangeable ,
                self.Stimulation_.setSingleAnalogChannelUnits_(i, strtrim(newValue)) ;
                self.DoesProtocolNeedSave_ = true ;
            end
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_();            
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setSingleAOChannelScale(self, i, newValue)
            isChangeableFull = ~(self.getNumberOfElectrodesClaimingAOChannel()==1) ;
            isChangeable = isChangeableFull(i) ;
            if isChangeable && isfinite(newValue) && newValue>0 ,
                self.Stimulation_.setSingleAnalogChannelScale_(i, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
            end
            self.notifyOtherSubsystemsThatDidSetAnalogChannelUnitsOrScales_() ;
            self.broadcast('UpdateChannels') ;
        end  % function
        
        function setSingleAOChannelTerminalID(self, i, newValue)
            if 1<=i && i<=self.NAOChannels && isnumeric(newValue) && isscalar(newValue) && isfinite(newValue) ,
                newValueAsDouble = double(newValue) ;
                if newValueAsDouble>=0 && newValueAsDouble==round(newValueAsDouble) ,
                    self.Stimulation_.setSingleAnalogTerminalID(i, newValueAsDouble) ;
                    self.DoesProtocolNeedSave_ = true ;
                end
            end
            self.syncIsAOChannelTerminalOvercommitted_() ;
            self.broadcast('UpdateChannels') ;
        end
        
        function set.IsDOChannelTimed(self, newValue)
            try
                wasSet = self.Stimulation_.setIsDigitalChannelTimed_(newValue) ;
            catch exception
                self.broadcast('EMDidSetIsDigitalOutputTimed') ;
                self.broadcast('UpdateChannels') ;
                rethrow(exception) ;
            end            
            self.broadcast('EMDidSetIsDigitalOutputTimed') ;
            self.broadcast('UpdateChannels') ;            
            if wasSet ,
                self.isDigitalChannelTimedWasSetInStimulationSubsystem() ;
                self.DoesProtocolNeedSave_ = true ;
            end
        end
        
        function result = get.IsDOChannelTimed(self) 
            result = self.Stimulation_.getIsDigitalChannelTimed_() ;
        end
        
        function set.DOChannelStateIfUntimed(self, newValue)
            try
                self.Stimulation_.setDigitalOutputStateIfUntimed_(newValue) ;
            catch exception
                self.broadcast('UpdateDigitalOutputStateIfUntimed') ;
                rethrow(exception) ;
            end
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateDigitalOutputStateIfUntimed') ;
            self.broadcast('DidMaybeChangeProtocol') ;            
            self.Looper_.digitalOutputStateIfUntimedWasSetInFrontend(self.DOChannelStateIfUntimed, self.IsDOChannelTimed) ;
            self.Refiller_.digitalOutputStateIfUntimedWasSetInFrontend(self.DOChannelStateIfUntimed) ;
        end  % function
        
        function out = get.DOChannelStateIfUntimed(self)
            out= self.Stimulation_.getDigitalOutputStateIfUntimed_() ;
        end
        
        function out = get.StimulationSampleRate(self)
            out= self.Stimulation_.getSampleRate_() ;
        end
        
        function set.StimulationSampleRate(self, newValue)
            if isscalar(newValue) && isnumeric(newValue) && isfinite(newValue) && newValue>0 ,                
                % Constrain value appropriately
                isValueValid = true ;
                newValue = double(newValue) ;
                sampleRate = ws.WavesurferModel.coerceSampleFrequencyToAllowedValue(self.PrimaryDeviceName, self.IsPrimaryDeviceAPXIDevice, newValue) ;
                self.Stimulation_.setSampleRate_(sampleRate) ;
                self.DoesProtocolNeedSave_ = true ;
                %self.didSetAcquisitionSampleRate(sampleRate);
            else
                isValueValid = false ;
            end
            self.broadcast('DidSetStimulationSampleRate');
            self.broadcast('DidMaybeChangeProtocol') ;            
            if ~isValueValid ,
                error('ws:invalidPropertyValue', ...
                      'StimulationSampleRate must be a positive finite numeric scalar');
            end                
        end  % function
        
        function result=get.IsAOChannelMarkedForDeletion(self)
            result =  self.Stimulation_.getIsAnalogChannelMarkedForDeletion_() ;
        end
        
        function set.IsAOChannelMarkedForDeletion(self, newValue)
            self.Stimulation_.setIsAnalogChannelMarkedForDeletion_(newValue) ;
            self.broadcast('UpdateChannels') ;
        end
        
        function result=get.IsDOChannelMarkedForDeletion(self)
            result = self.Stimulation_.getIsDigitalChannelMarkedForDeletion_() ;
        end
        
        function set.IsDOChannelMarkedForDeletion(self, newValue)
            self.Stimulation_.setIsDigitalChannelMarkedForDeletion_(newValue) ;
            self.broadcast('UpdateChannels') ;
        end
        
        function result=aoChannelUnitsFromName(self,channelName)
            if isempty(channelName) ,
                result = '' ;
            else
                iChannel=self.Stimulation_.aoChannelIndexFromName(channelName);
                if isnan(iChannel) ,
                    result='';
                else
                    result=self.AOChannelUnits{iChannel} ;
                end
            end
        end  % function
        
        function value = aoChannelScaleFromName(self, channelName)
            channelIndex = self.Stimulation_.aoChannelIndexFromName(channelName) ;
            if isnan(channelIndex) ,
                value = nan ;
            else
                value = self.AOChannelScales(channelIndex) ;
            end
        end  % function

        function result = isTestPulsingEnabled(self)
            electrodeIndex = self.TestPulseElectrodeIndex ;
            result = ...
                (self.isIdle() || self.isTestPulsing()) && ...
                ~isempty(electrodeIndex) && ...
                self.areTestPulseElectrodeChannelsValid() && ...
                ~self.areTestPulseElectrodeMonitorAndCommandChannelsOnDiffrentDevices() ;            
%                 self.areAllMonitorAndCommandChannelNamesDistinct() && ...
        end
        
        function startTestPulsing(self)
            if self.isTestPulsing() ,
                return
            end            
            
            try       
                % Takes some time to start...
                self.changeReadiness_(-1) ;

                % Update the smart electrode channel scales, if possible and
                % needed
                if self.DoTrodeUpdateBeforeRun ,
                    self.updateSmartElectrodeGainsAndModes();
                end

                % Check that we can start, and if not, return
                %electrodeIndex = self.TestPulseElectrodeIndex ;
                canStart = self.isTestPulsingEnabled() ;
                if ~canStart ,
                    return
                end

                % Free up resources we will need for test pulsing
                self.releaseTimedHardwareResourcesOfAllProcesses_();
                
                % Get things we need for test-pulsing
                fs = self.AcquisitionSampleRate ;
                isVCPerTestPulseElectrode = self.getIsVCPerTestPulseElectrode() ;
                isCCPerTestPulseElectrode = self.getIsCCPerTestPulseElectrode() ;            
                commandDeviceNamePerTestPulseElectrode = self.getCommandDeviceNamePerTestPulseElectrode_() ;
                monitorDeviceNamePerTestPulseElectrode = self.getMonitorDeviceNamePerTestPulseElectrode_() ;
                commandTerminalIDPerTestPulseElectrode = self.getCommandTerminalIDPerTestPulseElectrode_() ;
                monitorTerminalIDPerTestPulseElectrode = self.getMonitorTerminalIDPerTestPulseElectrode_() ;
                commandChannelScalePerTestPulseElectrode = self.getCommandChannelScalePerTestPulseElectrode_() ;
                monitorChannelScalePerTestPulseElectrode = self.getMonitorChannelScalePerTestPulseElectrode_() ;
                primaryDeviceName =self.PrimaryDeviceName ;
                isPrimaryDeviceAPXIDevice = ws.isDeviceAPXIDevice(primaryDeviceName) ;
                gainOrResistanceUnitsPerTestPulseElectrode = self.getGainOrResistanceUnitsPerTestPulseElectrode() ;
                
                % Make sure all the channels are on the same device
                allDeviceNames = unique([commandDeviceNamePerTestPulseElectrode monitorDeviceNamePerTestPulseElectrode]) ;
                if ~isscalar(allDeviceNames) ,
                    error('ws:allTestPulseChannelsMustBeOnSameDevice', ...
                          'All test pulse channels must be on the same device') ;
                end
                deviceName = allDeviceNames{1} ;
                
                % Call the main routine
                self.Ephys_.prepareForTestPulsing(fs, ...
                                                  isVCPerTestPulseElectrode, ...
                                                  isCCPerTestPulseElectrode, ...
                                                  commandTerminalIDPerTestPulseElectrode, ...
                                                  monitorTerminalIDPerTestPulseElectrode, ...
                                                  commandChannelScalePerTestPulseElectrode, ...
                                                  monitorChannelScalePerTestPulseElectrode, ...
                                                  deviceName, ...
                                                  primaryDeviceName, ...
                                                  isPrimaryDeviceAPXIDevice, ...
                                                  gainOrResistanceUnitsPerTestPulseElectrode, ...
                                                  self) ;

                % Change our state
                if isequal(self.State,'idle') ,
                    self.setState_('test_pulsing');
                end

                % OK, now we consider the TP no longer busy
                self.changeReadiness_(+1);

                % Actually start the test pulsing
                self.Ephys_.startTestPulsing() ;                
            catch exception
                self.abortTestPulsing_() ;
                self.changeReadiness_(+1) ;
                rethrow(exception) ;                
            end
        end

        function stopTestPulsing(self)
            if ~self.isTestPulsing() ,
                %fprintf('About to exit stop() via short-circuit...\n');                            
                return
            end
            
            try             
                self.changeReadiness_(-1) ;  % Takes some time to stop
                %fprintf('About to call self.Ephys_.stopTestPulsing()\n');
                self.Ephys_.stopTestPulsing() ;
                %fprintf('Done with call to self.Ephys_.stopTestPulsing()\n');
                self.changeReadiness_(+1) ;
                if isequal(self.State,'test_pulsing') ,
                    self.setState_('idle');
                end
            catch exception
                fprintf('Hit an exception while trying to stop test pulsing\n');
                self.abortTestPulsing_() ;
                self.changeReadiness_(+1) ;
                rethrow(exception) ;                                
            end            
        end  % function    
        
        function completingTestPulserSweep(self)
            %fprintf('Inside ws.WavesurferModel::completingTestPulserSweep()\n');
            self.Ephys_.completingTestPulserSweep() ;
            self.broadcast('TPUpdateTrace');
        end
    end
    
    methods (Access=protected)
%         function changeTestPulserReadiness_(self, delta)
%             self.Ephys_.changeTestPulserReadiness_(delta) ;
%         end

        function abortTestPulsing_(self)
            % This is called when a problem arises during test pulsing, and we
            % want to try very hard to get back to a known, sane, state.
            self.changeReadiness_(-1);
            self.Ephys_.abortTestPulsing() ;
            self.changeReadiness_(+1);
            if isequal(self.State,'test_pulsing') ,
                self.setState_('idle') ;
            end
        end  % function
    end

    methods
        function result = isTestPulsing(self)
            result = self.Ephys_.isTestPulsing() ;
        end
        
        function toggleIsTestPulsing(self)
            if self.isTestPulsing() , 
                self.stopTestPulsing() ;
            else
                self.startTestPulsing() ;
            end
        end

        function result = isIdle(self)
            result = isequal(self.State_, 'idle') ;
        end
        
        function result = isLoggingEnabled(self)
            result = self.Logging_.IsEnabled ;
        end
        
        function result = get.DoSubtractBaselineInTestPulseView(self)
            result = self.Ephys_.getDoSubtractBaselineInTestPulseView_() ;
        end
        
        function set.DoSubtractBaselineInTestPulseView(self, newValue)            
            self.Ephys_.setDoSubtractBaselineInTestPulseView(newValue) ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;            
        end
       
        function toggleIsGridOn(self)
            self.IsGridOn = ~(self.IsGridOn) ;
        end

        function set.IsGridOn(self, newValue)
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                self.Display_.IsGridOn = logical(newValue) ;
            else
                self.broadcast('UpdateMain');
                error('ws:invalidPropertyValue', ...
                      'IsGridOn must be a scalar, and must be logical, 0, or 1');
            end
            self.broadcast('UpdateMain');
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;            
        end
        
        function result = get.IsGridOn(self)
            result = self.Display_.IsGridOn ;
        end
        
        function toggleAreColorsNormal(self)
            self.AreColorsNormal = ~(self.AreColorsNormal) ;
        end

        function set.AreColorsNormal(self,newValue)
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                self.Display_.AreColorsNormal = logical(newValue) ;
            else
                self.broadcast('UpdateMain');
                error('ws:invalidPropertyValue', ...
                      'AreColorsNormal must be a scalar, and must be logical, 0, or 1');
            end
            self.broadcast('UpdateMain');
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;            
        end
        
        function result = get.AreColorsNormal(self)
            result = self.Display_.AreColorsNormal ;
        end
        
        function toggleDoShowZoomButtons(self)
            self.DoShowZoomButtons = ~(self.DoShowZoomButtons) ;
        end
        
        function set.DoShowZoomButtons(self, newValue)
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                self.Display_.DoShowZoomButtons = logical(newValue) ;
            else
                self.broadcast('UpdateMain');
                error('ws:invalidPropertyValue', ...
                      'DoShowZoomButtons must be a scalar, and must be logical, 0, or 1');
            end
            self.broadcast('UpdateMain');
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function toggleDoColorTraces(self)
            self.DoColorTraces = ~(self.DoColorTraces) ;
        end        
        
%         function toggleDoShowZoomButtons(self)
%             self.Display_.toggleDoShowZoomButtons_() ;
%             self.DoesProtocolNeedSave_ = true ;
%             self.broadcast('DidMaybeChangeProtocol') ;            
%         end
        
        function result = get.DoShowZoomButtons(self)
            result = self.Display_.DoShowZoomButtons ;
        end
        
%         function toggleDoColorTraces(self)
%             self.Display_.toggleDoColorTraces_() ;       
%             self.DoesProtocolNeedSave_ = true ;
%             self.broadcast('DidMaybeChangeProtocol') ;            
%         end        
        
        function result = get.DoColorTraces(self)
            result = self.Display_.DoColorTraces ;
        end
        
        function set.DoColorTraces(self,newValue)
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
                self.Display_.DoColorTraces = logical(newValue) ;
            else
                self.broadcast('UpdateMain');
                error('ws:invalidPropertyValue', ...
                      'DoColorTraces must be a scalar, and must be logical, 0, or 1');
            end
            self.broadcast('UpdateMain');
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function setPlotHeightsAndOrder(self, isDisplayed, plotHeights, rowIndexFromChannelIndex)
            self.Display_.setPlotHeightsAndOrder(isDisplayed, plotHeights, rowIndexFromChannelIndex) ;
            self.broadcast('UpdateMain') ;            
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;            
        end

        function result = getTestPulseElectrodeCommandUnits(self)
            channelName = self.getTestPulseElectrodeProperty('CommandChannelName') ;
            result = self.aoChannelUnitsFromName(channelName) ;            
        end
        
        function result = getTestPulseElectrodeMonitorUnits(self)
            channelName = self.getTestPulseElectrodeProperty('MonitorChannelName') ;
            result = self.aiChannelUnitsFromName(channelName) ;           
        end  % function

        function set.TestPulseYLimits(self, newValue)
            self.Ephys_.setTestPulseYLimits(newValue) ;
            self.broadcast('UpdateTestPulser') ;                        
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;            
        end
        
        function result = get.TestPulseYLimits(self)
            result = self.Ephys_.getTestPulseYLimits_() ;
        end
        
        function result = getCommandUnitsPerTestPulseElectrode(self)
            testPulseElectrodeIndex = self.Ephys_.TestPulseElectrodeIndex ;
            if isempty(testPulseElectrodeIndex) ,
                result = cell(1,0) ;
            else
                channelName = self.Ephys_.getElectrodeProperty(testPulseElectrodeIndex, 'CommandChannelName') ;
                channelUnits = self.aoChannelUnitsFromName(channelName) ;
                result = {channelUnits} ;
            end
%                 commandChannelNames = {commandChannelName} ;
% %             testPulseElectrodes = self.Ephys_.getTestPulseElectrodes_() ;
% %             commandChannelNames = cellfun(@(electrode)(electrode.CommandChannelName), ...
% %                                           testPulseElectrodes, ...
% %                                           'UniformOutput',false);
%                 result = cellfun(@(channelName)(self.aoChannelUnitsFromName(channelName)), ...
%                                  commandChannelNames, ...
%                                  'UniformOutput',false);
%             end
        end  % function
        
        function result = getMonitorUnitsPerTestPulseElectrode(self)
            testPulseElectrodeIndex = self.Ephys_.TestPulseElectrodeIndex ;
            if isempty(testPulseElectrodeIndex) ,
                result = cell(1,0) ;
            else
                channelName = self.Ephys_.getElectrodeProperty(testPulseElectrodeIndex, 'MonitorChannelName') ;
                channelUnits = self.aiChannelUnitsFromName(channelName) ;
                result = {channelUnits} ;
            end
%             testPulseElectrodes = self.Ephys_.getTestPulseElectrodes_() ;
%             monitorChannelNames = cellfun(@(electrode)(electrode.MonitorChannelName), ...
%                                           testPulseElectrodes, ...
%                                           'UniformOutput',false);
%             result = cellfun(@(channelName)(self.aiChannelUnitsFromName(channelName)), ...
%                              monitorChannelNames, ...
%                              'UniformOutput',false);
        end  % function
        
        function result = getIsVCPerTestPulseElectrode(self) 
            % Returns a logical row array indicated whether each trode is
            % in VC mode.  Note that to be in VC mode, from the Test
            % Pulser's point of view, is a different matter from being in
            % VC mode from the Electrode Manager's point of view.  The EM
            % mode just determines which channels get used as command and
            % monitor for the electrode.  The TP only considers an
            % electrode to be in VC if the command units are commensurable
            % (summable) with Volts, and the monitor units are
            % commensurable with Amps.
            commandUnitsPerElectrode = self.getCommandUnitsPerTestPulseElectrode() ;
            monitorUnitsPerElectrode = self.getMonitorUnitsPerTestPulseElectrode() ;
            result = ws.isVCFromMonitorAndCommandUnits(monitorUnitsPerElectrode, commandUnitsPerElectrode) ;
        end  % function

        function result = getIsCCPerTestPulseElectrode(self) 
            % Returns a logical row array indicated whether each trode is
            % in CC mode.  Note that to be in CC mode, from the Test
            % Pulser's point of view, is a different matter from being in
            % VC mode from the Electrode Manager's point of view.  The EM
            % mode just determines which channels get used as command and
            % monitor for the electrode.  The TP only considers an
            % electrode to be in CC if the command units are commensurable
            % (summable) with amps, and the monitor units are
            % commensurable with volts.
            commandUnitsPerElectrode = self.getCommandUnitsPerTestPulseElectrode() ;
            monitorUnitsPerElectrode = self.getMonitorUnitsPerTestPulseElectrode() ;
            result = ws.isCCFromMonitorAndCommandUnits(monitorUnitsPerElectrode, commandUnitsPerElectrode) ;
        end  % function
        
        function result = getGainOrResistanceUnitsPerTestPulseElectrode(self)
            if self.isTestPulsing() ,
                result = self.Ephys_.getGainOrResistanceUnitsPerTestPulseElectrodeCached_() ;
            else
                commandUnitsPerElectrode = self.getCommandUnitsPerTestPulseElectrode() ;
                monitorUnitsPerElectrode = self.getMonitorUnitsPerTestPulseElectrode() ;                
                resultIfCC = ws.divideUnits(monitorUnitsPerElectrode, commandUnitsPerElectrode) ;
                resultIfVC = ws.divideUnits(commandUnitsPerElectrode, monitorUnitsPerElectrode) ;
                isVCPerElectrode = self.getIsVCPerTestPulseElectrode() ;
                result = ws.fif(isVCPerElectrode, resultIfVC, resultIfCC) ;
            end
        end
        
        function value = getGainOrResistancePerTestPulseElectrode(self)
            value=self.Ephys_.getGainOrResistancePerTestPulseElectrode() ;
        end
        
        function [gainOrResistance, gainOrResistanceUnits] = getGainOrResistancePerTestPulseElectrodeWithNiceUnits(self)
            rawGainOrResistance = self.getGainOrResistancePerTestPulseElectrode() ;
            rawGainOrResistanceUnits = self.getGainOrResistanceUnitsPerTestPulseElectrode() ;
            % [gainOrResistanceUnits,gainOrResistance] = rawGainOrResistanceUnits.convertToEngineering(rawGainOrResistance) ;  
            [gainOrResistanceUnits,gainOrResistance] = ...
                ws.convertDimensionalQuantityToEngineering(rawGainOrResistanceUnits,rawGainOrResistance) ;
        end

        function result = areTestPulseElectrodeMonitorAndCommandChannelsOnDiffrentDevices(self)
            electrodeIndex = self.TestPulseElectrodeIndex ;
            result = self.areElectrodeMonitorAndCommandChannelsOnDifferentDevices(electrodeIndex) ;
        end
        
        function result = areElectrodeMonitorAndCommandChannelsOnDifferentDevices(self, electrodeIndex) 
            monitorDeviceName = self.getElectrodeMonitorChannelDeviceName(electrodeIndex) ;
            if isempty(monitorDeviceName) ,
                result = false ;
            else
                commandDeviceName = self.getElectrodeCommandChannelDeviceName(electrodeIndex) ;
                if isempty(commandDeviceName) ,
                    result = false ;
                else
                    result = ~isequal(monitorDeviceName, commandDeviceName) ;
                end
            end
        end
        
        function result = getElectrodeMonitorChannelDeviceName(self, electrodeIndex)
            monitorChannelName = self.Ephys_.getElectrodeProperty(electrodeIndex, 'MonitorChannelName') ;
            if isempty(monitorChannelName) ,
                result = '' ;
            else
                aiChannelIndex = self.Acquisition_.aiChannelIndexFromName(monitorChannelName) ;
                if isfinite(aiChannelIndex) ,
                    result = self.AIChannelDeviceNames{aiChannelIndex} ;
                else
                    result = '' ;
                end
            end
        end

        function result = getElectrodeCommandChannelDeviceName(self, electrodeIndex)
            commandChannelName = self.Ephys_.getElectrodeProperty(electrodeIndex, 'CommandChannelName') ;
            if isempty(commandChannelName) ,
                result = '' ;
            else
                aoChannelIndex = self.Stimulation_.aoChannelIndexFromName(commandChannelName) ;
                if isfinite(aoChannelIndex) ,
                    result = self.AOChannelDeviceNames{aoChannelIndex} ;
                else
                    result = '' ;
                end
            end
        end        
    end
    
    methods (Access=protected)
        function result = getTPElectrodeCommandChannelNames_(self)
            % this returns a cell array of length zero or one
            electrodeIndex = self.TestPulseElectrodeIndex ;
            if isempty(electrodeIndex) ,
                result = cell(1,0) ;
            else
                result = { self.getElectrodeProperty(electrodeIndex, 'CommandChannelName') } ; 
           end
        end
        
        function result = getTPElectrodeMonitorChannelNames_(self)
            % this returns a cell array of length zero or one
            electrodeIndex = self.TestPulseElectrodeIndex ;
            if isempty(electrodeIndex) ,
                result = cell(1,0) ;
            else
                result = { self.getElectrodeProperty(electrodeIndex, 'MonitorChannelName') } ; 
           end
        end
        
        function result = getCommandTerminalIDPerTestPulseElectrode_(self)
            commandChannelNames = self.getTPElectrodeCommandChannelNames_() ;
            stimulationSubsystem = self.Stimulation_ ;
            result = cellfun(@(channelName)(stimulationSubsystem.analogTerminalIDFromName(channelName)), ...
                             commandChannelNames) ;
        end  % function       

        function result = getMonitorTerminalIDPerTestPulseElectrode_(self)
            monitorChannelNames = self.getTPElectrodeMonitorChannelNames_() ;
            acquisition = self.Acquisition_ ;
            result = cellfun(@(channelName)(acquisition.analogTerminalIDFromName(channelName)), ...
                             monitorChannelNames) ;
        end        
        
        function result = getCommandChannelScalePerTestPulseElectrode_(self)
            commandChannelNames = self.getTPElectrodeCommandChannelNames_() ;
            result = cellfun(@(channelName)(self.aoChannelScaleFromName(channelName)), ...
                             commandChannelNames) ;
        end
        
        function result = getMonitorChannelScalePerTestPulseElectrode_(self)
            monitorChannelNames = self.getTPElectrodeMonitorChannelNames_() ;
            result = cellfun(@(channelName)(self.aiChannelScaleFromName(channelName)), ...
                             monitorChannelNames) ;
        end        
        
        function result = getCommandDeviceNamePerTestPulseElectrode_(self)
            commandChannelNames = self.getTPElectrodeCommandChannelNames_() ;
            stimulationSubsystem = self.Stimulation_ ;
            result = cellfun(@(channelName)(stimulationSubsystem.getDeviceNameFromChannelName(channelName)), ...
                             commandChannelNames, ...
                             'UniformOutput', false) ;
        end  % function       
        
        function result = getMonitorDeviceNamePerTestPulseElectrode_(self)
            monitorChannelNames = self.getTPElectrodeMonitorChannelNames_() ;
            acquisitionSubsystem = self.Acquisition_ ;
            result = cellfun(@(channelName)(acquisitionSubsystem.getDeviceNameFromChannelName(channelName)), ...
                             monitorChannelNames, ...
                             'UniformOutput', false) ;
        end        
        
    end  % protected methods block
    
    methods
        function zoomInTestPulseView(self)
            self.Ephys_.zoomInTestPulseView() ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function zoomOutTestPulseView(self)
            self.Ephys_.zoomOutTestPulseView() ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function scrollUpTestPulseView(self)
            self.Ephys_.scrollUpTestPulseView() ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function scrollDownTestPulseView(self)
            self.Ephys_.scrollDownTestPulseView() ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function result = get.TestPulseDuration(self) 
            result = self.Ephys_.getTestPulseDuration_() ;
        end
        
        function set.TestPulseDuration(self, newValue) 
            self.Ephys_.setTestPulseDuration(newValue) ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function result = get.IsAutoYInTestPulseView(self) 
            result = self.Ephys_.getIsAutoYInTestPulseView_() ;
        end
        
        function set.IsAutoYInTestPulseView(self, newValue) 
            self.Ephys_.setIsAutoYInTestPulseView_(newValue) ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
       
        function result = get.IsAutoYRepeatingInTestPulseView(self) 
            result = self.Ephys_.getIsAutoYRepeatingInTestPulseView_() ;
        end
        
        function set.IsAutoYRepeatingInTestPulseView(self, newValue) 
            self.Ephys_.setIsAutoYRepeatingInTestPulseView_(newValue) ;
            self.broadcast('UpdateTestPulser') ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function value = getUpdateRateInTestPulseView(self)
            value = self.Ephys_.getUpdateRateInTestPulseView_() ;
        end        

        function result = getTestPulseMonitorTraceTimeline(self)
            fs = self.AcquisitionSampleRate ;
            result = self.Ephys_.getTestPulseMonitorTraceTimeline_(fs) ;
        end  % function         
        
        function result = getTestPulseMonitorTrace(self)
            result = self.Ephys_.getTestPulseMonitorTrace() ;
        end  % function         
        
        function setElectrodeProperty(self, electrodeIndex, propertyName, newValue)
            switch propertyName ,
                case 'Type' ,
                    self.setElectrodeType_(electrodeIndex, newValue) ;
                    self.DoesProtocolNeedSave_ = true ;
                    self.broadcast('DidMaybeChangeProtocol') ;
                case 'IndexWithinType' ,
                    self.setElectrodeIndexWithinType_(electrodeIndex, newValue) ;
                    self.broadcast('UpdateElectrodeManager') ;
                    self.DoesProtocolNeedSave_ = true ;
                    self.broadcast('DidMaybeChangeProtocol') ;
                otherwise ,
                    % the common case
                    try
                        self.Ephys_.setElectrodeProperty(electrodeIndex, propertyName, newValue) ;
                        self.broadcast('UpdateElectrodeManager') ;
                        self.broadcast('UpdateTestPulser') ;
                        self.DoesProtocolNeedSave_ = true ;
                    catch exception
                        % deal with EPCMasterSocket exceptions,
                        % otherwise rethrow
                        indicesThatMatch=strfind(exception.identifier,'EPCMasterSocket:');                
                        if ~isempty(indicesThatMatch) && indicesThatMatch(1)==1 ,
                            % The error was an EPCMasterSocket error,
                            % so we just neglect to set the mode in the
                            % electrode, and make sure the view gets
                            % resynced (happens below now)
                            %self.electrodeMayHaveChanged(electrodeIndex, propertyName);
                        else
                            % There was some other kind of problem
                            rethrow(exception);
                        end
                    end                                        
                    % Notify other systems that an electrode may have changed in the Ephys
                    % subsystem.
                    isModeOrChannelNameOrScale = ...
                        isempty(propertyName) || ...
                        ismember(propertyName, ...
                                 {'Mode' ...
                                  'CommandChannelName' 'CommandScaling' ...
                                  'MonitorChannelName' 'MonitorScaling' ...
                                  'VoltageCommandChannelName' 'VoltageCommandScaling' ...
                                  'CurrentCommandChannelName' 'CurrentCommandScaling' ...
                                  'VoltageMonitorChannelName' 'VoltageMonitorScaling' ...
                                  'CurrentMonitorChannelName' 'CurrentMonitorScaling' }) ;
                    if isModeOrChannelNameOrScale ,
                        %self.Display_.didSetAnalogChannelUnitsOrScales() ;
                        %self.broadcast('UpdateTraces') ;
                        self.broadcast('UpdateMain') ;
                        self.broadcast('UpdateChannels') ;
                    end                                
                    self.broadcast('DidMaybeChangeProtocol') ;
            end
        end
        
%         function setElectrodeModeAndScalings(self,...
%                                              electrodeIndex, ...
%                                              newMode, ...
%                                              newCurrentMonitorScaling, ...
%                                              newVoltageMonitorScaling, ...
%                                              newCurrentCommandScaling, ...
%                                              newVoltageCommandScaling,...
%                                              newIsCommandEnabled)
%             self.Ephys_.setElectrodeModeAndScalings_(electrodeIndex, ...
%                                                      newMode, ...
%                                                      newCurrentMonitorScaling, ...
%                                                      newVoltageMonitorScaling, ...
%                                                      newCurrentCommandScaling, ...
%                                                      newVoltageCommandScaling,...
%                                                      newIsCommandEnabled) ;
%             self.DoesProtocolNeedSave_ = true ;
%             self.Display_.didSetAnalogChannelUnitsOrScales() ;
%             self.broadcast('UpdateChannels') ;
%         end  % function
        
        function result = areTestPulseElectrodeChannelsValid(self)
            aiChannelNames = self.AIChannelNames ;
            isAIChannelActive = self.Acquisition_.IsAnalogChannelActive ;            
            activeAIChannelNames = aiChannelNames(isAIChannelActive) ;
            aoChannelNames = self.AOChannelNames ;
            result = self.Ephys_.areTestPulseElectrodeChannelsValid(activeAIChannelNames, aoChannelNames) ;
        end  % function

        function updateSmartElectrodeGainsAndModes(self)
            self.changeReadiness_(-1) ;
            % Get the current mode and scaling from any smart electrodes
            smartElectrodeTypes = setdiff(ws.Electrode.Types,{'Manual'}) ;
            didSetSomething = false ;
            for k = 1:length(smartElectrodeTypes) , 
                smartElectrodeType = smartElectrodeTypes{k} ;                
                [areAnyOfThisType, ...
                 indicesOfThisTypeOfElectrodes, ...
                 overallError, ...
                 modes, ...
                 currentMonitorScalings, voltageMonitorScalings, currentCommandScalings, voltageCommandScalings, ...
                 isCommandEnabled] = ...
                    self.Ephys_.probeHardwareForSmartElectrodeModesAndScalings_(smartElectrodeType) ;
                if areAnyOfThisType && isempty(overallError) ,
                    nElectrodesOfThisType = length(indicesOfThisTypeOfElectrodes) ;
                    for j = 1:nElectrodesOfThisType ,
                        electrodeIndex = indicesOfThisTypeOfElectrodes(j) ;
                        % Even if there's was an error on the electrode
                        % and no new info could be gathered, those ones
                        % should just be nan's or empty's, which
                        % setModeAndScalings() knows to ignore.
%                         self.setElectrodeModeAndScalings(electrodeIndex, ...
%                                                          modes{j}, ...
%                                                          currentMonitorScalings(j), ...
%                                                          voltageMonitorScalings(j), ...
%                                                          currentCommandScalings(j), ...
%                                                          voltageCommandScalings(j), ...
%                                                          isCommandEnabled{j}) ;
                        self.Ephys_.setElectrodeModeAndScalings(electrodeIndex, ...
                                                                modes{j}, ...
                                                                currentMonitorScalings(j), ...
                                                                voltageMonitorScalings(j), ...
                                                                currentCommandScalings(j), ...
                                                                voltageCommandScalings(j), ...
                                                                isCommandEnabled{j}) ;                                                            
                        didSetSomething = true ;
                    end
                end
            end
            if didSetSomething ,
                %self.DoesProtocolNeedSave_ = true ;
                % Should maybe be smarter about this, but it's annoying to
                % have the protocol think it needs saving after each
                % press of the Update button.
                self.syncTraces_() ;                
                self.broadcast('UpdateElectrodeManager') ;
                self.broadcast('UpdateTestPulser') ;
                self.broadcast('UpdateTraces') ;
                self.broadcast('UpdateMain') ;
                self.broadcast('UpdateChannels') ;
            end
            self.changeReadiness_(+1) ;
            self.broadcast('UpdateElectrodes') ;
        end  % function
        
        function result = isElectrodeOfType(self, queryType)
            result = self.Ephys_.isElectrodeOfType(queryType) ;
        end  % function
        
        function reconnectWithSmartElectrodes(self)
            % Close and repoen the connection to any smart electrodes
            self.changeReadiness_(-1) ;
            self.Ephys_.reconnectWithSmartElectrodes_() ;
            self.updateSmartElectrodeGainsAndModes() ;
            self.changeReadiness_(+1) ;
            self.broadcast('UpdateElectrodes');
        end  % function
    end
    
    methods (Access=protected)
        function setElectrodeType_(self, electrodeIndex, newValue)
            % can only change the electrode type if softpanels are
            % enabled.  I.e. only when WS is _not_ in command of the
            % gain settings
            self.changeReadiness_(-1);  % may have to establish contact with the softpanel, which can take a little while
            doNeedToUpdateGainsAndModes = self.Ephys_.setElectrodeType(electrodeIndex, newValue) ;
            self.broadcast('UpdateElectrodeManager') ;            
            if doNeedToUpdateGainsAndModes, 
                self.updateSmartElectrodeGainsAndModes() ;
            end
            self.changeReadiness_(+1);
        end  % function
       
        function setElectrodeIndexWithinType_(self, electrodeIndex, newValue)
            doUpdateSmartElectrodeGainsAndModes = self.Ephys_.setElectrodeIndexWithinType(electrodeIndex, newValue) ;
            if doUpdateSmartElectrodeGainsAndModes ,
                self.updateSmartElectrodeGainsAndModes() ;
            end
        end
    end  % protected methods block
    
    methods
        function toggleIsInControlOfSoftpanelModeAndGains(self)
            currentValue = self.IsInControlOfSoftpanelModeAndGains ;
            self.IsInControlOfSoftpanelModeAndGains = ~currentValue ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end        
        
        function set.IsInControlOfSoftpanelModeAndGains(self, newValue)
            if self.areAnyElectrodesCommandable() ,
                try
                    doUpdateSmartElectrodeGainsAndModes = self.Ephys_.setIsInControlOfSoftpanelModeAndGains(newValue) ;
                catch err
                    self.broadcast('UpdateElectrodeManager') ;
                    rethrow(err) ;
                end
                self.broadcast('UpdateElectrodeManager') ;
                self.DoesProtocolNeedSave_ = true ;
                self.broadcast('DidMaybeChangeProtocol') ;
                if doUpdateSmartElectrodeGainsAndModes ,
                    self.updateSmartElectrodeGainsAndModes() ;
                end
            end
        end
        
        function electrodeIndex = addNewElectrode(self)
            try 
                electrodeIndex = self.Ephys_.addNewElectrode() ;
            catch err
                self.broadcast('UpdateElectrodeManager') ;
                self.broadcast('UpdateTestPulser') ;
                rethrow(err) ;
            end
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateElectrodeManager') ;
            self.broadcast('UpdateTestPulser') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
       
        function removeMarkedElectrodes(self)
            try
                self.Ephys_.removeMarkedElectrodes() ;
                self.syncTraces_() ;
                self.DoesProtocolNeedSave_ = true ;
            catch err
                self.broadcast('UpdateElectrodeManager');
                self.broadcast('UpdateTestPulser') ;
                self.broadcast('UpdateTraces') ;
                self.broadcast('UpdateMain') ;
                self.broadcast('UpdateChannels') ;
                rethrow(err) ;
            end
            self.broadcast('UpdateElectrodeManager');            
            self.broadcast('UpdateTestPulser') ;
            self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateChannels') ;
        end
        
        function set.DoTrodeUpdateBeforeRun(self, newValue)
            self.Ephys_.setDoTrodeUpdateBeforeRun(newValue) ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateElectrodeManager') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end        
       
        function result = get.DoTrodeUpdateBeforeRun(self)
            result = self.Ephys_.getDoTrodeUpdateBeforeRun() ;
        end

%         function setElectrodeModeOrScaling(self, electrodeIndex, propertyName, newValue)
%             try
%                 self.Ephys_.setElectrodeModeOrScaling_(electrodeIndex, propertyName, newValue) ;
%             catch exception
%                 % deal with EPCMasterSocket exceptions,
%                 % otherwise rethrow
%                 indicesThatMatch=strfind(exception.identifier,'EPCMasterSocket:');                
%                 if ~isempty(indicesThatMatch) && indicesThatMatch(1)==1 ,
%                     % The error was an EPCMasterSocket error,
%                     % so we just neglect to set the mode in the
%                     % electrode, and make sure the view gets
%                     % resynced
%                     self.electrodeMayHaveChanged(electrodeIndex, propertyName);
%                 else
%                     % There was some other kind of problem
%                     rethrow(exception);
%                 end
%             end                                        
%         end
        
%         function result = get.IsElectrodeMarkedForTestPulse(self)
%             result = self.Ephys_.getIsElectrodeMarkedForTestPulse() ;
%         end
        
%         function set.IsElectrodeMarkedForTestPulse(self, newValue)
%             % Don't want to allow this anymore, so just ignore new value
%             %self.Ephys_.setIsElectrodeMarkedForTestPulse_(newValue) ;
%             self.broadcast('UpdateElectrodes') ;
%         end        
        
        function result = get.IsElectrodeMarkedForRemoval(self)
            result = self.Ephys_.getIsElectrodeMarkedForRemoval_() ;
        end
        
        function set.IsElectrodeMarkedForRemoval(self, newValue)
            self.Ephys_.setIsElectrodeMarkedForRemoval_(newValue) ;
            %self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateElectrodes') ;
            %self.broadcast('DidMaybeChangeProtocol') ;
        end        
        
        function result = get.TestPulseElectrodeIndex(self)
            result = self.Ephys_.TestPulseElectrodeIndex ;
        end

        function set.TestPulseElectrodeIndex(self, newValue)
            self.Ephys_.TestPulseElectrodeIndex = newValue ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function setTestPulseElectrodeProperty(self, propertyName, newValue)
            testPulseElectrodeIndex = self.TestPulseElectrodeIndex ;
            if ~isempty(testPulseElectrodeIndex) ,
                self.setElectrodeProperty(testPulseElectrodeIndex, propertyName, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
                self.broadcast('DidMaybeChangeProtocol') ;
            end
        end
        
        function value=get.NextRunAbsoluteFileName(self)
            logging = self.Logging_ ;            
            firstSweepIndex = logging.NextSweepIndex ;
            numberOfSweeps = self.NSweepsPerRun ;
            fileName = logging.sweepSetFileNameFromNumbers(firstSweepIndex, numberOfSweeps) ;
            value = fullfile(logging.FileLocation, fileName);
        end  % function
        
        function set.UserClassName(self, newValue)
            if ws.isString(newValue) ,
                % If it's a string, we'll keep it, but we have to check if
                % it's a valid class name
                trimmedValue = strtrim(newValue) ;
                err = self.UserCodeManager_.setClassName_(trimmedValue) ;
                self.DoesProtocolNeedSave_ = true ;
                self.broadcast('DidMaybeSetUserClassName');
                self.broadcast('DidMaybeChangeProtocol') ;
                if isempty(err) ,
                    self.callUserMethod_('wake');  % wake the user object
                else
                  error('wavesurfer:errorWhileInstantiatingUserObject', ...
                        'Unable to instantiate user object: %s.',err.message);
                end
            else
                self.broadcast('DidMaybeSetUserClassName');  % replace the bad value with the old value in the view
                error('ws:invalidPropertyValue', ...
                      'Invalid value for property ''ClassName'' supplied.');
            end
        end
        
        function result = get.UserClassName(self) 
            result = self.UserCodeManager_.getClassName_() ;
        end
        
        function reinstantiateUserObject(self)
            err = self.UserCodeManager_.reinstantiateUserObject() ;
            if ~isempty(err) ,
                self.broadcast('UpdateUserCodeManager');
                throw(err) ;
            end
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateUserCodeManager');
            self.broadcast('DidMaybeChangeProtocol') ;
            self.callUserMethod_('wake');  % wake the user object
        end        
        
        function result = get.IsUserClassNameValid(self)
            result = self.UserCodeManager_.getIsClassNameValid_() ;
        end
        
        function result = get.DoesTheUserObjectMatchTheUserClassName(self)
            result = self.UserCodeManager_.getDoesTheObjectMatchClassName_() ;
        end               

        function result = get.TheUserObject(self)
            result = self.UserCodeManager_.getTheObject_() ;
        end               

        function result = getElectrodeProperty(self, electrodeIndex, propertyName)
            result = self.Ephys_.getElectrodeProperty(electrodeIndex, propertyName) ;
        end  % function

        function result = get.ElectrodeCount(self)
            result = self.Ephys_.getElectrodeCount() ;
        end               
        
%         function electrode = getElectrodeByIndex(self, electrodeIndex)
%             electrode = self.Ephys_.getElectrodeByIndex_(electrodeIndex) ;
%         end    
        
        function result = get.IsInControlOfSoftpanelModeAndGains(self)
            result = self.Ephys_.getIsInControlOfSoftpanelModeAndGains_() ;
        end

%         function set.IsInControlOfSoftpanelModeAndGains(self, newValue)
%             self.Ephys_.setIsInControlOfSoftpanelModeAndGains_(newValue) ;
%         end

        function result = areAnyElectrodesCommandable(self)
            result = self.Ephys_.areAnyElectrodesCommandable() ;
        end  % function
        
        function result = get.DidLastElectrodeUpdateWork(self)
            result = self.Ephys_.DidLastElectrodeUpdateWork ;
        end
        
        function result = getNumberOfElectrodesClaimingMonitorChannel(self, queryChannelNames)
            result = self.Ephys_.getNumberOfElectrodesClaimingMonitorChannel(queryChannelNames) ;
        end
        
        function result = getNumberOfElectrodesClaimingCommandChannel(self, queryChannelNames)
            result = self.Ephys_.getNumberOfElectrodesClaimingCommandChannel(queryChannelNames) ;
        end
        
        function result = get.AreSoftpanelsEnabled(self)
            result = self.Ephys_.AreSoftpanelsEnabled ;
        end

%         function set.AreSoftpanelsEnabled(self, newValue)
%             self.Ephys_.AreSoftpanelsEnabled = newValue ;
%         end
        
        function result = doesElectrodeHaveCommandOnOffSwitch(self)
            result = self.Ephys_.doesElectrodeHaveCommandOnOffSwitch() ;
        end        
        
        function result = get.IsDoTrodeUpdateBeforeRunSensible(self)
            result = self.Ephys_.IsDoTrodeUpdateBeforeRunSensible() ;
        end        
        
        function result = areAnyElectrodesSmart(self)
            result = self.Ephys_.areAnyElectrodesSmart() ;
        end        
        
%         function result = areAllMonitorAndCommandChannelNamesDistinct(self)
%             result = self.Ephys_.areAllMonitorAndCommandChannelNamesDistinct() ;
%         end  % function
        
        function result = getTestPulseElectrodeNames(self)
            result = self.Ephys_.getTestPulseElectrodeNames() ;
        end
        
%         function subscribeMeToEphysEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Ephys_.subscribeMe(subscriber,eventName,propertyName,methodName) ;
%         end
        
%         function subscribeMeToElectrodeManagerEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Ephys_.subscribeMeToElectrodeManagerEvent(subscriber,eventName,propertyName,methodName) ;
%         end
        
%         function subscribeMeToTestPulserEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Ephys_.subscribeMeToTestPulserEvent(subscriber,eventName,propertyName,methodName) ;
%         end
        
        function result = get.TestPulseElectrodesCount(self)
            result = self.Ephys_.TestPulseElectrodesCount ;
        end
        
        function result = getAllElectrodeNames(self)
            result = self.Ephys_.getAllElectrodeNames() ;
        end
        
%         function result = get.TestPulseElectrodeAmplitude(self)
%             result = self.Ephys_.TestPulseElectrodeAmplitude ;
%         end
%         
%         function set.TestPulseElectrodeAmplitude(self, newValue)  % in units of the electrode command channel
%             self.Ephys_.TestPulseElectrodeAmplitude = newValue ;
%         end  % function         
        
%         function set.TestPulseElectrodeName(self, newValue)
%             self.Ephys_.TestPulseElectrodeName = newValue ;
%         end        
        
%         function result = getIsTestPulserReady(self)
%             result = self.Ephys_.getIsTestPulserReady() ;
%         end
        
        function setTestPulseElectrodeByName(self, newValue)
            self.Ephys_.setTestPulseElectrodeByName(newValue) ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('UpdateTestPulser');
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function result = getTestPulseElectrodeProperty(self, propertyName)
            result = self.Ephys_.getTestPulseElectrodeProperty(propertyName) ;
        end        
        
        function result = get.IsStimulationEnabled(self)
            result = self.Stimulation_.IsEnabled ;
        end

        function set.IsStimulationEnabled(self, newValue)
            try
                self.Stimulation_.IsEnabled = newValue ;
                self.DoesProtocolNeedSave_ = true ;
            catch me
                self.broadcast('UpdateGeneral') ;
                rethrow(me) ;
            end                
            self.broadcast('UpdateGeneral') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end

        function result = get.IsLoggingEnabled(self)
            result = self.Logging_.IsEnabled ;
        end

        function result = get.IsDisplayEnabled(self)
            result = self.Display_.IsEnabled ;
        end

        function set.IsDisplayEnabled(self, newValue)
            self.Display_.IsEnabled = newValue ;
            self.DoesProtocolNeedSave_ = true ;
            %self.broadcast('UpdateTraces') ;
            self.broadcast('UpdateMain') ;
            self.broadcast('UpdateGeneral') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end

        function result = getNActiveAIChannels(self)
            result = self.Acquisition_.NActiveAnalogChannels ;
        end

        function result = getNActiveDIChannels(self)
            result = self.Acquisition_.NActiveDigitalChannels ;
        end
        
        function result = get.DataFileLocation(self)
            result = self.Logging_.FileLocation ;
        end

        function set.DataFileLocation(self, newValue)
            if ws.isString(newValue) ,
                self.Logging_.FileLocation = newValue ;
            else
                self.broadcast('UpdateLogging');
                error('ws:invalidPropertyValue', ...
                      'DataFileLocation must be a string');                    
            end
            self.broadcast('UpdateLogging');            
        end
        
        function result = get.DataFileBaseName(self)
            result = self.Logging_.FileBaseName ;
        end
        
        function set.DataFileBaseName(self, newValue)
            if ws.isString(newValue) ,
                self.Logging_.FileBaseName = newValue ;
            else
                self.broadcast('UpdateLogging') ;
                error('ws:invalidPropertyValue', ...
                      'DataFileBaseName must be a string');                    
            end
            self.broadcast('UpdateLogging');                        
        end
                
        function result = get.IsOKToOverwriteDataFile(self)
            result = self.Logging_.IsOKToOverwrite ;
        end
        
        function set.IsOKToOverwriteDataFile(self, newValue)            
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && isfinite(newValue))) ,
                self.Logging_.IsOKToOverwrite = newValue ;
            else
                self.broadcast('UpdateLogging');
                error('ws:invalidPropertyValue', ...
                      'IsOKToOverwriteDataFile must be a logical scalar, or convertable to one');                  
            end
            self.broadcast('UpdateLogging');                                    
        end
        
        function result = get.NPlots(self) 
            result = self.Display_.NPlots ;
        end
        
        function set.NextSweepIndex(self, newValue)
            if isnumeric(newValue) && isreal(newValue) && isscalar(newValue) && (newValue==round(newValue)) && newValue>=0 ,
                self.Logging_.NextSweepIndex = newValue ;
            else
                self.broadcast('UpdateLogging');
                error('ws:invalidPropertyValue', ...
                      'NextSweepIndex must be a (scalar) nonnegative integer');
            end
            self.broadcast('UpdateLogging');                        
        end
        
        function result = get.NextSweepIndex(self)
            result = self.Logging_.NextSweepIndex ;
        end
        
        function result = getStimulusLibraryCopy(self)
            result = self.Stimulation_.getStimulusLibraryCopy() ;
        end
        
        function value = get.DisplayUpdateRate(self)
            value = self.Display_.UpdateRate ;
        end
        
        function set.DisplayUpdateRate(self, newValue)
            if isnumeric(newValue) && isscalar(newValue) && isfinite(newValue) && newValue>0 ,
                self.Display_.UpdateRate = max(0.1,min(newValue,10)) ;
            else
                self.broadcast('DidSetUpdateRate');
                error('ws:invalidPropertyValue', ...
                      'UpdateRate must be a scalar finite positive number') ;
            end
            self.broadcast('DidSetUpdateRate');
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
%         function mimicStimulusLibrary_(self, newValue) 
%             self.Stimulation_.mimicStimulusLibrary_(newValue) ;
%         end
        
        function result = get.AIScalingCoefficients(self)
            result = self.Acquisition_.AnalogScalingCoefficients ;
        end
        
        function value = get.XOffset(self)
            value = self.Display_.XOffset ;
        end
           
        function set.XOffset(self, newValue)
            if isnumeric(newValue) && isscalar(newValue) && isfinite(newValue) ,
                self.Display_.XOffset = double(newValue) ;
            else
                self.broadcast('DidSetXOffset');
                error('ws:invalidPropertyValue', ...
                      'XOffset must be a scalar finite number') ;
            end
            self.broadcast('DidSetXOffset');
        end
        
        function value = get.YLimitsPerAIChannel(self)
            value = self.Display_.YLimitsPerAnalogChannel ;
        end

        function result = get.AreYLimitsLockedTightToDataForAIChannel(self)
            result = self.Display_.AreYLimitsLockedTightToDataForAnalogChannel ;
        end
        
        function result = get.ChannelIndexWithinTypeFromPlotIndex(self)
            result = self.Display_.ChannelIndexWithinTypeFromPlotIndex ;
        end

        function result = get.IsAnalogFromPlotIndex(self)
            result = self.Display_.IsAnalogFromPlotIndex ;
        end
        
        function result = get.ChannelIndexFromPlotIndex(self)
            result = self.Display_.ChannelIndexFromPlotIndex ;
        end

        function result = get.IsInputChannelInCacheFromInputChannelIndex(self)
            % Cache input channels are the input channels that were active
            % at the time of the last sweep.
            result = [ self.Acquisition_.IsInCacheFromAnalogChannelIndex self.Acquisition_.IsInCacheFromDigitalChannelIndex ] ;
        end        
        
        function result = get.CacheInputChannelIndexFromInputChannelIndex(self)
            % Cache input channels are the input channels that were active
            % at the time of the last sweep.
            result = [ self.Acquisition_.IndexInCacheFromAnalogChannelIndex self.Acquisition_.IndexInCacheFromDigitalChannelIndex ] ;
        end
        
        function result = get.PlotHeightFromPlotIndex(self)
            result = self.Display_.PlotHeightFromPlotIndex ;
        end
        
%         function subscribeMeToDisplayEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Display_.subscribeMe(subscriber,eventName,propertyName,methodName) ;
%         end
        
%         function subscribeMeToStimulationEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Stimulation_.subscribeMe(subscriber,eventName,propertyName,methodName) ;
%         end
        
%         function subscribeMeToLoggingEvent(self,subscriber,eventName,propertyName,methodName)
%             self.Logging_.subscribeMe(subscriber,eventName,propertyName,methodName) ;
%         end
        
        function result = get.PlotIndexFromChannelIndex(self)
            result = self.Display_.PlotIndexFromChannelIndex ;
        end
        
        function setAreYLimitsLockedTightToDataForSingleAIChannel_(self, aiChannelIndex, newValue)            
            % Underscore b/c doesn't trigger an update
            self.Display_.setAreYLimitsLockedTightToDataForSingleChannel_(aiChannelIndex, newValue) ;
        end        
        
%         function setYLimitsForSingleAIChannel_(self, aiChannelIndex, newValue)
%             % Underscore b/c doesn't trigger an update
%             self.Display_.setYLimitsForSingleAIChannel_(aiChannelIndex, newValue) ;
%         end

%         function subscribeMeToUserCodeManagerEvent(self,subscriber,eventName,propertyName,methodName)
%             self.UserCodeManager_.subscribeMe(subscriber,eventName,propertyName,methodName) ;
%         end
        
        function result = get.AIChannelTerminalNames(self)
            result = self.Acquisition_.AnalogTerminalNames ;
        end
        
        function result = get.DIChannelTerminalNames(self)
            result = self.Acquisition_.DigitalTerminalNames ;
        end
        
        function result = get.AOChannelTerminalNames(self)
            result = self.Stimulation_.AnalogTerminalNames ;
        end
        
        function result = get.DOChannelTerminalNames(self)
            result = self.Stimulation_.DigitalTerminalNames ;
        end
        
        function result = isStimulusLibrarySelfConsistent(self)
            result = self.Stimulation_.isStimulusLibrarySelfConsistent() ;
        end
        
        function result = get.DoRepeatStimulusSequence(self)
            result= self.Stimulation_.DoRepeatSequence ;
        end
        
        function set.DoRepeatStimulusSequence(self, newValue)
            try
                self.Stimulation_.DoRepeatSequence = newValue ;
                self.DoesProtocolNeedSave_ = true ;
            catch me
                self.broadcast('UpdateGeneral') ;
                rethrow(me) ;
            end
            self.broadcast('UpdateGeneral') ;            
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function setYLimitsForSinglePlot(self, plotIndex, newValue)
            if isnumeric(newValue) && isequal(size(newValue),[1 2]) && newValue(1)<=newValue(2) ,
                aiChannelIndex = self.Display_.setYLimitsForSinglePlot(plotIndex, double(newValue)) ;
                self.broadcast('DidSetYAxisLimits', plotIndex, aiChannelIndex);
                self.DoesProtocolNeedSave_ = true ;
                self.broadcast('DidMaybeChangeProtocol') ;
            else
                error('ws:invalidPropertyValue', ...
                      'y limits must be 2 element numeric row vector, with the first element less than or equal to the second') ;                
            end
        end

        function scrollUp(self, plotIndex)  % works on analog channels only
            channelIndex = self.Display_.scrollUp(plotIndex) ;
            self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
        
        function scrollDown(self, plotIndex)  % works on analog channels only
            channelIndex = self.Display_.scrollDown(plotIndex) ;
            self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
                
        function zoomIn(self, plotIndex)  % works on analog channels only
            channelIndex = self.Display_.zoomIn(plotIndex) ;
            self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end  % function
                
        function zoomOut(self, plotIndex)  % works on analog channels only
            channelIndex = self.Display_.zoomOut(plotIndex) ;
            self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex) ;
            self.DoesProtocolNeedSave_ = true ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end        
        
        function result = get.DoIncludeDateInDataFileName(self)
            result = self.Logging_.DoIncludeDate ;
        end
        
        function set.DoIncludeDateInDataFileName(self, newValue)            
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && ~isnan(newValue))) ,
                self.Logging_.DoIncludeDate = newValue ;
            else
                self.broadcast('UpdateLogging') ;
                error('ws:invalidPropertyValue', ...
                      'DoIncludeDateInDataFileName must be a logical scalar, or convertable to one') ;                  
            end
            self.broadcast('UpdateLogging') ;                        
        end

        function result = get.DoIncludeSessionIndexInDataFileName(self)
            result = self.Logging_.DoIncludeSessionIndex ;
        end
        
        function set.DoIncludeSessionIndexInDataFileName(self, newValue)
            self.Logging_.DoIncludeSessionIndex = newValue ;
            
            if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && ~isnan(newValue))) ,
                self.Logging_.DoIncludeSessionIndex = newValue ;
            else
                self.broadcast('UpdateDoIncludeSessionIndexInDataFileName');
                error('ws:invalidPropertyValue', ...
                      'DoIncludeSessionIndexInDataFileName must be a logical scalar, or convertable to one');                  
            end
            self.broadcast('UpdateDoIncludeSessionIndexInDataFileName');            
            
        end
        
        function result = get.SessionIndex(self)
            result = self.Logging_.SessionIndex ;
        end        
        
        function set.SessionIndex(self, newValue)
            if self.DoIncludeSessionIndexInDataFileName ,
                if isnumeric(newValue) && isscalar(newValue) && round(newValue)==newValue && newValue>=1 ,
                    self.Logging_.SessionIndex = newValue ;
                else
                    self.broadcast('UpdateLogging');
                    error('ws:invalidPropertyValue', ...
                          'SessionIndex must be an integer greater than or equal to one');
                end
            else
                self.broadcast('UpdateLogging');
                error('ws:invalidPropertyValue', ...
                      'Can''t set SessionIndex when DoIncludeSessionIndexInDataFileName is false');
            end
            self.broadcast('UpdateLogging');            
        end

%         function result = stimulusLibrary(self)
%             % Note that this returns a *copy* of the internal stimulus library
%             result = self.Stimulation_.getStimulusLibraryCopy() ;
%         end
        
        function result = get.CurrentRunAbsoluteFileName(self)
            result = self.Logging_.CurrentRunAbsoluteFileName ;
        end                
        
        function result = isIdleSensuLato(self)
            state = self.State ;
            result = isequal(state,'idle') || isequal(state,'no_device') ;
        end

        function result = get.PlotHeightFromAIChannelIndex(self)
            result = self.Display_.PlotHeightFromAnalogChannelIndex ;
        end
        
        function result = get.PlotHeightFromDIChannelIndex(self)
            result = self.Display_.PlotHeightFromDigitalChannelIndex ;
        end
        
        function result = get.RowIndexFromAIChannelIndex(self)
            result = self.Display_.RowIndexFromAnalogChannelIndex ;
        end
        
        function result = get.RowIndexFromDIChannelIndex(self)
            result = self.Display_.RowIndexFromDigitalChannelIndex ;
        end
        
        function result = get.NRunsCompleted(self)
            result = self.NRunsCompleted_ ;
        end
        
        function result = get.AIChannelDeviceNames(self)
            result = self.Acquisition_.AnalogDeviceNames ;
        end
        
        function result = get.AOChannelDeviceNames(self)
            result = self.Stimulation_.AnalogDeviceNames ;
        end
        
        function result = get.DIChannelDeviceNames(self)
            result = self.Acquisition_.DigitalDeviceNames ;
        end
        
        function result = get.DOChannelDeviceNames(self)
            result = self.Stimulation_.DigitalDeviceNames ;
        end
        
        function result = get.AIChannelTerminalIDs(self)
            result = self.Acquisition_.AnalogTerminalIDs ;
        end
        
        function result = get.AOChannelTerminalIDs(self)
            result = self.Stimulation_.AnalogTerminalIDs ;
        end
        
        function result = get.DIChannelTerminalIDs(self)
            result = self.Acquisition_.DigitalTerminalIDs ;
        end
        
        function result = get.DOChannelTerminalIDs(self)
            result = self.Stimulation_.DigitalTerminalIDs ;
        end
        
        function setSingleAIChannelDeviceName(self, i, newValue)
            if 1<=i && i<=self.NAIChannels && i==round(i) && ws.isString(newValue) && ismember(newValue, self.AllDeviceNames) ,
                self.Acquisition_.setSingleAnalogDeviceName(i, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
                self.syncIsAIChannelTerminalOvercommitted_() ;
            else
                self.broadcast('UpdateChannels') ;
                error('ws:invalidPropertyValue', ...
                      'The AI channel index must be an integer between 1 and %d, and the value must be a valid device name', self.NAIChannels);
            end
            self.broadcast('UpdateChannels') ;
            self.broadcast('UpdateElectrodes') ;  % Changing the device name can put channels for same trode on diff devices, which is a problem
            self.broadcast('UpdateTestPulser') ;            
        end  % function
        
%         function setSingleDIChannelDeviceName(self, i, newValue)
%             if 1<=i && i<=self.NDIChannels && i==round(i) && ws.isString(newValue) && ismember(newValue, self.AllDeviceNames) ,
%                 self.Acquisition_.setSingleDigitalDeviceName(i, newValue) ;
%                 self.syncIsDIOChannelTerminalOvercommitted_() ;
%             else
%                 self.broadcast('UpdateChannels') ;
%                 error('ws:invalidPropertyValue', ...
%                       'The DI channel index must be an integer between 1 and %d, and the value must be a valid device name', self.NDIChannels);
%             end                
%             self.broadcast('UpdateChannels') ;
%         end  % function
        
        function setSingleAOChannelDeviceName(self, i, newValue)
            if 1<=i && i<=self.NAOChannels && i==round(i) && ws.isString(newValue) && ismember(newValue, self.AllDeviceNames) ,
                self.Stimulation_.setSingleAnalogDeviceName(i, newValue) ;
                self.DoesProtocolNeedSave_ = true ;
                self.syncIsAOChannelTerminalOvercommitted_() ;
            else                
                self.broadcast('UpdateChannels') ;
                error('ws:invalidPropertyValue', ...
                      'The AO channel index must be an integer between 1 and %d, and the value must be a valid device name', self.NAOChannels);
            end                
            self.broadcast('UpdateChannels') ;
            self.broadcast('UpdateElectrodes') ;  % Changing the device name can put channels for same trode on diff devices, which is a problem
            self.broadcast('UpdateTestPulser') ;
        end  % function                
        
%         function setSingleDOChannelDeviceName(self, i, newValue)
%             if 1<=i && i<=self.NDOChannels && i==round(i) && ws.isString(newValue) && ismember(newValue, self.AllDeviceNames) ,
%                 self.Stimulation_.setSingleDigitalDeviceName(i, newValue) ;
%                 self.syncIsDIOChannelTerminalOvercommitted_() ;
%             else
%                 self.broadcast('UpdateChannels') ;
%                 error('ws:invalidPropertyValue', ...
%                       'The DO channel index must be an integer between 1 and %d, and the value must be a valid device name', self.NDOChannels);
%             end                
%             self.broadcast('UpdateChannels') ;
%         end  % function

        function result = getPrimaryDeviceIndexMaybe(self)
            % The index of self.PrimaryDeviceName in self.AllDeviceNames.  Returns empty if
            % self.PrimaryDeviceName is not in self.AllDeviceNames.
            result = self.getDeviceIndexFromName(self.PrimaryDeviceName) ;
        end
        
        function result = getDeviceIndexFromName(self, deviceName)
            % The index of deviceName in self.AllDeviceNames.  Returns empty if
            % deviceName is not in self.AllDeviceNames.
            allDeviceNames = self.AllDeviceNames ;
            isAMatch = strcmpi(deviceName, allDeviceNames) ;  % DAQMX device names are not case-sensitive
            result = find(isAMatch,1) ;
        end
        
        function setSingleAOChannelTerminalName(self, iChannel, terminalName)
            terminalID = ws.aoTerminalIDFromName(terminalName) ;
            self.setSingleAOChannelTerminalID(iChannel, terminalID) ;
        end
        
        function setSingleDOChannelTerminalName(self, iChannel, terminalName)
            terminalID = ws.dioTerminalIDFromName(terminalName) ;
            self.setSingleDOChannelTerminalID(iChannel, terminalID) ;
        end
        
        function result = nextFreeAITerminal(self)
            % Tries to find a free AI terminal, and returns it in a scalar struct with
            % fields deviceName and terminalID if found.  Otherwise, returns an empty
            % struct.  This is normally called when adding an AI channel.
            deviceNameForEachAIChannel = self.AIChannelDeviceNames ;
            terminalIDForEachAIChannel = self.Acquisition_.AnalogTerminalIDs ;
            allDeviceNames = self.AllDeviceNames ;
            nAITerminalsPerDevice = self.NAITerminalsPerDevice ;
            
            nAIChannels = length(deviceNameForEachAIChannel) ;
            nDevices = length(allDeviceNames) ;
            for iDevice = 1:nDevices ,
                deviceName = allDeviceNames{iDevice} ;
                nAITerminals = nAITerminalsPerDevice(iDevice) ;
                aiTerminalsOnDevice = ws.differentialAITerminalIDsGivenCount(nAITerminals) ;
                for terminalID = aiTerminalsOnDevice ,
                    isTerminalFree = true ;  % optimism!
                    for iChannel = 1:nAIChannels ,
                        channelDeviceName = deviceNameForEachAIChannel{iChannel} ;
                        channelTerminalID = terminalIDForEachAIChannel(iChannel) ;
                        if isequal(deviceName, channelDeviceName) && isequal(terminalID, channelTerminalID) ,
                            isTerminalFree = false ;
                            break
                        end
                    end
                    if isTerminalFree ,
                        result = struct('deviceName', {deviceName}, ...
                                        'terminalID', {terminalID} ) ;
                        return
                    end
                end
            end
            % if get here, no free terminal
            result = struct('deviceName', cell(1,0), ...
                            'terminalID', cell(1,0) ) ;
        end
        
        function result = nextFreeAOTerminal(self)
            % Tries to find a free AO terminal, and returns it in a scalar struct with
            % fields deviceName and terminalID if found.  Otherwise, returns an empty
            % struct.  This is normally called when adding an AO channel.
            deviceNameForEachAOChannel = self.AOChannelDeviceNames ;
            terminalIDForEachAOChannel = self.Stimulation_.AnalogTerminalIDs ;
            allDeviceNames = self.AllDeviceNames ;
            nAOTerminalsPerDevice = self.NAOTerminalsPerDevice ;
            
            nAOChannels = length(deviceNameForEachAOChannel) ;
            nDevices = length(allDeviceNames) ;
            for iDevice = 1:nDevices ,
                deviceName = allDeviceNames{iDevice} ;
                nAOTerminals = nAOTerminalsPerDevice(iDevice) ;
                for terminalID = 0:(nAOTerminals-1) ,
                    isTerminalFree = true ;  % optimism!
                    for iChannel = 1:nAOChannels ,
                        channelDeviceName = deviceNameForEachAOChannel{iChannel} ;
                        channelTerminalID = terminalIDForEachAOChannel(iChannel) ;
                        if isequal(deviceName, channelDeviceName) && isequal(terminalID, channelTerminalID) ,
                            isTerminalFree = false ;
                            break
                        end
                    end
                    if isTerminalFree ,
                        result = struct('deviceName', {deviceName}, ...
                                        'terminalID', {terminalID} ) ;
                        return
                    end
                end
            end
            % if get here, no free terminal
            result = struct('deviceName', cell(1,0), ...
                            'terminalID', cell(1,0) ) ;
        end
        
        function result = nextFreeDIOTerminal(self)
            % Tries to find a free DIO terminal, and returns it in a scalar struct with
            % fields deviceName and terminalID if found.  Otherwise, returns an empty
            % struct.  This is normally called when adding a DI or DO channel.
            deviceNameForEachDIChannel = self.DIChannelDeviceNames ;  % all elements should be equal to the primary device name
            terminalIDForEachDIChannel = self.Acquisition_.DigitalTerminalIDs ;
            deviceNameForEachDOChannel = self.DOChannelDeviceNames ;  % all elements should be equal to the primary device name
            terminalIDForEachDOChannel = self.Stimulation_.DigitalTerminalIDs ;
            deviceNameForEachDIOChannel = horzcat(deviceNameForEachDIChannel, deviceNameForEachDOChannel) ;  
                % all elements should be equal to the primary device name
            terminalIDForEachDIOChannel = horzcat(terminalIDForEachDIChannel, terminalIDForEachDOChannel) ;
            %allDeviceNames = self.AllDeviceNames ;
            nDIOTerminalsPerDevice = self.NDIOTerminalsPerDevice ;
            primaryDeviceName = self.PrimaryDeviceName ;
            
            nDIOChannels = length(deviceNameForEachDIOChannel) ;
            %nDevices = length(allDeviceNames) ;
            %for iDevice = 1:nDevices ,
            %deviceName = primaryDeviceName ;
            primaryDeviceIndexMaybe = self.getPrimaryDeviceIndexMaybe() ;
            if isempty(primaryDeviceIndexMaybe) ,
                nDIOTerminals = 0 ;
            else
                primaryDeviceIndex = primaryDeviceIndexMaybe(1) ;
                nDIOTerminals = nDIOTerminalsPerDevice(primaryDeviceIndex) ;
            end
            for terminalID = 0:(nDIOTerminals-1) ,
                isTerminalFree = true ;  % optimism!
                for iChannel = 1:nDIOChannels ,
                    channelDeviceName = deviceNameForEachDIOChannel{iChannel} ;
                    channelTerminalID = terminalIDForEachDIOChannel(iChannel) ;
                    if isequal(primaryDeviceName, channelDeviceName) && isequal(terminalID, channelTerminalID) ,
                        isTerminalFree = false ;
                        break
                    end
                end
                if isTerminalFree ,
                    result = struct('deviceName', {primaryDeviceName}, ...
                                    'terminalID', {terminalID} ) ;
                    return
                end
            end
            %end
            % if get here, no free terminal
            result = struct('deviceName', cell(1,0), ...
                            'terminalID', cell(1,0) ) ;
        end
        
        function value=get.IsReady(self)
            value=(self.DegreeOfReadiness_>0);
        end               
        
%         function encoding = encodeForHeader(self)
%             % Get the default header encoding
%             encoding = encodeForHeader@ws.Encodable(self) ;            
%             
%             % Add custom field
%             thisPropertyValue = self.getStimulusLibraryCopy() ;
%             encodingOfPropertyValue = ws.encodeAnythingForHeader(thisPropertyValue) ;
%             encoding.StimulusLibrary = encodingOfPropertyValue ;
%         end
        
        function result = get.DoUsePreferences(self)
            result = self.DoUsePreferences_ ;
        end
        
        function set.DoUsePreferences(self, rawNewValue)
            if (islogical(rawNewValue) || isnumeric(rawNewValue)) && isscalar(rawNewValue) && isfinite(rawNewValue) ,
                newValue = logical(rawNewValue) ;
                self.DoUsePreferences_ = newValue ;
            else
                error('ws:invalidPropertyValue', ...
                      'DoUsePreferences must be a scalar, and must be logical or numeric and finite') ;
            end               
        end        

        function [aoData, doData] = getStimulationData(self, indexOfEpisodeWithinRun)
            % Get the current stimulus map
            stimulusMapIndex = self.Stimulation_.getCurrentStimulusMapIndex(indexOfEpisodeWithinRun) ;
            %stimulusMap = self.getCurrentStimulusMap_(indexOfEpisodeWithinSweep);

            % Set the channel data in the tasks
            aoData = self.getAnalogChannelData_(stimulusMapIndex, indexOfEpisodeWithinRun) ;
            doData = self.getDigitalChannelData_(stimulusMapIndex, indexOfEpisodeWithinRun) ;
        end
        
    end  % public methods block
    
    methods (Access = protected)
        function changeReadiness_(self, delta)
            if ~( isnumeric(delta) && isscalar(delta) && (delta==-1 || delta==0 || delta==+1 || (isinf(delta) && delta>0) ) ),
                return
            end
                    
            newDegreeOfReadinessRaw = self.DegreeOfReadiness_ + delta ;
            self.setReadiness_(newDegreeOfReadinessRaw) ;
        end  % function        
        
        function resetReadiness_(self)
            % Used during error handling to reset model back to the ready
            % state.  (NB: But only if called via the do() method!)
            self.setReadiness_(1) ;
        end  % function        
        
        function setReadiness_(self, newDegreeOfReadinessRaw)
            %fprintf('Inside setReadiness_(%d)\n', newDegreeOfReadinessRaw) ;
            %dbstack
            isReadyBefore = self.IsReady ;
            
            self.DegreeOfReadiness_ = ...
                ws.fif(newDegreeOfReadinessRaw<=1, ...
                       newDegreeOfReadinessRaw, ...
                       1) ;
                        
            isReadyAfter = self.IsReady ;
            
            if isReadyAfter ~= isReadyBefore ,
                %fprintf('Inside setReadiness_(%d), about to broadcast UpdateReadiness\n', newDegreeOfReadinessRaw) ;
                self.broadcast('UpdateReadiness');
            end            
        end  % function                
        
        function aoDataScaledAndLimited = getAnalogChannelData_(self, stimulusMapIndex, episodeIndexWithinSweep)
            % Get info about which analog channels are in the task
            isInTaskForEachAOChannel = ~self.IsAOChannelTerminalOvercommitted ;
            %isInTaskForEachAnalogChannel = self.TheFiniteAnalogOutputTask_.IsChannelInTask ;
            nAnalogChannelsInTask = sum(isInTaskForEachAOChannel) ;

            % Calculate the signals
            if isempty(stimulusMapIndex) ,
                aoData = zeros(0,nAnalogChannelsInTask) ;
                %nChannelsWithStimulus = 0 ;
            else
                channelNamesInTask = self.AOChannelNames(isInTaskForEachAOChannel) ;
                isChannelAnalog = true(1,nAnalogChannelsInTask) ;
%                 [aoData, nChannelsWithStimulus] = ...
%                     stimulusMap.calculateSignals(self.StimulationSampleRate_, channelNamesInTask, isChannelAnalog, episodeIndexWithinSweep) ;  
                aoData = ...
                    self.Stimulation_.calculateSignalsForMap(stimulusMapIndex, ...
                                                             self.StimulationSampleRate, ...
                                                             channelNamesInTask, ...
                                                             isChannelAnalog, ...
                                                             episodeIndexWithinSweep) ;
                  % each signal of aoData is in native units
            end
            
            % % Want to return the number of scans in the stimulus data
            % nScans = size(aoData,1);
            
            % If any channel scales are problematic, deal with this
            analogChannelScales = self.AOChannelScales(isInTaskForEachAOChannel) ;  % (native units)/V
            inverseAnalogChannelScales=1./analogChannelScales;  % e.g. V/(native unit)
            sanitizedInverseAnalogChannelScales = ...
                ws.fif(isfinite(inverseAnalogChannelScales), inverseAnalogChannelScales, zeros(size(inverseAnalogChannelScales)));            

            % scale the data by the channel scales
            if isempty(aoData) ,
                aoDataScaled=aoData;
            else
                aoDataScaled=bsxfun(@times,aoData,sanitizedInverseAnalogChannelScales);
            end
            % all signals in aoDataScaled are in V
            
            % limit the data to [-10 V, +10 V]
            aoDataScaledAndLimited=max(-10,min(aoDataScaled,+10));  % also eliminates nan, sets to +10

%             % Finally, assign the stimulation data to the the relevant part
%             % of the output task
%             if isempty(self.TheFiniteAnalogOutputTask_) ,
%                 error('Adam is dumb');
%             else
%                 self.TheFiniteAnalogOutputTask_.setChannelData(aoDataScaledAndLimited) ;
%             end
        end  % function

        function doDataLimited = getDigitalChannelData_(self, stimulusMapIndex, episodeIndexWithinRun)
            %import ws.*
            
            % % Calculate the episode index
            % episodeIndexWithinRun=self.NEpisodesCompleted_+1;
            
            % Calculate the signals
            isTimedForEachDOChannel = self.IsDOChannelTimed ;
            isInTaskForEachDOChannel = isTimedForEachDOChannel & ~self.IsDOChannelTerminalOvercommitted ;
            %isInTaskForEachDigitalChannel = self.TheFiniteDigitalOutputTask_.IsChannelInTask ;
            nDigitalChannelsInTask = sum(isInTaskForEachDOChannel) ;
            if isempty(stimulusMapIndex) ,
                doData=zeros(0,nDigitalChannelsInTask);  
            else
                isChannelAnalogForEachDigitalChannelInTask = false(1,nDigitalChannelsInTask) ;
                namesOfDigitalChannelsInTask = self.DOChannelNames(isInTaskForEachDOChannel) ;                
                doData = ...
                    self.Stimulation_.calculateSignalsForMap(stimulusMapIndex, ...
                                                             self.StimulationSampleRate, ...
                                                             namesOfDigitalChannelsInTask, ...
                                                             isChannelAnalogForEachDigitalChannelInTask, ...
                                                             episodeIndexWithinRun) ;
            end
            
            % % Want to return the number of scans in the stimulus data
            % nScans = size(doData,1) ;
            
            % limit the data to {false,true}
            doDataLimited = logical(doData) ;

%             % Finally, assign the stimulation data to the the relevant part
%             % of the output task
%             self.TheFiniteDigitalOutputTask_.setChannelData(doDataLimited) ;
        end  % function        
        
        function performOneLooperIterationDuringOngoingSweep_(self, ...
                                                              timeSinceSweepStart, ...
                                                              fromRunStartTicId)
                                                    
            % Acquire data, update soft real-time outputs
            [didReadFromTasks, rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData, areTasksDone] = ...
                self.Looper_.pollAcquisition(timeSinceSweepStart, fromRunStartTicId, self.SweepDuration) ;
            
%             % DEBUG: This is for debugging purposes only!
%             if timeSinceRunStartAtStartOfData > 5 ,
%                 error('ws:fakeerror', 'Stuff went bad.  Real bad.') ;
%             end
            
            % Deal with the acquired samples
            if didReadFromTasks ,
                %self.NTimesSamplesAcquiredCalledSinceRunStart_ = self.NTimesSamplesAcquiredCalledSinceRunStart_ + 1 ;
                %self.TimeOfLastSamplesAcquired_ = timeSinceRunStartAtStartOfData ;
                nScans=size(rawAnalogData,1);
                %nChannels=size(data,2);
                %assert(nChannels == numel(expectedChannelNames));

                if (nScans>0)
                    % update the current time
                    %dt = 1/acquisitionSampleRate ;
                    %self.t_ = self.t_ + nScans*dt ;  % Note that this is the time stamp of the sample just past the most-recent sample

                    % Add data to the user cache
                    %isSweepBased = isfinite(sweepDuration) ;
                    %self.Looper_.addDataToUserCache(rawAnalogData, rawDigitalData, isSweepBased) ;
                    self.samplesAcquired_(rawAnalogData, ...
                                          rawDigitalData, ...
                                          timeSinceRunStartAtStartOfData ) ;
                end
                
                if areTasksDone ,
                    self.Looper_.completeTheOngoingSweep() ;
                    %self.acquisitionSweepComplete() ;
                end
            end                        
            
%             % We'll use this in a sanity-check
%             didAcquireNonzeroScans = (size(rawAnalogData,1)>0) ;
        end  % function
        
        function performOneRefillerIterationDuringOngoingRun_(self, isStimulationTriggerIdenticalToAcquisitionTrigger)
            % Action in a run depends on whether we are also in an
            % episode, or are in-between episodes
            % Check the finite outputs, refill them if
            % needed.
            if self.Refiller_.IsPerformingEpisode ,
                didCompleteEpisode = self.Refiller_.checkIfTasksAreDoneAndEndEpisodeIfSo() ;
                if didCompleteEpisode ,
                    self.callUserMethod_('completingEpisode');           
                end
            else
                % If we're not performing an episode, see if
                % we need to start one.
                if self.Refiller_.NEpisodesCompletedSoFarThisRun < self.Refiller_.NEpisodesPerRun ,
                    %isStimulationTriggerIdenticalToAcquisitionTrigger = self.Frontend_.isStimulationTriggerIdenticalToAcquisitionTrigger() ;
                    if isStimulationTriggerIdenticalToAcquisitionTrigger ,
                        % do nothing.
                        % if they're identical, startEpisode_()
                        % is called from the startingSweep()
                        % req-rep method.
                    else
                        try
                            self.callUserMethod_('startingEpisode') ;
                            [aoData, doData] = self.getStimulationData(self.Refiller_.NEpisodesCompletedSoFarThisRun+1) ;
                            self.Refiller_.startEpisode(aoData, doData) ;
                        catch err
                            % Something went wrong
                            self.abortOngoingRun_();
                            self.changeReadiness_(+1);
                            rethrow(err);
                        end
                    end
                end
            end
        end  % function        
        
    end  % protected methods block        
    
    methods
        function result = get.IsPerformingRun(self)
            result = self.IsPerformingRun_ ;
        end
        
        function result = get.IsPerformingSweep(self)
            result = self.IsPerformingSweep_ ;
        end
        
        function result = get.DataCacheDurationWhenContinuous(self)
            result = self.Acquisition_.DataCacheDurationWhenContinuous ;
        end
        
        function waitForRunToComplete(self)
            while self.IsPerformingRun ,
                pause(0.2) ;
            end                
        end
        
        function playAndBlock(self)
            self.play() ;
            self.waitForRunToComplete() ;
        end
        
        function recordAndBlock(self)
            self.record() ;
            self.waitForRunToComplete() ;
        end            
        
        function result = get.DoesProtocolNeedSave(self)
            result = self.DoesProtocolNeedSave_ ;
        end
    end  % public methods block
    
    methods (Static)
        function [sampleFrequency, referenceClockRate] = ...
                coerceSampleFrequencyToAllowedValue(primaryDeviceName, isPrimaryDeviceAPXIDevice, desiredSampleFrequency)
            % Make sure desiredSampleFrequency is a scalar
            defaultSampleFrequency  = 20000 ;
            if isempty(desiredSampleFrequency) ,
                temp1 = defaultSampleFrequency ;  % the default value
            else
                temp1 = desiredSampleFrequency(1) ;  % make it a scalar
            end
            
            % Make sure a double
            if isnumeric(temp1) ,
                temp2 = double(temp1) ;
            else
                temp2 = defaultSampleFrequency ;  % the default value
            end
            
            % Make sure finite
            if isfinite(temp2) ,
                sanitizedDesiredSampleFrequency = temp2 ;
            else
                sanitizedDesiredSampleFrequency = defaultSampleFrequency ;  % the default value
            end
            
            % Limit to the allowed range of sampling frequencies
            %primaryDeviceName = self.PrimaryDeviceName ;
            %isPrimaryDeviceAPXIDevice = self.IsPrimaryDeviceAPXIDevice ;
            [~, referenceClockRate] = ws.getReferenceClockSourceAndRate(primaryDeviceName, primaryDeviceName, isPrimaryDeviceAPXIDevice) ;  
                % the rate only depends on the primary device
            %referenceClockRate = self.ReferenceClockRate ;  % Hz
            desiredTimebaseTicksPerSample = referenceClockRate/sanitizedDesiredSampleFrequency ;  
            integralTimebaseTicksPerSample = floor(desiredTimebaseTicksPerSample);  % err on the side of sampling faster
            maximumTimebaseTicksPerSample = 2^32-1 ;  % Note that this sets the *minimum* frequency
            minimumTimebaseTicksPerSample = 1 ;  % Note that this sets the *maximum* frequency (Although usually this isn't achievable in practice)
              % See here:
              % http://digital.ni.com/public.nsf/allkb/4BBE1409700F6CE686256E9200652F6B
              % Entitled "What Sample Rates Is my DAQ Board Actually
              % Capable of Achieving?"
              %
            actualTimebaseTicksPerSample = ...
                ws.limit(minimumTimebaseTicksPerSample, integralTimebaseTicksPerSample, maximumTimebaseTicksPerSample) ;  % enforce min, max
              % Note that if actualTimebaseTicksPerSample is 1, then the
              % sampleFrequency is equal to the timebaseFrequency.  You can
              % set it to this in the hardware, and the board will try to
              % do it, but it will generally fail to keep up once the sweep
              % starts, because for instance the default timebase for X series cards is
              % 100 MHz.
            sampleFrequency = referenceClockRate/actualTimebaseTicksPerSample ;            
        end  % method        
    end  % static methods block
    
    methods
        function result = get.IsWavesurferMainFigureVisible(self)  %#ok<MANU>
            result = true ;
        end

        function result = get.IsGeneralSettingsFigureVisible(self)
            result = self.IsGeneralSettingsFigureVisible_ ;
        end
        
        function result = get.IsChannelsFigureVisible(self)
            result = self.IsChannelsFigureVisible_ ;
        end
        
        function result = get.IsStimulusLibraryFigureVisible(self)
            result = self.IsStimulusLibraryFigureVisible_ ;
        end
        
        function result = get.IsStimulusPreviewFigureVisible(self)
            result = self.IsStimulusPreviewFigureVisible_ ;
        end
        
        function result = get.IsTriggersFigureVisible(self)
            result = self.IsTriggersFigureVisible_ ;
        end
        
        function result = get.IsUserCodeManagerFigureVisible(self)
            result = self.IsUserCodeManagerFigureVisible_ ;
        end
        
        function result = get.IsElectrodeManagerFigureVisible(self)
            result = self.IsElectrodeManagerFigureVisible_ ;
        end
        
        function result = get.IsTestPulserFigureVisible(self)
            result = self.IsTestPulserFigureVisible_ ;
        end
        
        function result = get.IsFastProtocolsFigureVisible(self)
            result = self.IsFastProtocolsFigureVisible_ ;
        end
        
        function set.IsGeneralSettingsFigureVisible(self, newValue)
            oldValue = self.IsGeneralSettingsFigureVisible_ ;
            self.IsGeneralSettingsFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'GeneralSettings', oldValue) ;
            self.broadcast('UpdateGeneral') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsChannelsFigureVisible(self, newValue)
            oldValue = self.IsChannelsFigureVisible_ ;
            self.IsChannelsFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'Channels', oldValue) ;
            self.broadcast('UpdateChannels') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsStimulusLibraryFigureVisible(self, newValue)
            oldValue = self.IsStimulusLibraryFigureVisible_ ;
            self.IsStimulusLibraryFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'StimulusLibrary', oldValue) ;
            self.broadcast('UpdateStimulusLibrary') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsStimulusPreviewFigureVisible(self, newValue)
            oldValue = self.IsStimulusPreviewFigureVisible_ ;
            self.IsStimulusPreviewFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'StimulusPreview', oldValue) ;
            self.broadcast('UpdateStimulusPreview') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsTriggersFigureVisible(self, newValue)
            oldValue = self.IsTriggersFigureVisible_ ;
            self.IsTriggersFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'Triggers', oldValue) ;
            self.broadcast('UpdateTriggering') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsUserCodeManagerFigureVisible(self, newValue)
            oldValue = self.IsUserCodeManagerFigureVisible_ ;
            self.IsUserCodeManagerFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'UserCodeManager', oldValue) ;
            self.broadcast('UpdateUserCodeManager') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsElectrodeManagerFigureVisible(self, newValue)
            oldValue = self.IsElectrodeManagerFigureVisible_ ;
            self.IsElectrodeManagerFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'ElectrodeManager', oldValue) ;
            self.broadcast('UpdateElectrodeManager') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsTestPulserFigureVisible(self, newValue)
            oldValue = self.IsTestPulserFigureVisible_ ;
            self.IsTestPulserFigureVisible_ = newValue ;
            self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'TestPulser', oldValue) ;
            self.broadcast('UpdateTestPulser') ;
            self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function set.IsFastProtocolsFigureVisible(self, newValue)
            oldValue = self.IsFastProtocolsFigureVisible_ ;
            self.IsFastProtocolsFigureVisible_ = newValue ;
            %self.DoesProtocolNeedSave_ = self.DoesProtocolNeedSave_ || (newValue ~= oldValue) ;
            self.broadcast('DidSetSingleFigureVisibility', 'FastProtocols', oldValue) ;
            self.broadcast('UpdateFastProtocols') ;
            %self.broadcast('DidMaybeChangeProtocol') ;
        end
        
        function result = get.MainFigurePosition(self)
            result = self.MainFigurePosition_ ;
        end
        
        function result = get.GeneralSettingsFigurePosition(self)
            result = self.GeneralSettingsFigurePosition_ ;
        end
        
        function result = get.ChannelsFigurePosition(self)
            result = self.ChannelsFigurePosition_ ;
        end
        
        function result = get.StimulusLibraryFigurePosition(self)
            result = self.StimulusLibraryFigurePosition_ ;
        end
        
        function result = get.StimulusPreviewFigurePosition(self)
            result = self.StimulusPreviewFigurePosition_ ;
        end
        
        function result = get.TriggersFigurePosition(self)
            result = self.TriggersFigurePosition_ ;
        end
        
        function result = get.UserCodeManagerFigurePosition(self)
            result = self.UserCodeManagerFigurePosition_ ;
        end
        
        function result = get.ElectrodeManagerFigurePosition(self)
            result = self.ElectrodeManagerFigurePosition_ ;
        end
        
        function result = get.TestPulserFigurePosition(self)
            result = self.TestPulserFigurePosition_ ;
        end
        
        function result = get.FastProtocolsFigurePosition(self)
            result = self.FastProtocolsFigurePosition_ ;
        end        
        
        function set.MainFigurePosition(self, newValue)
            self.MainFigurePosition_ = newValue ;
        end
        
        function set.GeneralSettingsFigurePosition(self, newValue)
            self.GeneralSettingsFigurePosition_ = newValue ;
        end
        
        function set.ChannelsFigurePosition(self, newValue)
            self.ChannelsFigurePosition_ = newValue ;
        end
        
        function set.StimulusLibraryFigurePosition(self, newValue)
            self.StimulusLibraryFigurePosition_ = newValue ;
        end
        
        function set.StimulusPreviewFigurePosition(self, newValue)
            self.StimulusPreviewFigurePosition_ = newValue ;
        end
        
        function set.TriggersFigurePosition(self, newValue)
            self.TriggersFigurePosition_ = newValue ;
        end
        
        function set.UserCodeManagerFigurePosition(self, newValue)
            self.UserCodeManagerFigurePosition_ = newValue ;
        end
        
        function set.ElectrodeManagerFigurePosition(self, newValue)
            self.ElectrodeManagerFigurePosition_ = newValue ;
        end
        
        function set.TestPulserFigurePosition(self, newValue)
            self.TestPulserFigurePosition_ = newValue ;
        end
        
        function set.FastProtocolsFigurePosition(self, newValue)
            self.FastProtocolsFigurePosition_ = newValue ;
        end        
    end    
    
    methods
        % These are intended for getting/setting *public* properties.
        % I.e. they are for general use, not restricted to special cases like
        % encoding or ugly hacks.
        function result = get(self, propertyName) 
            result = self.(propertyName) ;
        end
        
        function set(self, propertyName, newValue)
            self.(propertyName) = newValue ;
        end           
    end  % public methods block
    
    methods
        function callUserObjectMethod(self, methodName, varargin)
            self.UserCodeManager_.callUserObjectMethod(methodName, varargin{:}) ;
        end
        
        function setUserObjectProperty(self, propertyName, newValue)
            self.UserCodeManager_.setUserObjectProperty(propertyName, newValue) ;
        end

        function result = getUserObjectProperty(self, propertyName)
            result = self.UserCodeManager_.getUserObjectProperty(propertyName) ;
        end
    end  % public methods block
    
    methods
        function result = get.LastProtocolFilePath(self)
            result = self.LastProtocolFilePath_ ;
        end
        
        function result = get.CurrentProfileName(self)
            result = self.CurrentProfileName_ ;            
        end

        function set.CurrentProfileName(self, newProfileName)
            profileNames = self.ProfileNames_ ;
            if ws.isString(newProfileName) && ismember(newProfileName, profileNames) ,
                self.savePreferences_(self.CurrentProfileName) ;
                try
                    self.loadPreferences_(newProfileName) ;
                catch err
                    self.broadcast('UpdateMain') ;
                    error('Unable to load preferences for new profile') ;
                end                                                        
                self.CurrentProfileName_ = newProfileName ;
            else
                self.broadcast('UpdateMain') ;
                error('ws:invalidPropertyValue', ...
                      'CurrentProfileName must be one of the ProfileNames');
            end
            self.broadcast('UpdateMain') ;
            self.broadcast('UpdateFastProtocols') ;
        end
        
        function result = get.ProfileNames(self)
            result = self.ProfileNames_ ;            
        end        

        function result = packagePreferences(self)  % constant method
            lastProtocolFilePath = self.LastProtocolFilePath_ ;
            %if self.HasUserSpecifiedProtocolFileName ,
            %    lastProtocolFilePath = self.AbsoluteProtocolFileName ;
            %else
            %    lastProtocolFilePath = '' ;
            %end
            fastProtocolsAsStruct = struct([]) ;
            for i=1:self.NFastProtocols ,
                fastProtocolsAsStruct(i).ProtocolFileName = self.FastProtocols_{i}.ProtocolFileName ;
                fastProtocolsAsStruct(i).AutoStartType = self.FastProtocols_{i}.AutoStartType ;
            end
            result = struct('LastProtocolFilePath', {lastProtocolFilePath}, ...
                            'FastProtocols', {fastProtocolsAsStruct}) ;
        end        
        
        function newProfileName = createNewProfile(self)
            self.changeReadiness_(-1) ;            
            didFindAvailableName = false ;
            for i = 1:10 ,
                if i==1 ,
                    putativeNewProfileName = 'New Profile' ;
                else
                    putativeNewProfileName = sprintf('New Profile %d', i) ;
                end                    
                %putativePreferencesFilePath = fullfile(preferencesFolderPath, sprintf('%s.mat', putativeProfileName)) ;
                profileNames = self.ProfileNames_ ;
                if ~ismember(putativeNewProfileName, profileNames) ,
                    newProfileName = putativeNewProfileName ;
                    didFindAvailableName = true ;
                    break
                end
            end
            if didFindAvailableName ,
                % Write the current preferences to the current profile
                self.savePreferences_(self.CurrentProfileName) ;
                % Change state to accord with the new profile
                self.CurrentProfileName_ = newProfileName ;
                newProfileNames = sort(horzcat(self.ProfileNames_, {newProfileName})) ;
                self.ProfileNames_ = newProfileNames ;
            else
               error('Unable to find an available name for the new profile') ; 
            end
            self.changeReadiness_(+1) ;            
            self.broadcast('UpdateMain') ;
            self.broadcast('UpdateFastProtocols') ;
        end
        
        function deleteCurrentProfile(self)
            self.changeReadiness_(-1) ;            
            profileNameToDelete = self.CurrentProfileName ;
            if isequal(profileNameToDelete, 'Default') ,
               error('Sorry, you can''t delete the default profile') ;
            end
            self.CurrentProfileName = 'Default' ;  % this will update the view
            preferencesFilePath = ws.preferencesFileNameFromProfileName(profileNameToDelete) ;
            try
                delete(preferencesFilePath) ;
            catch err
               error('Unable to delete the profile preferences file from disk') ;                                 
            end
            % If get here, the file was successfully deleted
            self.ProfileNames_ = setdiff(self.ProfileNames_, profileNameToDelete) ;
            self.changeReadiness_(+1) ;            
            self.broadcast('UpdateMain') ;            
            self.broadcast('UpdateFastProtocols') ;
        end
        
        function renameCurrentProfile(self, newProfileName)
            self.changeReadiness_(-1) ;            
            % Get current values out of self
            oldProfileName = self.CurrentProfileName ;            
            oldProfileNames = self.ProfileNames_ ;
            
            % Check for name collision
            if ismember(newProfileName, oldProfileNames) ,
                error('There is already a profile named "%s"', newProfileName) ;
            end
            
            % Delete the old-profile-name preference file from disk, if it exists
            oldPreferencesFilePath = ws.preferencesFileNameFromProfileName(oldProfileName) ;
            try 
                if exist(oldPreferencesFilePath, 'file') ,
                    delete(oldPreferencesFilePath) ;
                end
            catch err
                error('Unable to delete preferences file for profile "%s", so not renaming', oldProfileName) ;
            end            
            
            % Make sure the new name is valid by trying to write out the preferences under
            % the new name
            try
                self.savePreferences_(newProfileName) ;
            catch err
                error('%s is not an allowed profile name', newProfileName) ;
            end            
            
            % Replace the old profile name with the new in profile names
            newProfileNames = ws.renameInCellString(oldProfileNames, oldProfileName, newProfileName) ;
            
            % Commit things to self            
            % The renamed profile preferences will get written to disk in delete(), or when
            % user switches to a new profile.
            self.CurrentProfileName_ = newProfileName ;
            self.ProfileNames_ = newProfileNames ;
            
            % Finally, update the view
            self.changeReadiness_(+1) ;            
            self.broadcast('UpdateMain') ;                        
            self.broadcast('UpdateFastProtocols') ;
        end
        
    end  % public methods block
    
    methods (Access = protected)        
        function loadProfileNameAndNames_(self)
            if self.DoUsePreferences_ ,
                profileName = ws.loadLastProfileName() ;
                profileNames = ws.loadProfileNames() ;  % read from disk; these will be sorted
            else
                profileName = 'Default' ;
                profileNames = { profileName } ;
            end
            self.CurrentProfileName_ = profileName ;
            self.ProfileNames_ = profileNames ;
        end
        
        function saveLastProfileName_(self) 
            if self.DoUsePreferences_ ,
                profileName = self.CurrentProfileName_ ;
                ws.saveLastProfileName(profileName) ;
            end            
        end
        
        function savePreferences_(self, profileName)
            if self.DoUsePreferences_ ,
                preferences = self.packagePreferences() ;
                ws.saveProfilePreferences(profileName, preferences) ;
            end
            % Notify SI whether or not we're *really* saving the preferences
            self.notifyScanImageThatSavingPreferencesIfYoked_(profileName) ;
        end

        function loadPreferences_(self, profileName)
            if self.DoUsePreferences_ ,
                % Load the preferences from disk
                preferences = ws.loadProfilePreferences(profileName) ;

                % Set the state of self to match the given preferences
                lastProtocolFilePath = preferences.LastProtocolFilePath ;                
                self.LastProtocolFilePath_ = lastProtocolFilePath ; 

                % Restore the fast protocols from the profile preferences
                fastProtocolsAsStruct = preferences.FastProtocols ;
                nFastProtocolsToSet = min(self.NFastProtocols, length(fastProtocolsAsStruct)) ;
                for i = 1:nFastProtocolsToSet ,
                    self.FastProtocols_{i}.setPropertyValue_('ProtocolFileName_', fastProtocolsAsStruct(i).ProtocolFileName) ;
                    self.FastProtocols_{i}.setPropertyValue_('AutoStartType_'   , fastProtocolsAsStruct(i).AutoStartType   ) ;
                end
            end
            % Notify SI whether or not we're *really* loading the preferences
            self.notifyScanImageThatLoadingPreferencesIfYoked_(profileName) ;
        end        
        
        function addData_(self, t, recentScaledAnalogData, recentRawDigitalData)
            % t is a scalar, the time stamp of the scan *just after* the
            % most recent scan.  (I.e. it is one dt==1/fs into the future.
            % Queue Doctor Who music.)
            
            nActiveDIChannels = self.getNActiveDIChannels() ;
            sampleRate = self.AcquisitionSampleRate ;
            xSpan = self.XSpan ;
            self.Display_.addData(t, recentScaledAnalogData, recentRawDigitalData, nActiveDIChannels, sampleRate, xSpan) ;
            
            self.setYAxisLimitsInModelTightToDataIfAreYLimitsLockedTightToData_() ;            
            %plotIndicesNeedingYLimitUpdate = self.PlotIndexFromChannelIndex(indicesOfAIChannelsNeedingYLimitUpdate) ;
            
            self.broadcast('UpdateAfterDataAdded') ;
        end  % function        
        
        function indicesOfAIChannelsNeedingYLimitUpdate = setYAxisLimitsInModelTightToDataIfAreYLimitsLockedTightToData_(self)
            isChannelDisplayed = self.IsAIChannelDisplayed ;
            areYLimitsLockedTightToData = self.AreYLimitsLockedTightToDataForAIChannel ;
            doesAIChannelNeedYLimitUpdate = isChannelDisplayed & areYLimitsLockedTightToData ;
            indicesOfAIChannelsNeedingYLimitUpdate = find(doesAIChannelNeedYLimitUpdate) ;
            plotIndexFromChannelIndex = self.PlotIndexFromChannelIndex ;  % for AI channels, the channel index is equal to the AI channel index
            plotIndicesOfAIChannelsNeedingYLimitUpdate = plotIndexFromChannelIndex(indicesOfAIChannelsNeedingYLimitUpdate) ;
            for i = 1:length(plotIndicesOfAIChannelsNeedingYLimitUpdate) ,
                %channelIndex = indicesOfAIChannelsNeedingYLimitUpdate(i) ;
                plotIndex = plotIndicesOfAIChannelsNeedingYLimitUpdate(i) ;
                %self.setYAxisLimitsInModelTightToData_(plotIndex, channelIndex) ;
                self.setYLimitsTightToDataForSinglePlot(plotIndex) ;
            end               
        end  % function        
        
        function syncTraces_(self)
            scaledAnalogData = self.getAIDataFromCache() ;
            [digitalDataAsUint, cachedDigitalSignalCount] = self.getDIDataFromCache() ;
            t = self.getTimestampsForDataInCache() ;
            xSpan = self.XSpan ;
            sampleRate = self.AcquisitionSampleRate ;
            self.Display_.updateTraces(scaledAnalogData, digitalDataAsUint, cachedDigitalSignalCount, t, xSpan, sampleRate) ;
            self.setYAxisLimitsInModelTightToDataIfAreYLimitsLockedTightToData_() ;            
        end
    end  % protected methods block
    
    methods
        function set.WidthOfPlotsInPixels(self, newValue)
            self.Display_.XSpanInPixels = newValue ;
            self.syncTraces_() ;            
            self.broadcast('UpdateTraces') ;  % Is this what we want, really?
        end
        
        function result = get.XDataForDisplay(self)
            result = self.Display_.XData ;
        end
        
        function result = get.YDataForDisplay(self)
            result = self.Display_.YData ;
        end
        
        function yMinAndMax = plottedDataYMinAndMax(self, aiChannelIndex)
            % Min and max of the data, across all plotted channels.
            % Returns a 1x2 array.
            % If all channels are empty, returns [+inf -inf].
            cacheChannelIndexFromChannelIndex = self.CacheInputChannelIndexFromInputChannelIndex ;
            indexWithinData = cacheChannelIndexFromChannelIndex(aiChannelIndex) ;
            if isfinite(indexWithinData) ,
                yData = self.Display_.YData ;
                y = yData(:,indexWithinData) ;
                yMinRaw = min(y) ;
                yMin = ws.fif(isempty(yMinRaw),+inf,yMinRaw) ;
                yMaxRaw = max(y) ;
                yMax = ws.fif(isempty(yMaxRaw),-inf,yMaxRaw) ;
                yMinAndMax = double([yMin yMax]) ;
            else
                yMinAndMax = [-inf +inf] ;
            end
        end                        
        
        function setAreYLimitsLockedTightToDataForSinglePlot(self, plotIndex) 
            channelIndex = self.ChannelIndexWithinTypeFromPlotIndex(plotIndex) ;
            currentValue = self.AreYLimitsLockedTightToDataForAIChannel(channelIndex) ;
            newValue = ~currentValue ;
            self.setAreYLimitsLockedTightToDataForSingleAIChannel_(channelIndex, newValue) ;
            if newValue ,
                self.setYLimitsTightToDataForSinglePlot(plotIndex) ;
            end
            self.broadcast('UpdateMain') ;  % need to update the button
        end  % method       
        
        function setYLimitsTightToDataForSinglePlot(self, plotIndex) 
            aiChannelIndex = self.ChannelIndexWithinTypeFromPlotIndex(plotIndex) ;
            yMinAndMax = self.plottedDataYMinAndMax(aiChannelIndex) ;
            if any( ~isfinite(yMinAndMax) ) ,
                return
            end
            yCenter = mean(yMinAndMax) ;
            yRadius = 0.5*diff(yMinAndMax) ;
            if yRadius == 0 ,
                yRadius = 0.001 ;
            end
            newYLimits = yCenter + 1.05*yRadius*[-1 +1] ;
            self.setYLimitsForSinglePlot(plotIndex, newYLimits) ;  % this will broadcast          
        end  % method          
        
        function result = get.WidthOfPlotsInPixels(self)
            result = self.Display_.XSpanInPixels ;
        end
    end  % public methods block
    
    methods (Access = protected)        
        function callUserMethod_(self, eventName, varargin)
            % Handle user functions.  It would be possible to just make the UserCodeManager
            % subsystem a regular listener of these events.  Handling it
            % directly removes at 
            % least one layer of function calls and allows for user functions for 'events'
            % that are not formally events on the model.
            self.UserCodeManager_.invoke(self, eventName, varargin{:});
            
            % Handle as standard event if applicable.
            %self.broadcast(eventName);
        end  % function     
    end  % protected methods block    
end  % classdef
