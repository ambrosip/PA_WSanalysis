classdef DOTask < handle
    properties (Access = protected)
        SampleRate_ = 20000
        ChannelCount_ = 0
        ChannelData_
        IsOutputBufferSyncedToChannelData_ = false
        DAQmxTaskHandle_ = []  % Can be empty if there are zero channels
    end
    
    methods
        function self = DOTask(taskName, primaryDeviceName, isPrimaryDeviceAPXIDevice, terminalIDs, ...
                               sampleRate, ...
                               keystoneTaskType, keystoneTaskDeviceName, ...
                               triggerDeviceNameIfKeystone, triggerPFIIDIfKeystone, triggerEdgeIfKeystone)
                                    
            % Create the channels, set the timing mode (has to be done
            % after adding channels)
            nChannels=length(terminalIDs);
            if nChannels>0 ,
                %self.DAQmxTaskHandle_ = ws.dabs.ni.daqmx.Task(taskName);
                self.DAQmxTaskHandle_ = ws.ni('DAQmxCreateTask', taskName) ;
                for i=1:nChannels ,
                    deviceName = primaryDeviceName ;
                    terminalID = terminalIDs(i) ;                    
                    %lineName = sprintf('line%d',terminalID) ;                        
                    %self.DAQmxTaskHandle_.createDOChan(deviceName, lineName) ;
                    lineSpecification = sprintf('%s/line%d', deviceName, terminalID) ;
                    ws.ni('DAQmxCreateDOChan', self.DAQmxTaskHandle_, lineSpecification, 'DAQmx_Val_ChanForAllLines')
                end
                [referenceClockSource, referenceClockRate] = ...
                    ws.getReferenceClockSourceAndRate(primaryDeviceName, primaryDeviceName, isPrimaryDeviceAPXIDevice) ;                
                %set(self.DAQmxTaskHandle_, 'refClkSrc', referenceClockSource) ;                
                %set(self.DAQmxTaskHandle_, 'refClkRate', referenceClockRate) ;                
                ws.ni('DAQmxSetRefClkSrc', self.DAQmxTaskHandle_, referenceClockSource) ;
                ws.ni('DAQmxSetRefClkRate', self.DAQmxTaskHandle_, referenceClockRate) ;
                %self.DAQmxTaskHandle_.cfgSampClkTiming(sampleRate, 'DAQmx_Val_FiniteSamps');
                source = '';  % means to use default sample clock
                scanCount = 2 ;
                  % 2 is the minimum value advisable when using FiniteSamps mode. If this isn't
                  % set, then Error -20077 can occur if writeXXX() operation precedes configuring
                  % sampQuantSampPerChan property to a non-zero value -- a bit strange,
                  % considering that it is allowed to buffer/write more data than specified to
                  % generate.
                ws.ni('DAQmxCfgSampClkTiming', self.DAQmxTaskHandle_, source, sampleRate, 'DAQmx_Val_Rising', 'DAQmx_Val_FiniteSamps', scanCount);
                try
                    %self.DAQmxTaskHandle_.control('DAQmx_Val_Task_Verify');
                    ws.ni('DAQmxTaskControl', self.DAQmxTaskHandle_, 'DAQmx_Val_Task_Verify');
                catch me
                    error('There was a problem setting up the finite output task');
                end
                try
                    %taskSampleClockRate = self.DAQmxTaskHandle_.sampClkRate ;
                    taskSampleClockRate = ws.ni('DAQmxGetSampClkRate', self.DAQmxTaskHandle_) ;
                catch cause
                    rawException = MException('ws:errorGettingTaskSampleClockRate', ...
                                              'There was a problem getting the sample clock rate of the DAQmx task in order to check it') ;
                    exception = addCause(rawException, cause) ;
                    throw(exception) ;
                end                                    
                if taskSampleClockRate ~= sampleRate ,
                    error('The DAQmx task sample rate is not equal to the desired sampling rate');
                end
                
                % Figure out the trigger terminal and edge
                if isequal(keystoneTaskType,'ai') ,
                    triggerTerminalName = sprintf('/%s/ai/StartTrigger', keystoneTaskDeviceName) ;
                    triggerEdge = 'rising' ;
                elseif isequal(keystoneTaskType,'di') ,
                    triggerTerminalName = sprintf('/%s/di/StartTrigger', keystoneTaskDeviceName) ;
                    triggerEdge = 'rising' ;
                elseif isequal(keystoneTaskType,'ao') ,
                    triggerTerminalName = sprintf('/%s/ao/StartTrigger', keystoneTaskDeviceName) ;
                    triggerEdge = 'rising' ;
                elseif isequal(keystoneTaskType,'do') ,
                    triggerTerminalName = sprintf('/%s/PFI%d', triggerDeviceNameIfKeystone, triggerPFIIDIfKeystone) ;
                    triggerEdge = triggerEdgeIfKeystone ;
                else
                    % This is here mostly to allow for easier testing
                    triggerTerminalName = '' ;
                    triggerEdge = [] ;
                end
                
                % Set up triggering
%                 if isempty(triggerTerminalName) ,
%                     self.DAQmxTaskHandle_.disableStartTrig() ;
%                 else
%                     dabsTriggerEdge = ws.dabsEdgeTypeFromEdgeType(triggerEdge) ;
%                     self.DAQmxTaskHandle_.cfgDigEdgeStartTrig(triggerTerminalName, dabsTriggerEdge) ;
%                 end                
                if isempty(triggerTerminalName) ,
                    % This is mostly here for testing
                    ws.ni('DAQmxDisableStartTrig', self.DAQmxTaskHandle_) ;
                else
                    daqmxTriggerEdge = ws.daqmxEdgeTypeFromEdgeType(triggerEdge) ;
                    ws.ni('DAQmxCfgDigEdgeStartTrig', self.DAQmxTaskHandle_, triggerTerminalName, daqmxTriggerEdge);
                end                
            else
                % if no channels
                self.DAQmxTaskHandle_ = [] ;
            end        
            
            % Store stuff
            %self.PrimaryDeviceName_ = primaryDeviceName ;
            self.ChannelCount_ = nChannels ;
            self.SampleRate_ = sampleRate ;
            
            % Init the buffer
            self.clearChannelData() ;
        end  % function
        
        function delete(self)
            %ws.deleteIfValidHandle(self.DAQmxTaskHandle_) ;  % have to explicitly delete, b/c ws.dabs.ni.daqmx.System has refs to, I guess
            if ~isempty(self.DAQmxTaskHandle_) ,
                if ~ws.ni('DAQmxIsTaskDone', self.DAQmxTaskHandle_) ,
                    ws.ni('DAQmxStopTask', self.DAQmxTaskHandle_) ;
                end
                ws.ni('DAQmxClearTask', self.DAQmxTaskHandle_) ;
            end            
            self.DAQmxTaskHandle_ = [] ;  % not really necessary...
        end  % function
        
        function start(self)
            if ~isempty(self.DAQmxTaskHandle_) ,
                %self.DAQmxTaskHandle_.start();
                ws.ni('DAQmxStartTask', self.DAQmxTaskHandle_) ;
            end
        end  % function
        
        function stop(self)
            if ~isempty(self.DAQmxTaskHandle_) ,
                %self.DAQmxTaskHandle_.stop();
                ws.ni('DAQmxStopTask', self.DAQmxTaskHandle_) ;
            end
        end  % function
        
        function disableTrigger(self)
            % This is used in the process of setting all outputs to zero when stopping
            % a run.  The DOTask is mostly useless after, because no was to restore the
            % triggers.
            if ~isempty(self.DAQmxTaskHandle_) ,
                %self.DAQmxTaskHandle_.disableStartTrig() ;
                ws.ni('DAQmxDisableStartTrig', self.DAQmxTaskHandle_) ;
            end
        end
        
        function clearChannelData(self)
            % This is used to just get rid of any pre-existing channel
            % data.  Typically used at the start of a run, to clear out any
            % old channel data.              
            nChannels = self.ChannelCount_ ;
            self.ChannelData_ = false(0,nChannels);
            self.IsOutputBufferSyncedToChannelData_ = false ;  % we don't sync up the output buffer to no data
        end  % function
        
        function zeroChannelData(self)
            % This is used to replace the channel data with a small number
            % of all-zero scans.
            nChannels = self.ChannelCount_ ;
            nScans = 2 ;  % Need at least 2...
            self.setChannelData(false(nScans,nChannels)) ;  % N.B.: Want to use public setter, so output gets sync'ed
        end  % function       
        
        function setChannelData(self, value)
            nChannels = self.ChannelCount_ ;
            requiredType = 'logical' ;
            if isa(value,requiredType) && ismatrix(value) && (size(value,2)==nChannels) ,
                self.ChannelData_ = value;
                self.IsOutputBufferSyncedToChannelData_ = false ;
                self.syncOutputBufferToChannelData_();
            else
                error('ws:invalidPropertyValue', ...
                      'ChannelData must be an NxR matrix, R the number of channels, of the appropriate type.');
            end
        end  % function        
        
        function debug(self) %#ok<MANU>
            keyboard
        end  % function        
        
        function result = isDone(self)
            if isempty(self.DAQmxTaskHandle_) ,
                % This means there are no channels, so nothing to do
                result = true ;  % things work out better if you use this convention
            else
                %result = self.DAQmxTaskHandle_.isTaskDoneQuiet() ;
                result = ws.ni('DAQmxIsTaskDone', self.DAQmxTaskHandle_) ;
            end
        end  % function
    end  % public methods
    
    methods (Access = protected)
        function syncOutputBufferToChannelData_(self)
            % If already up-to-date, do nothing
            if self.IsOutputBufferSyncedToChannelData_ ,
                return
            end
            
            % Actually set up the task, if present
            if isempty(self.DAQmxTaskHandle_) ,
                % do nothing
            else            
                % Get the channel data into a local
                channelData = self.ChannelData_ ;

                % The outputData is the channelData, unless the channelData is
                % very short, in which case the outputData is just long
                % enough, and all zeros
                nScansInChannelData = size(channelData,1) ;            
                if nScansInChannelData<2 ,
                    nChannels = self.ChannelCount_ ;
                    outputData = false(2,nChannels) ;
                else
                    outputData = self.ChannelData_ ;
                end
                
                % Resize the output buffer to the number of scans in outputData
                %self.DAQmxTaskHandle_.control('DAQmx_Val_Task_Unreserve');  
                % this is needed b/c we might have filled the buffer before bur never outputed that data
                nScansInOutputData = size(outputData,1) ;
                %nScansInBuffer = self.DAQmxTaskHandle_.get('bufOutputBufSize') ;
                nScansInBuffer = ws.ni('DAQmxGetBufOutputBufSize', self.DAQmxTaskHandle_) ;
                if nScansInBuffer ~= nScansInOutputData ,
                    %self.DAQmxTaskHandle_.cfgOutputBuffer(nScansInOutputData) ;
                    ws.ni('DAQmxCfgOutputBuffer', self.DAQmxTaskHandle_, nScansInOutputData) ;
                end

                % Configure the the number of scans in the finite-duration output
                sampleRate = self.SampleRate_ ;
                %self.DAQmxTaskHandle_.cfgSampClkTiming(sampleRate, 'DAQmx_Val_FiniteSamps', nScansInOutputData) ;
                ws.ni('DAQmxCfgSampClkTiming', self.DAQmxTaskHandle_, '', sampleRate, 'DAQmx_Val_Rising', 'DAQmx_Val_FiniteSamps', nScansInOutputData);
                  % We validated the sample rate when we created the
                  % FiniteOutputTask, so this should be OK, but check
                  % anyway.
                sampleRateAsRead = ws.ni('DAQmxGetSampClkRate', self.DAQmxTaskHandle_) ;  
                if sampleRateAsRead ~= sampleRate ,
                    error('The DAQmx task sample rate is not equal to the desired sampling rate');
                end

                % Write the data to the output buffer
                outputData(end,:) = false ;  % don't want to end on nonzero value
                %self.DAQmxTaskHandle_.reset('writeRelativeTo');
                %self.DAQmxTaskHandle_.reset('writeOffset');
                ws.ni('DAQmxResetWriteRelativeTo', self.DAQmxTaskHandle_) ;
                ws.ni('DAQmxResetWriteOffset', self.DAQmxTaskHandle_) ;
                %self.DAQmxTaskHandle_.writeDigitalData(outputData) ;
                autoStart = false ;  % Don't automatically start the task.  This is typically what you want for a timed task.
                timeout = -1 ;  % wait indefinitely
                ws.ni('DAQmxWriteDigitalLines', self.DAQmxTaskHandle_, autoStart, timeout, outputData) ;
            end
            
            % Note that we are now synched
            self.IsOutputBufferSyncedToChannelData_ = true ;
        end  % function
    end  % protected methods block    
end  % classdef
