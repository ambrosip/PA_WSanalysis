classdef Stimulation < ws.Model
    % Stimulation subsystem
    
    properties (Dependent = true)
        IsEnabled
    end
    
    properties (Access = protected)
        IsEnabled_ = false
    end
    
    properties (Dependent = true)
        DoRepeatSequence  % should really be named DoRepeatOutputable, since it applies to 'naked' maps also
        AnalogTerminalIDs
        DigitalTerminalIDs
    end
    
    properties (Dependent = true, SetAccess = immutable)  % N.B.: it's not settable, but it can change over the lifetime of the object
        AnalogDeviceNames
        DigitalDeviceNames
        AnalogTerminalNames % the physical channel name for each analog channel, e.g. 'AO0'
        DigitalTerminalNames  % the physical channel name for each digital channel, e.g. 'line0'
        TerminalNames
        AnalogChannelNames
        DigitalChannelNames
        ChannelNames
        NAnalogChannels
        NDigitalChannels
        %NTimedDigitalChannels        
        NChannels
        IsChannelAnalog
    end
    
    properties (Access = protected)
        SampleRate_ = 20000  % Hz
        AnalogDeviceNames_ = cell(1,0)  % the device name for each analog channel
        DigitalDeviceNames_ = cell(1,0)  % the device name for each digital channel        
        AnalogChannelNames_ = cell(1,0)  % the (user) channel name for each analog channel
        DigitalChannelNames_ = cell(1,0)  % the (user) channel name for each digital channel        
        AnalogChannelScales_ = zeros(1,0)  % Store for the current AnalogChannelScales values, but values may be "masked" by ElectrodeManager
        AnalogChannelUnits_ = cell(1,0)  % Store for the current AnalogChannelUnits values, but values may be "masked" by ElectrodeManager
        StimulusLibrary_ 
        DoRepeatSequence_ = true  % If true, the stimulus sequence will be repeated ad infinitum
        IsDigitalChannelTimed_ = false(1,0)
        DigitalOutputStateIfUntimed_ = false(1,0)
        AnalogTerminalIDs_ = zeros(1,0)
        DigitalTerminalIDs_ = zeros(1,0)
    end
    
    properties (Access = protected, Transient=true)
        IsAnalogChannelMarkedForDeletion_ = false(1,0)
        IsDigitalChannelMarkedForDeletion_ = false(1,0)        
    end    
    
    methods
        function self = Stimulation()
            %self@ws.Subsystem() ;
            self.StimulusLibrary_ = ws.StimulusLibrary();  % create a StimulusLibrary, which doesn't need to know its parent
        end
        
        function startingRun(self) %#ok<MANU>
        end

        function debug(self) %#ok<MANU>
            keyboard
        end
        
    end  % methods block
        
    methods 
        % Allows access to protected and protected variables from ws.Encodable.
        function out = getPropertyValue_(self, name)
            out = self.(name);
        end
        
        % Allows access to protected and protected variables from ws.Encodable.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end        
    end  % protected methods block
    
%     methods (Access = protected)
%         function syncTasksToChannelMembership_(self)             %#ok<MANU>
%             % Clear the timed digital output task, will be recreated when acq is
%             % started.  Have to do this b/c the channels used for the timed digital output task has changed.
%             % And have to do it first to avoid a temporary collision.
%             
%             % Subclasses override this as appropriate
%         end  % function        
%     end  % protected methods

    methods 
        function wasSet = setSingleDigitalTerminalID_(self, i, newValue)
            % This should only be called by the parent WavesurferModel, to
            % ensure the self-consistency of the WavesurferModel.
            %wasSet = setSingleDigitalTerminalID_@ws.StimulationSubsystem(self, i, newValue) ;            
            if 1<=i && i<=self.NDigitalChannels && isnumeric(newValue) && isscalar(newValue) && isfinite(newValue) ,
                newValueAsDouble = double(newValue) ;
                if newValueAsDouble>=0 && newValueAsDouble==round(newValueAsDouble) ,
                    self.DigitalTerminalIDs_(i) = newValueAsDouble ;
                    wasSet = true ;
                else
                    wasSet = false ;
                end
            else
                wasSet = false ;
            end
%             self.Parent.didSetDigitalOutputTerminalID() ;
%             if wasSet ,
%                 self.Parent.singleDigitalOutputTerminalIDWasSetInStimulationSubsystem(i) ;
%             end
        end
        
        function newChannelIndex = addDigitalChannel(self, deviceNameForNewChannel, newTerminalID)
            newChannelName = sprintf('P0.%d',newTerminalID) ;
            
            self.DigitalDeviceNames_ = [self.DigitalDeviceNames_ {deviceNameForNewChannel} ] ;
            self.DigitalTerminalIDs_ = [self.DigitalTerminalIDs_ newTerminalID] ;
            self.DigitalChannelNames_ = [self.DigitalChannelNames_ {newChannelName}] ;
            self.IsDigitalChannelTimed_ = [  self.IsDigitalChannelTimed_ true  ];
            self.DigitalOutputStateIfUntimed_ = [  self.DigitalOutputStateIfUntimed_ false ];
            self.IsDigitalChannelMarkedForDeletion_ = [  self.IsDigitalChannelMarkedForDeletion_ false ];
            
            newChannelIndex = length(self.DigitalChannelNames_) ;
        end  % function
        
        function deleteMarkedDigitalChannels_(self, isToBeDeleted)
            % Should only be called from parent.
            
%             % Do some accounting
%             isToBeDeleted = self.IsDigitalChannelMarkedForDeletion_ ;
%             %indicesOfChannelsToDelete = find(isToBeDeleted) ;
           isKeeper = ~isToBeDeleted ;
%             
%             % Turn off any untimed DOs that are about to be deleted
%             digitalOutputStateIfUntimed = self.DigitalOutputStateIfUntimed_ ;
%             self.DigitalOutputStateIfUntimed = digitalOutputStateIfUntimed & isKeeper ;

            % Now do the real deleting
            if all(isToBeDeleted)
                % Keep everything a row vector
                self.DigitalDeviceNames_ = cell(1,0) ;
                self.DigitalTerminalIDs_ = zeros(1,0) ;
                self.DigitalChannelNames_ = cell(1,0) ;
                self.IsDigitalChannelTimed_ = true(1,0) ;
                self.DigitalOutputStateIfUntimed_ = false(1,0) ;
                self.IsDigitalChannelMarkedForDeletion_ = false(1,0) ;                
            else
                self.DigitalDeviceNames_ = self.DigitalDeviceNames_(isKeeper) ;
                self.DigitalTerminalIDs_ = self.DigitalTerminalIDs_(isKeeper) ;
                self.DigitalChannelNames_ = self.DigitalChannelNames_(isKeeper) ;
                self.IsDigitalChannelTimed_ = self.IsDigitalChannelTimed_(isKeeper) ;
                self.DigitalOutputStateIfUntimed_ = self.DigitalOutputStateIfUntimed_(isKeeper) ;
                self.IsDigitalChannelMarkedForDeletion_ = self.IsDigitalChannelMarkedForDeletion_(isKeeper) ;
            end
            %self.syncIsDigitalChannelTerminalOvercommitted_() ;

%             % Notify others of what we have done
%             self.Parent.didDeleteDigitalOutputChannels() ;
%             self.notifyLibraryThatDidChangeNumberOfOutputChannels_() ;
        end  % function
       
        function wasSet = setIsDigitalChannelTimed_(self, newValue)
            nDigitalChannels = length(self.IsDigitalChannelTimed_) ;
            if isequal(size(newValue),[1 nDigitalChannels]) && (islogical(newValue) || (isnumeric(newValue) && ~any(isnan(newValue)))) ,
                coercedNewValue = logical(newValue) ;
                if any(self.IsDigitalChannelTimed_ ~= coercedNewValue) ,
                    self.IsDigitalChannelTimed_=coercedNewValue;
                    wasSet = true ;
                else
                    wasSet = false ;
                end
            else
                %self.Parent.didSetIsDigitalOutputTimed();
                error('ws:invalidPropertyValue', ...
                      'IsDigitalChannelTimed must be a logical 1x%d vector, or convertable to one',nDigitalChannels);
            end
%             self.Parent.didSetIsDigitalOutputTimed();
%             
%             if wasSet ,
%                 self.Parent.isDigitalChannelTimedWasSetInStimulationSubsystem() ;
%             end  
        end  % function        
    end  % public methods
        
    methods
        function setDigitalOutputStateIfUntimed_(self, newValue)
            if isequal(size(newValue),size(self.DigitalOutputStateIfUntimed_)) && ...
                    (islogical(newValue) || (isnumeric(newValue) && ~any(isnan(newValue)))) ,
                coercedNewValue = logical(newValue) ;
                self.DigitalOutputStateIfUntimed_ = coercedNewValue ;
            else
                %self.Parent.didSetDigitalOutputStateIfUntimed() ;
                error('ws:invalidPropertyValue', ...
                      'DigitalOutputStateIfUntimed must be a logical row vector, or convertable to one, of the proper size');
            end
            %self.Parent.didSetDigitalOutputStateIfUntimed() ;            
            %self.Parent.digitalOutputStateIfUntimedWasSetInStimulationSubsystem() ;
        end  % function
    end  % protected methods

%     methods (Access=protected)    
%         function disableAllBroadcastsDammit_(self)
%             self.disableBroadcasts() ;
%             self.StimulusLibrary_.disableBroadcasts() ;
%         end
%         
%         function enableBroadcastsMaybeDammit_(self)
%             self.StimulusLibrary_.enableBroadcastsMaybe() ;            
%             self.enableBroadcastsMaybe() ;
%         end
%     end  % protected methods block
    
    methods
%         function self = StimulationSubsystem(parent)
%             self@ws.Subsystem(parent) ;            
%             self.StimulusLibrary_ = ws.StimulusLibrary(self);  % create a StimulusLibrary
%         end
        
%         function value=get.StimulusLibrary(self)
%             value=self.StimulusLibrary_;
%         end
%         
%         function set.StimulusLibrary(self,newValue)
%             if isempty(newValue) ,
%                 if isempty(self.StimulusLibrary_) ,
%                     % do nothing
%                 else
%                     self.StimulusLibrary_ = [] ;
%                 end
%             elseif isa(newValue, 'ws.StimulusLibrary') && isscalar(newValue) ,
%                 if isempty(self.StimulusLibrary_) || self.StimulusLibrary_ ~= newValue ,
%                     self.StimulusLibrary_ = newValue.copy() ;
%                     %self.StimulusLibrary_.Parent = self ;
%                 end
%             end
%             %self.broadcast('DidSetStimulusLibrary');
%         end

        function result = getStimulusLibraryCopy(self)
            result = ws.copy(self.StimulusLibrary_) ;
        end

        function out = getSampleRate_(self)
            out= self.SampleRate_ ;
        end
        
        function setSampleRate_(self, newValue)
            self.SampleRate_ = newValue ;
        end  % function
                
        function out = getIsDigitalChannelTimed_(self)
            out= self.IsDigitalChannelTimed_ ;
        end
        
        function out = getDigitalOutputStateIfUntimed_(self)
            out= self.DigitalOutputStateIfUntimed_ ;
        end
        
        function out = get.DoRepeatSequence(self)
            out= self.DoRepeatSequence_ ;
        end
        
        function set.DoRepeatSequence(self, value)
            if (islogical(value) || isnumeric(value)) && isscalar(value) ,
                self.DoRepeatSequence_ = logical(value);
            end
            %self.broadcast('DidSetDoRepeatSequence');
        end
        
        function result = get.AnalogTerminalNames(self)
            %result = self.AnalogTerminalNames_ ;
            terminalIDs = self.AnalogTerminalIDs_ ;
            function name = terminalNameFromID(id)
                name = sprintf('AO%d',id);
            end            
            result = arrayfun(@terminalNameFromID,terminalIDs,'UniformOutput',false);
        end
    
        function result = get.DigitalTerminalNames(self)
            %result = self.DigitalTerminalNames_ ;
            terminalIDs = self.DigitalTerminalIDs_ ;
            function name = terminalNameFromID(id)
                name = sprintf('P0.%d',id);
            end            
            result = arrayfun(@terminalNameFromID,terminalIDs,'UniformOutput',false);
        end

        function result = get.TerminalNames(self)
            result = [self.AnalogTerminalNames self.DigitalTerminalNames] ;
        end
        
        function result = get.AnalogChannelNames(self)
            result = self.AnalogChannelNames_ ;
        end
    
        function result = get.DigitalChannelNames(self)
            result = self.DigitalChannelNames_ ;
        end
    
        function result = get.ChannelNames(self)
            result = [self.AnalogChannelNames self.DigitalChannelNames] ;
        end
    
%         function result = get.DeviceNamePerAnalogChannel(self)
%             result = ws.deviceNamesFromTerminalNames(self.AnalogTerminalNames);
%         end
                
        function value = get.NAnalogChannels(self)
            value = length(self.AnalogChannelNames_);
        end
        
        function value = get.NDigitalChannels(self)
            value = length(self.DigitalChannelNames_);
        end

%         function value = get.NTimedDigitalChannels(self)
%             value = sum(self.IsDigitalChannelTimed);
%         end

        function value = get.NChannels(self)
            value = self.NAnalogChannels + self.NDigitalChannels ;
        end
        
        function value=isAnalogChannelName(self,name)
            value=any(strcmp(name,self.AnalogChannelNames));
        end
        
        function value = get.IsChannelAnalog(self)
            % Boolean array indicating, for each channel, whether is is analog or not
            value = [true(1,self.NAnalogChannels) false(1,self.NDigitalChannels)];
        end
        
%         function output = get.TriggerScheme(self)
%             triggering = self.Parent.getTriggeringEvenThoughThisIsDangerous() ;
%             output = triggering.StimulationTriggerScheme ;
%         end
        
%         function output = get.DeviceNames(self)
%             output = [self.AnalogDeviceNames_ self.DigitalDeviceNames_] ;
%         end
        
%         function out = get.AnalogDeviceNames(self)
%             %out = self.AnalogDeviceNames_ ;
%             deviceName = self.Parent.DeviceName ;
%             out = repmat({deviceName}, size(self.AnalogChannelNames)) ;             
%         end  % function
        
%         function digitalDeviceNames = get.DigitalDeviceNames(self)
%             %out = self.DigitalDeviceNames_ ;
%             deviceName = self.Parent.DeviceName ;
%             digitalDeviceNames = repmat({deviceName}, size(self.DigitalChannelNames)) ;
%         end  % function

        function result = get.AnalogTerminalIDs(self)
            result = self.AnalogTerminalIDs_;
        end
        
        function result = get.DigitalTerminalIDs(self)
            result = self.DigitalTerminalIDs_;
        end
        
        function result=getIsAnalogChannelMarkedForDeletion_(self)
            result =  self.IsAnalogChannelMarkedForDeletion_ ;
        end
        
        function setIsAnalogChannelMarkedForDeletion_(self, newValue)
            if islogical(newValue) && isequal(size(newValue),size(self.IsAnalogChannelMarkedForDeletion_)) ,
                self.IsAnalogChannelMarkedForDeletion_ = newValue;
            end
            %self.Parent.didSetIsInputChannelMarkedForDeletion() ;
        end
        
        function result=getIsDigitalChannelMarkedForDeletion_(self)
            result = self.IsDigitalChannelMarkedForDeletion_ ;
        end
        
        function setIsDigitalChannelMarkedForDeletion_(self, newValue)
            if islogical(newValue) && isequal(size(newValue),size(self.IsDigitalChannelMarkedForDeletion_)) ,
                self.IsDigitalChannelMarkedForDeletion_ = newValue;
            end
            %self.Parent.didSetIsInputChannelMarkedForDeletion() ;
        end
    end  % methods block
    
    methods
        function terminalID=analogTerminalIDFromName(self,channelName)
            % Get the channel ID, given the name.
            % This returns a channel ID, e.g. if the channel is ao2,
            % it returns 2.
            channelIndex = self.aoChannelIndexFromName(channelName) ;
            if isnan(channelIndex) ,
                terminalID = nan ;
            else
                terminalID = self.AnalogTerminalIDs_(channelIndex) ;
                %terminalName = self.AnalogTerminalNames_{iChannel};
                %terminalID = ws.terminalIDFromTerminalName(terminalName);
            end
        end  % function

        function result = getDeviceNameFromChannelName(self, channelName)
            % Get the DeviceName, given the channel name.
            iChannel = self.aoChannelIndexFromName(channelName) ;
            result = self.AnalogDeviceNames{iChannel} ;
        end  % function
        
%         function value = channelScaleFromName(self,channelName)
%             channelIndex = self.aoChannelIndexFromName(channelName) ;
%             if isnan(channelIndex) ,
%                 value = nan ;
%             else
%                 value = self.AnalogChannelScales(channelIndex) ;
%             end
%         end  % function

        function result=aoChannelIndexFromName(self,channelName)            
            iChannel=find(strcmp(channelName,self.AnalogChannelNames),1);
            if isempty(iChannel) ,
                result = nan ;
            else
                result = iChannel ;
            end
        end  % function

%         function result=channelUnitsFromName(self,channelName)
%             if isempty(channelName) ,
%                 result = '' ;
%             else
%                 iChannel=self.aoChannelIndexFromName(channelName);
%                 if isnan(iChannel) ,
%                     result='';
%                 else
%                     result=self.AnalogChannelUnits{iChannel};
%                 end
%             end
%         end  % function
        
        function result = getAnalogChannelUnits_(self)
            result = self.AnalogChannelUnits_ ;
        end  % function

        function result = getAnalogChannelScales_(self)
            result = self.AnalogChannelScales_ ;
        end

%         function set.AnalogChannelUnits(self,newValue)
%             newValue = cellfun(@strtrim,newValue,'UniformOutput',false);
%             oldValue=self.AnalogChannelUnits_;
%             isChangeable= ~(self.getNumberOfElectrodesClaimingAnalogChannel()==1);
%             editedNewValue=ws.fif(isChangeable,newValue,oldValue);
%             self.AnalogChannelUnits_=editedNewValue;
%             self.Parent.didSetAnalogChannelUnitsOrScales();            
%             %self.broadcast('DidSetAnalogChannelUnitsOrScales');
%         end  % function
        
        function setAnalogChannelScales_(self, newValue)
            self.AnalogChannelScales_ = newValue ;
        end  % function
        
%         function setAnalogChannelUnitsAndScales(self,newUnits,newScales)
%             newUnits = cellfun(@strtrim,newUnits,'UniformOutput',false);
%             isChangeable= ~(self.getNumberOfElectrodesClaimingAnalogChannel()==1);
%             oldUnits=self.AnalogChannelUnits_;
%             editedNewUnits=ws.fif(isChangeable,newUnits,oldUnits);
%             oldScales=self.AnalogChannelScales_;
%             editedNewScales=fif(isChangeable,newScales,oldScales);
%             self.AnalogChannelUnits_=editedNewUnits;
%             self.AnalogChannelScales_=editedNewScales;
%             self.Parent.didSetAnalogChannelUnitsOrScales();            
%             %self.broadcast('DidSetAnalogChannelUnitsOrScales');
%         end  % function
        
        function setSingleAnalogChannelUnits_(self, i, newValue)
            self.AnalogChannelUnits_{i} = newValue ;
        end  % function
        
        function setSingleAnalogChannelScale_(self, i, newValue)
            self.AnalogChannelScales_(i) = newValue ;
        end  % function
        
        function newChannelIndex = addAnalogChannel(self, deviceNameForNewChannel, newTerminalID)
            newChannelName = sprintf('AO%d',newTerminalID) ;
            
            self.AnalogDeviceNames_ = [self.AnalogDeviceNames_ {deviceNameForNewChannel} ] ;
            self.AnalogTerminalIDs_ = [self.AnalogTerminalIDs_ newTerminalID] ;
            self.AnalogChannelNames_ = [self.AnalogChannelNames_ {newChannelName}] ;
            self.AnalogChannelScales_ = [ self.AnalogChannelScales_ 1 ] ;
            self.AnalogChannelUnits_ = [ self.AnalogChannelUnits_ {'V'} ] ;
            self.IsAnalogChannelMarkedForDeletion_ = [  self.IsAnalogChannelMarkedForDeletion_ false ];
            
            newChannelIndex = length(self.AnalogChannelNames_) ;
        end  % function

        function wasDeleted = deleteMarkedAnalogChannels(self)
            % This has to be public so that the parent can call it, but it
            % should not be called by anyone but the parent.
            isToBeDeleted = self.IsAnalogChannelMarkedForDeletion_ ;
            %channelNamesToDelete = self.AnalogChannelNames_(isToBeDeleted) ;
            if all(isToBeDeleted)
                % Want everything to still be a row vector
                self.AnalogDeviceNames_ = cell(1,0) ;
                self.AnalogTerminalIDs_ = zeros(1,0) ;
                self.AnalogChannelNames_ = cell(1,0) ;
                self.AnalogChannelScales_ = zeros(1,0) ;
                self.AnalogChannelUnits_ = cell(1,0) ;
                self.IsAnalogChannelMarkedForDeletion_ = false(1,0) ;
            else
                isKeeper = ~isToBeDeleted ;
                self.AnalogDeviceNames_ = self.AnalogDeviceNames_(isKeeper) ;
                self.AnalogTerminalIDs_ = self.AnalogTerminalIDs_(isKeeper) ;
                self.AnalogChannelNames_ = self.AnalogChannelNames_(isKeeper) ;
                self.AnalogChannelScales_ = self.AnalogChannelScales_(isKeeper) ;
                self.AnalogChannelUnits_ = self.AnalogChannelUnits_(isKeeper) ;
                self.IsAnalogChannelMarkedForDeletion_ = self.IsAnalogChannelMarkedForDeletion_(isKeeper) ;
            end
            %self.syncIsAnalogChannelTerminalOvercommitted_() ;

            %self.Parent.didDeleteAnalogOutputChannels(channelNamesToDelete) ;
            %self.notifyLibraryThatDidChangeNumberOfOutputChannels_() ;
            
            wasDeleted = isToBeDeleted() ;
        end  % function
        
        function settingPrimaryDeviceName(self, newPrimaryDeviceName)            
            % All DI channels must use the primary device
            self.DigitalDeviceNames_(:) = {newPrimaryDeviceName} ;            
        end

        function setSingleAnalogChannelName(self, i, newValue)
            oldValue = self.AnalogChannelNames_{i} ;
            self.AnalogChannelNames_{i} = newValue ;
            self.StimulusLibrary_.renameChannel(oldValue, newValue) ;
        end
        
        function setSingleDigitalChannelName(self, i, newValue)
            oldValue = self.DigitalChannelNames_{i} ;
            self.DigitalChannelNames_{i} = newValue ;
            self.StimulusLibrary_.renameChannel(oldValue, newValue) ;
            %self.Parent.didSetDigitalOutputChannelName(didSucceed,oldValue,newValue);
        end
        
        function setSingleAnalogTerminalID(self, i, newValue)
            self.AnalogTerminalIDs_(i) = newValue ;
        end
    end  % methods block
    
    methods
        function mimic(self, other)
            % Cause self to resemble other.
            
            % Get the list of property names for this file type
            propertyNames = ws.listPropertiesForPersistence(self);
            
            % Set each property to the corresponding one
            for i = 1:length(propertyNames) ,
                thisPropertyName=propertyNames{i};
                if any(strcmp(thisPropertyName,{'StimulusLibrary_'})) ,                    
                    source = other.(thisPropertyName) ;  % source as in source vs target, not as in source vs destination                    
                    target = self.(thisPropertyName) ;
                    if isempty(target) ,
                        self.setPropertyValue_(thisPropertyName, ws.copy(source)) ;
                    else
                        target.mimic(source);
                    end
                else
                    if isprop(other,thisPropertyName) ,
                        source = other.getPropertyValue_(thisPropertyName) ;
                        self.setPropertyValue_(thisPropertyName, source) ;
                    end
                end
            end
        end  % function
    end  % public methods block

    methods
%         function set.IsDigitalChannelTimed(self,newValue)
%             self.setIsDigitalChannelTimed_(newValue) ;
%         end  % function
        
%         function set.DigitalOutputStateIfUntimed(self,newValue)
%             self.setDigitalOutputStateIfUntimed_(newValue) ;  % want to be able to override setter
%         end  % function

        function synchronizeTransientStateToPersistedStateHelper(self)
            nAnalogChannels = length(self.AnalogChannelNames_) ;
            self.IsAnalogChannelMarkedForDeletion_ = false(1, nAnalogChannels) ;
            nDigitalChannels = length(self.DigitalChannelNames_) ;
            self.IsDigitalChannelMarkedForDeletion_ = false(1, nDigitalChannels) ;
        end
    end

    methods 
        function sanitizePersistedState_(self)
            % This method should perform any sanity-checking that might be
            % advisable after loading the persistent state from disk.
            % This is often useful to provide backwards compatibility.
            
            % the length of AnalogChannelNames_ is the "true" number of AI
            % channels
            nAOChannels = length(self.AnalogChannelNames_) ;
            self.AnalogDeviceNames_ = ws.sanitizeRowVectorLength(self.AnalogDeviceNames_, nAOChannels, {'Dev1'}) ;
            self.AnalogTerminalIDs_ = ws.sanitizeRowVectorLength(self.AnalogTerminalIDs_, nAOChannels, 0) ;
            self.AnalogChannelScales_ = ws.sanitizeRowVectorLength(self.AnalogChannelScales_, nAOChannels, 1) ;
            self.AnalogChannelUnits_ = ws.sanitizeRowVectorLength(self.AnalogChannelUnits_, nAOChannels, {'V'}) ;
            %self.IsAnalogChannelActive_ = ws.sanitizeRowVectorLength(self.IsAnalogChannelActive_, nAOChannels, true) ;
            self.IsAnalogChannelMarkedForDeletion_ = ws.sanitizeRowVectorLength(self.IsAnalogChannelMarkedForDeletion_, nAOChannels, false) ;
            
            nDOChannels = length(self.DigitalChannelNames_) ;
            self.DigitalDeviceNames_ = ws.sanitizeRowVectorLength(self.DigitalDeviceNames_, nDOChannels, {'Dev1'}) ;
            self.DigitalTerminalIDs_ = ws.sanitizeRowVectorLength(self.DigitalTerminalIDs_, nDOChannels, 0) ;
            self.IsDigitalChannelTimed_ = ws.sanitizeRowVectorLength(self.IsDigitalChannelTimed_, nDOChannels, true) ;
            self.DigitalOutputStateIfUntimed_ = ws.sanitizeRowVectorLength(self.DigitalOutputStateIfUntimed_, nDOChannels, false) ;            
            self.IsDigitalChannelMarkedForDeletion_ = ws.sanitizeRowVectorLength(self.IsDigitalChannelMarkedForDeletion_, nDOChannels, false) ;
        end
    end  % protected methods block
    
    methods  % expose stim library methods at stimulation subsystem level
        function clearStimulusLibrary(self)
            self.StimulusLibrary_.clear() ;
        end  % function
        
        function setSelectedStimulusLibraryItemByClassNameAndIndex(self, className, index)
            self.StimulusLibrary_.setSelectedItemByClassNameAndIndex(className, index) ;
        end  % function
        
        function index = addNewStimulusSequence(self)
            index = self.StimulusLibrary_.addNewSequence() ;
        end  % function
        
        function duplicateSelectedStimulusLibraryItem(self)
            self.StimulusLibrary_.duplicateSelectedItem() ;
        end  % function
        
        function bindingIndex = addBindingToSelectedStimulusLibraryItem(self)
            bindingIndex = self.StimulusLibrary_.addBindingToSelectedItem() ;
        end  % function
        
        function bindingIndex = addBindingToStimulusLibraryItem(self, className, itemIndex)
            bindingIndex = self.StimulusLibrary_.addBindingToItem(className, itemIndex) ;
        end  % function
        
        function deleteMarkedBindingsFromSequence(self)
            self.StimulusLibrary_.deleteMarkedBindingsFromSequence() ;
        end  % function
        
        function mapIndex = addNewStimulusMap(self)
            mapIndex = self.StimulusLibrary_.addNewMap() ;
        end  % function
        
%         function addChannelToSelectedStimulusLibraryItem(self)
%             self.StimulusLibrary_.addChannelToSelectedItem() ;
%         end  % function
        
        function deleteMarkedChannelsFromSelectedStimulusLibraryItem(self)
            self.StimulusLibrary_.deleteMarkedChannelsFromSelectedItem() ;
        end  % function
        
        function stimulusIndex = addNewStimulus(self)
            stimulusIndex = self.StimulusLibrary_.addNewStimulus() ;
        end  % function        
        
        function result = isSelectedStimulusLibraryItemInUse(self)
            result = self.StimulusLibrary_.isSelectedItemInUse() ;
        end  % function        
        
        function deleteSelectedStimulusLibraryItem(self)
            self.StimulusLibrary_.deleteSelectedItem() ;
        end  % function        
        
        function result = selectedStimulusLibraryItemClassName(self)
            result = self.StimulusLibrary_.SelectedItemClassName ;
        end  % function        
        
        function result = selectedStimulusLibraryItemIndexWithinClass(self)
            result = self.StimulusLibrary_.SelectedItemIndexWithinClass ;
        end  % function        
        
        function didSetOutputableName = setSelectedStimulusLibraryItemProperty(self, propertyName, newValue)
            didSetOutputableName = self.StimulusLibrary_.setSelectedItemProperty(propertyName, newValue) ;
        end  % function        
        
        function setSelectedStimulusAdditionalParameter(self, iParameter, newString)
            self.StimulusLibrary_.setSelectedStimulusAdditionalParameter(iParameter, newString) ;
        end  % function        

        function setBindingOfSelectedSequenceToNamedMap(self, indexOfElementWithinSequence, newMapName)
            self.StimulusLibrary_.setBindingOfSelectedSequenceToNamedMap(indexOfElementWithinSequence, newMapName) ;
        end  % function                   
        
%         function setIsMarkedForDeletionForElementOfSelectedSequence(self, indexOfElementWithinSequence, newValue)
%             self.StimulusLibrary_.setIsMarkedForDeletionForElementOfSelectedSequence(indexOfElementWithinSequence, newValue) ;
%         end  % function                
        
        function setBindingOfSelectedMapToNamedStimulus(self, bindingIndex, newTargetName)
            self.StimulusLibrary_.setBindingOfSelectedMapToNamedStimulus(bindingIndex, newTargetName) ;
        end  % function                   
        
%         function setPropertyForElementOfSelectedMap(self, indexOfElementWithinMap, propertyName, newValue)
%             self.StimulusLibrary_.setPropertyForElementOfSelectedMap(indexOfElementWithinMap, propertyName, newValue) ;            
%         end  % function                
        
        function setSelectedStimulusLibraryItemWithinClassBindingProperty(self, className, bindingIndex, propertyName, newValue)
            self.StimulusLibrary_.setSelectedItemWithinClassBindingProperty(className, bindingIndex, propertyName, newValue) ;
        end  % method        
        
%         function plotSelectedStimulusLibraryItem(self, figureGH, samplingRate, channelNames, isChannelAnalog)
%             self.StimulusLibrary_.plotSelectedItemBang(figureGH, samplingRate, channelNames, isChannelAnalog) ;
%         end  % function            

        function [y, t] = previewStimulus(self, stimulusIndex, sampleRate)
            [y, t] = self.StimulusLibrary_.previewStimulus(stimulusIndex, sampleRate) ;
        end  % function            
        
        function [y, t] = previewStimulusMap(self, mapIndex, sampleRate, channelNames, isChannelAnalog)
            [y, t] = self.StimulusLibrary_.previewMap(mapIndex, sampleRate, channelNames, isChannelAnalog) ;
        end  % function            
        
        function result = selectedStimulusLibraryItemProperty(self, propertyName)
            result = self.StimulusLibrary_.selectedItemProperty(propertyName) ;
        end  % method        
        
        function result = indexOfStimulusLibraryClassSelection(self, className)
            result = self.StimulusLibrary_.indexOfClassSelection(className) ;
        end  % method                    
        
        function result = propertyFromEachStimulusLibraryItemInClass(self, className, propertyName) 
            % Result is a cell array, even it seems like it could/should be another kind of array
            result = self.StimulusLibrary_.propertyFromEachItemInClass(className, propertyName) ;
        end  % function
        
        function result = stimulusLibraryClassSelectionProperty(self, className, propertyName)
            result = self.StimulusLibrary_.classSelectionProperty(className, propertyName) ;
        end  % method        
        
%         function result = propertyForElementOfSelectedStimulusLibraryItem(self, indexOfElementWithinItem, propertyName)
%             result = self.StimulusLibrary_.propertyForElementOfSelectedItem(indexOfElementWithinItem, propertyName) ;
%         end  % function        
        
        function result = stimulusLibrarySelectedItemProperty(self, propertyName)
            result = self.StimulusLibrary_.selectedItemProperty(propertyName) ;
        end

        function result = stimulusLibrarySelectedItemBindingProperty(self, bindingIndex, propertyName)
            result = self.StimulusLibrary_.selectedItemBindingProperty(bindingIndex, propertyName) ;
        end
        
        function result = stimulusLibrarySelectedItemBindingTargetProperty(self, bindingIndex, propertyName)
            result = self.StimulusLibrary_.selectedItemBindingTargetProperty(bindingIndex, propertyName) ;
        end
        
        function result = stimulusLibraryItemProperty(self, className, index, propertyName)
            result = self.StimulusLibrary_.itemProperty(className, index, propertyName) ;
        end  % function                        

        function didSetOutputableName = setStimulusLibraryItemProperty(self, className, index, propertyName, newValue)
            didSetOutputableName = self.StimulusLibrary_.setItemProperty(className, index, propertyName, newValue) ;
        end  % function                        
        
        function result = stimulusLibraryItemBindingProperty(self, className, itemIndex, bindingIndex, propertyName)
            result = self.StimulusLibrary_.itemBindingProperty(className, itemIndex, bindingIndex, propertyName) ;
        end  % function                        
        
        function setStimulusLibraryItemBindingProperty(self, className, itemIndex, bindingIndex, propertyName, newValue)
            self.StimulusLibrary_.setItemBindingProperty(className, itemIndex, bindingIndex, propertyName, newValue) ;
        end  % function                        
        
        function result = stimulusLibraryItemBindingTargetProperty(self, className, itemIndex, bindingIndex, propertyName)
            result = self.StimulusLibrary_.itemBindingTargetProperty(className, itemIndex, bindingIndex, propertyName) ;
        end  % function             
        
        function result = isStimulusLibraryItemBindingTargetEmpty(self, className, itemIndex, bindingIndex)
            result = self.StimulusLibrary_.isItemBindingTargetEmpty(className, itemIndex, bindingIndex) ;
        end  % function
        
        function result = isStimulusLibrarySelectedItemBindingTargetEmpty(self, bindingIndex)
            result = self.StimulusLibrary_.isSelectedItemBindingTargetEmpty(bindingIndex) ;
        end  % function
        
        function result = isStimulusLibraryEmpty(self)
            result = self.StimulusLibrary_.isEmpty() ;
        end  % function                
        
        function result = isAStimulusLibraryItemSelected(self)
            result = self.StimulusLibrary_.isAnItemSelected() ;
        end  % function                
        
        function result = isAnyBindingMarkedForDeletionForStimulusLibrarySelectedItem(self)
            result = self.StimulusLibrary_.isAnyBindingMarkedForDeletionForSelectedItem() ;            
        end  % function        
        
        function setStimulusLibraryToSimpleLibraryWithUnitPulse(self, outputChannelNames)
            self.StimulusLibrary_.setToSimpleLibraryWithUnitPulse(outputChannelNames) ;            
        end
        
        function setSelectedOutputableByIndex(self, index)            
            self.StimulusLibrary_.setSelectedOutputableByIndex(index) ;
        end  % method
        
        function setSelectedOutputableByClassNameAndIndex(self, className, indexWithinClass)            
            self.StimulusLibrary_.setSelectedOutputableByClassNameAndIndex(className, indexWithinClass) ;
        end  % method
        
        function overrideStimulusLibraryMapDuration(self, sweepDuration)
            self.StimulusLibrary_.overrideMapDuration(sweepDuration) ;
        end  % function
        
        function releaseStimulusLibraryMapDuration(self)
            self.StimulusLibrary_.releaseMapDuration() ;
        end  % function
        
        function result = stimulusLibraryOutputableNames(self)
            result = self.StimulusLibrary_.outputableNames() ;            
        end  % function
        
        function result = stimulusLibrarySelectedOutputableProperty(self, propertyName)
            result = self.StimulusLibrary_.selectedOutputableProperty(propertyName) ;            
        end  % function
        
        function result = areStimulusLibraryMapDurationsOverridden(self)
            result = self.StimulusLibrary_.areMapDurationsOverridden() ;            
        end  % function
        
        function result = isStimulusLibraryItemInUse(self, className, itemIndex)
            result = self.StimulusLibrary_.isItemInUse(className, itemIndex) ;
        end  % function        
        
        function deleteStimulusLibraryItem(self, className, itemIndex)
            self.StimulusLibrary_.deleteItem(className, itemIndex) ;
        end  % function        
        
        function populateStimulusLibraryForTesting(self)
            self.StimulusLibrary_.populateForTesting() ;
        end  % function
        
%         function mimicStimulusLibrary_(self, newValue) 
%             self.StimulusLibrary_.mimic(newValue) ;
%         end
        
        function result = isStimulusLibrarySelfConsistent(self)
            result = self.StimulusLibrary_.isSelfConsistent() ;
        end
        
        function result = get.AnalogDeviceNames(self)
            result = self.AnalogDeviceNames_ ;
        end  % function
        
        function result = get.DigitalDeviceNames(self)
            result = self.DigitalDeviceNames_ ;
        end  % function        
        
        function result = getCurrentStimulusMapIndex(self, episodeIndexWithinSweep)
            doRepeatSequence = self.DoRepeatSequence_ ;
            result = self.StimulusLibrary_.getCurrentStimulusMapIndex(episodeIndexWithinSweep, doRepeatSequence) ;
        end        
        
        function data = ...
                calculateSignalsForMap(self, mapIndex, sampleRate, channelNames, isChannelAnalog, sweepIndexWithinSet)
            data = ...
                self.StimulusLibrary_.calculateSignalsForMap(mapIndex, sampleRate, channelNames, isChannelAnalog, sweepIndexWithinSet) ;
        end
        
    end  % public methods block    
    
    methods
        function setSingleAnalogDeviceName(self, i, newValue)
            % Checking done by parent
            self.AnalogDeviceNames_{i} = newValue ;
        end  % function        
        
        function setSingleDigitalDeviceName(self, i, newValue)
            % Checking done by parent
            self.DigitalDeviceNames_{i} = newValue ;
        end  % function                
    end  % public methods block
    
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
        function result = get.IsEnabled(self)
            result = self.IsEnabled_ ;
        end
        
        function set.IsEnabled(self, value)
            self.IsEnabled_ = value ;
        end
    end  % public methods block        
    
end  % classdef
