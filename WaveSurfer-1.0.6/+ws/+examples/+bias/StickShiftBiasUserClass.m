classdef StickShiftBiasUserClass < ws.UserClass

    properties (Constant=true)
        serverIPAddressAsString = '127.0.0.1';                % bias listens on this ip
        serverFirstCameraPort = 5010 ;
        serverPortStride = 10 ;  % We assume bias listens at port serverFirstCameraPort for the 1st camera, serverFirstCameraPort+serverPortStride for 2nd, etc
    end  % constant properties
            
    properties (Dependent=true)
        cameraCount
    end
    
    properties (Access=protected, Transient=true)
        areCameraInterfacesInitialized_ = false
        cameraCount_
        biasCameraInterfaces_    % cell array of handles to bias camera interface object(s)
        isIInFrontend_ 
    end
    
    methods
        function self = StickShiftBiasUserClass()
            fprintf('Creating the BIAS user object\n') ;
        end

        function wake(self, rootModel)
            fprintf('Waking the BIAS user object\n') ;
            self.isIInFrontend_ = ( isa(rootModel,'ws.WavesurferModel') && rootModel.IsAwake ) ;
        end
        
        
        function delete(self)            
            fprintf('Deleting the BIAS user object\n') ;
            if self.areCameraInterfacesInitialized_ ,
                for i=1:self.cameraCount_ ,
                    try
                        %self.biasCameraInterfaces_{i}.disconnect();
                    catch me  %#ok<NASGU>
                        % ignore
                    end
                    delete(self.biasCameraInterfaces_{i});
                    self.biasCameraInterfaces_{i} = [] ;  % set the cell to empty, but don't change length of cell array
                end
                self.biasCameraInterfaces_ = cell(1,0) ;  % set to zero-length cell array
                self.areCameraInterfacesInitialized_ = false ;
            end
        end
        
        function willSaveToProtocolFile(self, wsModel)  %#ok<INUSD>
        end        
        
        function startingRun(self,~)
            fprintf('Starting a run.\n');
            if self.isIInFrontend_ ,
                if ~self.areCameraInterfacesInitialized_ ,
                    cameraCount = 0 ; %#ok<PROPLC>
                    for iCamera = 1:16 ,  % no one will have more cameras than this
                        portNumber = self.serverFirstCameraPort + (iCamera-1)*self.serverPortStride ;
                        cameraInterface = ws.examples.bias.SimpleBiasCameraInterface(self.serverIPAddressAsString, portNumber) ;
                        try
                            %cameraInterface.connectAndGetConfiguration() ;
                            response = cameraInterface.getStatus() ; %#ok<NASGU>
                        catch me
                            if isequal(me.identifier, 'SimpleBiasCameraInterface:unableToConnectToBIASServer') ,
                                break ;
                            else
                                me.rethrow() ;
                            end
                        end
                        % if get here, must have successfully connected to camera i
                        cameraCount = iCamera ;  %#ok<PROPLC>
                        self.biasCameraInterfaces_{1,iCamera} = cameraInterface ;
                    end
                    self.cameraCount_ = cameraCount ;  %#ok<PROPLC>
                    fprintf('Number of cameras found: %d\n', cameraCount) ;  %#ok<PROPLC>
                    % Finally, set the semaphore
                    self.areCameraInterfacesInitialized_ = true ;
                end
            end
        end
        
        function startingSweep(self ,~, ~)
            % Called just before each trial
            fprintf('Starting a sweep.\n');
            if self.isIInFrontend_ ,
                % Get status from each camera, to make sure BIAS is really
                % ready to go
                for i=1:self.cameraCount_ ,
                    response = self.biasCameraInterfaces_{i}.getStatus() ; 
                    if ~isfield(response,'value') || ~isfield(response.value,'connected') || ~response.value.connected ,
                        error('BIAS is not connected to camera %d', i-1) ;
                    end
                end
                for i=1:self.cameraCount_ ,
                    self.biasCameraInterfaces_{i}.startCapture() ;
                end

                % Wait for bias to be ready
                checkInterval = 0.1 ;  % s
                maxNumberOfChecks = 50 ;
                for i=1:maxNumberOfChecks ,
                    areAllCamerasCapturing = true ;
                    for j=1:self.cameraCount_ ,
                        response = self.biasCameraInterfaces_{j}.getStatus() ;   % call this just to make sure BIAS is capturing
                        if ~response.value.capturing ,
                            areAllCamerasCapturing = false ;
                            break ;
                        end
                    end
                    if areAllCamerasCapturing ,
                        break ;
                    else
                        pause(checkInterval) ;    % have to wait a bit for both cams to be capturing
                    end
                end                
            end
            
            % Warn if any cameras not capturing
            if ~areAllCamerasCapturing ,
                fprintf('Warning: Gave up waiting for all camera to be capturing\n') ;
            else
                fprintf('For sweep start, all cameras say they''re capturing\n') ;
            end 
        end
        
        function completingSweep(self, ~, ~)
            fprintf('Completing a sweep.\n');
            self.completingOrAbortingOrStoppingASweep_();
        end
        
        function abortingSweep(self,~)
            fprintf('Oh noes!  A sweep aborted.\n');
            self.completingOrAbortingOrStoppingASweep_();
        end
        
        function stoppingSweep(self,~)
            fprintf('A sweep was stopped.\n');
            self.completingOrAbortingOrStoppingASweep_();
        end
        
        function completingRun(self,~) %#ok<INUSD>
            % Called just after each set of trials (a.k.a. each
            % "experiment")
            fprintf('Completing a run.\n');
        end
        
        function abortingRun(self,~) %#ok<INUSD>
            % Called if a trial set goes wrong, after the call to
            % trialDidAbort()
            fprintf('Oh noes!  A run aborted.\n');
        end
        
        function stoppingRun(self,~) %#ok<INUSD>
            fprintf('A run was stopped.\n');
        end
        
        function dataAvailable(~,~)
        end

        function startingEpisode(~,~)
        end
        
        function completingEpisode(~,~)
        end
        
        function abortingEpisode(~,~)
        end
        
        function stoppingEpisode(~,~)
        end
        
        function result = get.cameraCount(self)
            result = self.cameraCount_ ;
        end
    end  % public methods


    methods (Access=private)
        function completingOrAbortingOrStoppingASweep_(self)
            % wait for bias to be done
            if self.isIInFrontend_ ,
                % Have to wait for a bit for all cameras to be done, we
                % discovered through a great deal of painful trial and
                % error.
                pause(0.1) ;
                
                % Tell bias to stop capturing
                for i=1:self.cameraCount_ ,
                    self.biasCameraInterfaces_{i}.stopCapture() ;
                end

                % Wait for it to stop capturing
                checkInterval = 0.1 ;  % s
                maxNumberOfChecks = 50 ;
                for i=1:maxNumberOfChecks ,
                    isACameraCapturing = false ;
                    for j=1:self.cameraCount_
                        response = self.biasCameraInterfaces_{j}.getStatus() ;   % call this just to make sure BIAS is done
                        if numel(response)==1 && isfield(response,'value') && isfield(response.value,'capturing') ,
                            isThisCameraCapturing = response.value.capturing ;
                            if numel(isThisCameraCapturing)==1 ,
                                if isThisCameraCapturing ,
                                    % A camera is still capturing, so we
                                    % can exit the for
                                    % j=1:self.cameraCount_ loop
                                    isACameraCapturing = true ;
                                    break ;
                                else
                                    % Communication is fine, and camera is
                                    % not capturing, so go on to check the
                                    % next camera.
                                end
                            else
                                fprintf('Problem communicating with camera %d at end of a sweep: We''ll assume it is done capturing.\n', j-1) ;                                    
                            end
                        else
                            fprintf('Problem communicating with camera %d at end of a sweep: We''ll assume it is done capturing.\n', j-1) ;                            
                        end
                    end
                    if isACameraCapturing ,
                        pause(checkInterval) ;    % have to wait a bit for all cams to be done
                    else
                        break ;
                    end
                end     
                
                if isACameraCapturing ,
                    fprintf('Warning: Gave up waiting for all cameras to not be capturing\n') ;
                else
                    fprintf('For sweep end, all cameras say they''re not capturing\n') ;
                end 
            end
        end
    end  % private methods block
    
    methods 
        % Allows access to protected and protected variables for encoding.
        function out = getPropertyValue_(self, name)
            out = self.(name);
        end
        
        % Allows access to protected and protected variables for encoding.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end        
    end  % protected methods block
    
    methods
        function mimic(self, other)
            ws.mimicBang(self, other) ;
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
        
end  % classdef

