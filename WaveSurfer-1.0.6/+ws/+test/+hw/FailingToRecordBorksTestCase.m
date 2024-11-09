classdef FailingToRecordBorksTestCase < matlab.unittest.TestCase
    % To run these tests, need to have an NI daq attached, pointed to by
    % the MDF.  (Can be a simulated daq board.)  Also, the MDF must be on the current path, 
    % and be named Machine_Data_File_WS_Test.m.
    
    methods (TestMethodSetup)
        function setup(self) %#ok<MANU>
            ws.clearDuringTests
        end
    end

    methods (TestMethodTeardown)
        function teardown(self) %#ok<MANU>
            ws.clearDuringTests
        end
    end

    methods (Test)
        function theTest(self)
            wsModel=wavesurfer('--nogui', '--noprefs');

            wsModel.addAIChannel() ;
            wsModel.addAIChannel() ;
            wsModel.addAIChannel() ;
            wsModel.addAOChannel() ;
                           
            % Turn on stimulation (there's a single pulse output by
            % default)
            wsModel.IsStimulationEnabled=true;

            % Turn on logging
            %wsModel.IsLoggingEnabled=true;

            % set the data file name
            thisFileName=mfilename();
            [~,dataFileBaseName]=fileparts(thisFileName);
            wsModel.DataFileBaseName=dataFileBaseName;

            % Want to make sure there's a pre-existing file by that name
            nextRunAbsoluteFileName=wsModel.NextRunAbsoluteFileName;
            if ~exist(nextRunAbsoluteFileName,'file');
                fid=fopen(nextRunAbsoluteFileName,'w');
                fclose(fid);
            end

            % wait a spell
            pause(0.1);
            
            % start the acq, which should error
            try
                wsModel.recordAndBlock();
            catch me
                if isequal(me.identifier,'wavesurfer:logFileAlreadyExists') ,
                    % ignore error
                else
                    rethrow(me);
                end
            end

            % Now check the "OK to overwrite" box.
            wsModel.IsOKToOverwriteDataFile=true;
            
            % wait a spell
            pause(0.1);
            
            % start the acq, which should work this time
            wsModel.recordAndBlock();  % this blocks
            
%             % Wait for acq to complete
%             dtBetweenChecks=0.1;  % s
%             maxTimeToWait=2;  % s
%             nTimesToCheck=ceil(maxTimeToWait/dtBetweenChecks);
%             for i=1:nTimesToCheck ,
%                 pause(dtBetweenChecks);
%                 if wsModel.NSweepsCompletedInThisRun>=1 ,
%                     break
%                 end
%             end                   

            % Delete the data file
            dataDirNameAbsolute=wsModel.DataFileLocation;
            dataFilePatternAbsolute=fullfile(dataDirNameAbsolute,[dataFileBaseName '*']);
            delete(dataFilePatternAbsolute);
            
            self.verifyEqual(wsModel.NSweepsCompletedInThisRun,1);            
        end  % function
        
    end  % test methods

 end  % classdef
