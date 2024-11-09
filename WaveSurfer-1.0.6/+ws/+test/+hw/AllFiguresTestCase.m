classdef AllFiguresTestCase < matlab.unittest.TestCase
    % To run these tests, need to have an NI daq attached.  
    % (Can be a simulated daq board.)
    
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
            %thisDirName = fileparts(mfilename('fullpath')) ;
            [wsModel,wsController] = wavesurfer('--noprefs') ;  %#ok<ASGLU>

            % Launch some windows
            wsController.GeneralSettingsMenuItemActuated([],[]) ;  % Launch some windows
            wsController.ChannelsMenuItemActuated([],[]) ;            
            wsController.TriggersMenuItemActuated([],[]) ;
            wsController.StimulusLibraryMenuItemActuated([],[]) ;        
            wsController.StimulusPreviewMenuItemActuated([],[]) ;        
            wsController.UserCodeManagerMenuItemActuated([],[]) ;        
            wsController.ElectrodesMenuItemActuated([],[]) ;        
            wsController.TestPulseMenuItemActuated([],[]) ;        
            wsController.ManageFastProtocolsButtonActuated([],[]) ;
            
            generalSettingsFigure = wsController.GeneralSettingsController ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'ContinuousRadiobutton') ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'SweepBasedRadiobutton') ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'OverwriteCheckbox', true) ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'NextSweepEdit', '4') ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'SessionIndexCheckbox', true) ;
            ws.fakeControlActuationInTestBang(generalSettingsFigure, 'SessionIndexEdit', '7') ;
            
            wsController.quit() ;
            self.verifyTrue(true) ;
        end  % function
        
    end  % test methods
    
end  % classdef
