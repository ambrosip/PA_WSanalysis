classdef CounterTriggerHangTestCase < matlab.unittest.TestCase
    
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
            wsModel = wavesurfer('--nogui', '--noprefs') ;
            wsModel.NSweepsPerRun = 3 ;
            wsModel.IsStimulationEnabled = true ;
            wsModel.addCounterTrigger() ;
            wsModel.setTriggerProperty('counter', 1, 'RepeatCount', 2) ;
            wsModel.setTriggerProperty('counter', 1, 'Interval', 1.5) ;            
            wsModel.StimulationUsesAcquisitionTrigger = false ;
            wsModel.StimulationTriggerIndex = 2 ;  % this should be the newly-defined counter trigger
            didTimerCallbackFire = false ;            
            
            function timerCallback(source, event)  %#ok<INUSD>
                % The timer callback.  If WS is working properly, the timer
                % will be stopped before this fires at all.
                fprintf('timerCallback() fired.\n') ;
                didTimerCallbackFire = true ;
                wsModel.stop() ;
            end            
            
            timerToStopWavesurfer = timer('ExecutionMode', 'fixedDelay', ...
                                          'TimerFcn',@timerCallback, ...
                                          'StartDelay',20, ...
                                          'Period', 20);  % do this repeatedly in case first is missed
            start(timerToStopWavesurfer) ;
            wsModel.playAndBlock() ;  % this hangs at present, but if things work properly, this should return after ~ 3 seconds
            stop(timerToStopWavesurfer) ;  % stop the timer
            delete(timerToStopWavesurfer) ;            
            self.verifyFalse(didTimerCallbackFire) ;
        end  % function
    end  % test methods

end  % classdef
