classdef ResolvedAnalogInputConflictTestCase < matlab.unittest.TestCase
    
    methods (TestMethodSetup)
        function setup(self) %#ok<MANU>
            %ws.reset() ;
        end
    end

    methods (TestMethodTeardown)
        function teardown(self) %#ok<MANU>
            %ws.reset() ;
        end
    end

    methods (Test)
        function theTest(self)
            wsModel=wavesurfer('--nogui', '--noprefs') ;

            wsModel.addAIChannel() ;
            wsModel.addAIChannel() ;
            wsModel.setSingleAIChannelTerminalID(2,0) ;  % this introduces a conflict
            try
                wsModel.playAndBlock() ;
            catch me
                if isequal(me.identifier, 'ws:ni:DAQmxError:n200489') ,
                    % we expect this error, so ignore it
                else
                    rethrow(me) ;
                end
            end
            wsModel.setSingleAIChannelTerminalID(2,1) ;  % this resolves the conflict
            wsModel.playAndBlock() ;  % this once errored, even though it shouldn't...
            self.verifyTrue(true) ;
        end  % function
    end  % test methods

end  % classdef
