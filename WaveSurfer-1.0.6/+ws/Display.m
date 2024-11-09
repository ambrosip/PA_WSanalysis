classdef Display < ws.Model
    % Display manages the display and update of one or more Scope objects.
    
    properties (Dependent = true)
        IsEnabled
    end
    
    properties (Access = protected)
        IsEnabled_ = true
    end
    
    properties (Dependent = true)
        IsGridOn
        AreColorsNormal        
        DoShowZoomButtons        
        DoColorTraces
        UpdateRate  % the rate at which the scopes are updated, in Hz
        XOffset  % the x coord at the left edge of the scope windows
        %XSpan  % the trace duration shown in the scope windows
        IsAnalogChannelDisplayed  % 1 x nAIChannels
        IsDigitalChannelDisplayed  % 1 x nDIChannels
        AreYLimitsLockedTightToDataForAnalogChannel  % 1 x nAIChannels
        YLimitsPerAnalogChannel  % 2 x nAIChannels, 1st row is the lower limit, 2nd is the upper limit
        NPlots
        %XData
        %YData
        PlotHeightFromAnalogChannelIndex  % 1 x nAIChannels
        PlotHeightFromDigitalChannelIndex  % 1 x nDIChannels
        PlotHeightFromChannelIndex  % 1 x nChannels
        RowIndexFromAnalogChannelIndex  % 1 x nAIChannels
        RowIndexFromDigitalChannelIndex  % 1 x nDIChannels
        ChannelIndexWithinTypeFromPlotIndex  % 1 x NPlots
        IsAnalogFromPlotIndex  % 1 x NPlots
        ChannelIndexFromPlotIndex  % 1 x NPlots       
        PlotIndexFromChannelIndex  % 1 x nChannels
        PlotHeightFromPlotIndex  % 1 x NPlots
        IsXSpanSlavedToAcquistionDuration
        XSpanInPixels       
        XData
        YData
    end

    properties (Access = protected)
        IsGridOn_ = true
        AreColorsNormal_ = true  % if false, colors are inverted, approximately
        DoShowZoomButtons_ = true % if false, don't show buttons in the figure
        DoColorTraces_ = true % if false, traces are black/white
        XSpan_ 
        UpdateRate_
        IsXSpanSlavedToAcquistionDuration_
          % if true, the x span for all the scopes is set to the acquisiton
          % sweep duration
        IsAnalogChannelDisplayed_  % 1 x nAIChannels
        IsDigitalChannelDisplayed_  % 1 x nDIChannels
        AreYLimitsLockedTightToDataForAnalogChannel_  % 1 x nAIChannels
        YLimitsPerAnalogChannel_  % 2 x nAIChannels, 1st row is the lower limit, 2nd is the upper limit
        PlotHeightFromAnalogChannelIndex_  % 1 x nAIChannels
        PlotHeightFromDigitalChannelIndex_  % 1 x nDIChannels
        RowIndexFromAnalogChannelIndex_  % 1 x nAIChannels
        RowIndexFromDigitalChannelIndex_  % 1 x nDIChannels
    end
    
    properties (Access = protected, Transient=true)
        XAutoScroll_   % if true, x limits of all scopes will change to accomodate the data as it is acquired
        XOffset_
        ClearOnNextData_
        CachedDisplayXSpan_
        %XData_
        %YData_  % analog and digital together, all as doubles, but only for the *active* channels
        ChannelIndexWithinTypeFromPlotIndex_  % 1 x NPlots
        IsAnalogFromPlotIndex_  % 1 x NPlots
        ChannelIndexFromPlotIndex_  % 1 x NPlots (the channel index is in the list of all analog, then all digital, channels)
        PlotIndexFromChannelIndex_ % 1 x nChannels (this has nan's for channels that are not displayed)
        
        % The (downsampled for display) data currently being shown.        
        XSpanInPixels_
        XData_
        YData_
    end
    
    methods
        function self = Display()
            %self@ws.Subsystem() ;
            self.XOffset_ = 0 ;  % s
            self.XSpan_ = 1 ;  % s
            self.UpdateRate_ = 10 ;  % Hz
            self.XAutoScroll_ = false ;
            self.IsXSpanSlavedToAcquistionDuration_ = true ;
            self.IsAnalogChannelDisplayed_ = true(1,0) ; % 1 x nAIChannels
            self.IsDigitalChannelDisplayed_  = true(1,0) ; % 1 x nDIChannels
            self.AreYLimitsLockedTightToDataForAnalogChannel_ = false(1,0) ; % 1 x nAIChannels
            self.YLimitsPerAnalogChannel_ = zeros(2,0) ; % 2 x nAIChannels, 1st row is the lower limit, 2nd is the upper limit            
            self.XSpanInPixels_ = 400 ;  % for when we're running headless, this is a reasonable fallback value
            self.PlotHeightFromAnalogChannelIndex_ = zeros(1,0) ;
            self.PlotHeightFromDigitalChannelIndex_ = zeros(1,0) ;
            self.RowIndexFromAnalogChannelIndex_ = zeros(1,0) ;  % 1 x nAIChannels
            self.RowIndexFromDigitalChannelIndex_ = zeros(1,0) ;  % 1 x nDIChannels     
            self.updateMappingsFromPlotIndices_() ;
        end
        
        function delete(self)  %#ok<INUSD>
        end
        
        function result = get.NPlots(self)
            result = length(self.IsAnalogChannelDisplayed_) + length(self.IsDigitalChannelDisplayed_) ;
        end
        
%         function result = get.XData(self)
%             result = self.XData_ ;
%         end
%         
%         function result = get.YData(self)
%             result = self.YData_ ;
%         end
        
        function result = get.AreYLimitsLockedTightToDataForAnalogChannel(self)
            result = self.AreYLimitsLockedTightToDataForAnalogChannel_ ;
        end
        
        function hereIsXSpanInPixels(self, xSpanInPixels)
            self.XSpanInPixels_ = xSpanInPixels ;
        end        
        
        function result = get.IsAnalogChannelDisplayed(self)
            result = self.IsAnalogChannelDisplayed_ ;
        end
        
        function toggleIsAnalogChannelDisplayed(self, aiChannelIndex, nAIChannels) 
            if isnumeric(aiChannelIndex) && isscalar(aiChannelIndex) && isreal(aiChannelIndex) && (aiChannelIndex==round(aiChannelIndex)) ,
                %nAIChannels = self.Parent.Acquisition.NAnalogChannels ;
                if 1<=aiChannelIndex && aiChannelIndex<=nAIChannels ,
                    currentValue = self.IsAnalogChannelDisplayed_(aiChannelIndex) ;
                    self.IsAnalogChannelDisplayed_(aiChannelIndex) = ~currentValue ;
                    self.updateMappingsFromPlotIndices_() ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            %self.broadcast('Update');
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to toggleIsAnalogChannelDisplayed must be a valid AI channel index') ;
            end                
        end
        
        function toggleIsDigitalChannelDisplayed(self, diChannelIndex, nDIChannels) 
            if isnumeric(diChannelIndex) && isscalar(diChannelIndex) && isreal(diChannelIndex) && (diChannelIndex==round(diChannelIndex))
                %nDIChannels = self.Parent.Acquisition.NDigitalChannels ;
                if 1<=diChannelIndex && diChannelIndex<=nDIChannels ,
                    currentValue = self.IsDigitalChannelDisplayed_(diChannelIndex) ;
                    self.IsDigitalChannelDisplayed_(diChannelIndex) = ~currentValue ;
                    self.updateMappingsFromPlotIndices_() ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            %self.broadcast('Update');
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to toggleIsDigitalChannelDisplayed must be a valid DI channel index') ;
            end                
        end
        
        function result = get.IsDigitalChannelDisplayed(self)
            result = self.IsDigitalChannelDisplayed_ ;
        end

        function result = get.PlotHeightFromAnalogChannelIndex(self)
            result = self.PlotHeightFromAnalogChannelIndex_ ;
        end
        
        function result = get.PlotHeightFromDigitalChannelIndex(self)
            result = self.PlotHeightFromDigitalChannelIndex_ ;
        end
        
        function result = get.PlotHeightFromChannelIndex(self)
            result = horzcat(self.PlotHeightFromAnalogChannelIndex_, self.PlotHeightFromDigitalChannelIndex_) ;
        end

        function result = get.PlotHeightFromPlotIndex(self)
            plotHeightFromChannelIndex = self.PlotHeightFromChannelIndex ;
            result = plotHeightFromChannelIndex(self.ChannelIndexFromPlotIndex) ;
        end
        
        function result = get.RowIndexFromAnalogChannelIndex(self)
            result = self.RowIndexFromAnalogChannelIndex_ ;
        end
        
        function result = get.RowIndexFromDigitalChannelIndex(self)
            result = self.RowIndexFromDigitalChannelIndex_ ;
        end
        
        function result = get.ChannelIndexWithinTypeFromPlotIndex(self)
            result = self.ChannelIndexWithinTypeFromPlotIndex_ ;
        end

        function result = get.ChannelIndexFromPlotIndex(self)
            result = self.ChannelIndexFromPlotIndex_ ;
        end

        function result = get.PlotIndexFromChannelIndex(self)
            result = self.PlotIndexFromChannelIndex_ ;
        end
        
        function result = get.IsAnalogFromPlotIndex(self)
            result = self.IsAnalogFromPlotIndex_ ;
        end
        
        function value = get.UpdateRate(self)
            value = self.UpdateRate_;
        end
        
        function set.UpdateRate(self, newValue)
            self.UpdateRate_ = newValue;
        end
        
        function value = getXSpan(self)
            value = self.XSpan_ ;
        end
        
        function setXSpan(self, newValue)            
            self.XSpan_ = newValue ;
        end  % function
                
        function value = get.XOffset(self)
            value = self.XOffset_;
        end
                
        function set.XOffset(self, newValue)
            self.XOffset_ = newValue ;
        end
        
        function value = get.YLimitsPerAnalogChannel(self)
            value = self.YLimitsPerAnalogChannel_ ;
        end

        function channelIndex = setYLimitsForSinglePlot(self, plotIndex, newValue)
            if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
                isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
                nPlots = length(isAnalogFromPlotIndex) ;
                if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
                    channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
                    self.YLimitsPerAnalogChannel_(:,channelIndex) = double(newValue(:)) ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'First argument to setYLimitsForSinglePlot() must be a valid AI plot index') ;
            end                            
        end
        
%         function setYLimitsForSingleAIChannel_(self, aiChannelIndex, newValue)
%             % This has an underscore b/c it doesn't do an update
%             if isnumeric(newValue) && isequal(size(newValue),[1 2]) && newValue(1)<=newValue(2) ,
%                 self.YLimitsPerAnalogChannel_(:,aiChannelIndex) = double(newValue') ;
%             end
%         end
        
        function value = get.IsXSpanSlavedToAcquistionDuration(self)
            value = self.IsXSpanSlavedToAcquistionDuration_;
        end  % function
        
        function set.IsXSpanSlavedToAcquistionDuration(self, newValue)
            self.IsXSpanSlavedToAcquistionDuration_ = newValue ;
        end
        
%         function value = get.IsXSpanSlavedToAcquistionDurationSettable(self)
%             value = self.Parent.AreSweepsFiniteDuration ;
%         end  % function       

%         function didRemoveElectrodes(self)
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%             self.broadcast('Update') ;
%         end

%         function didSetAnalogChannelUnitsOrScales(self)
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%             self.broadcast('Update') ;
%         end
        
        function startingRun(self, xSpan, sweepDuration)
            %self.XOffset = 0 ;
            %self.XSpan = self.XSpan;  % in case user has zoomed in on one or more scopes, want to reset now
            %self.XAutoScroll_ = (self.Parent.AreSweepsContinuous) ;
            self.XAutoScroll_ = (xSpan<sweepDuration) ;
        end  % function
        
        function completingRun(self)
            self.completingOrStoppingOrAbortingRun_();
        end
        
        function stoppingRun(self)
            self.completingOrStoppingOrAbortingRun_();
        end
        
        function abortingRun(self)
            self.completingOrStoppingOrAbortingRun_();
        end
        
        function didAddAnalogInputChannel(self)
            self.IsAnalogChannelDisplayed_ = horzcat(self.IsAnalogChannelDisplayed_, true) ;
            self.AreYLimitsLockedTightToDataForAnalogChannel_ = horzcat(self.AreYLimitsLockedTightToDataForAnalogChannel_, false) ;
            self.YLimitsPerAnalogChannel_ = horzcat(self.YLimitsPerAnalogChannel_, [-10 +10]') ;
            self.PlotHeightFromAnalogChannelIndex_ = horzcat(self.PlotHeightFromAnalogChannelIndex_, 1) ;
            nRowsBefore = length(self.RowIndexFromAnalogChannelIndex_) + length(self.RowIndexFromDigitalChannelIndex_) ;
            self.RowIndexFromAnalogChannelIndex_ = horzcat(self.RowIndexFromAnalogChannelIndex_, nRowsBefore+1) ;
            self.updateMappingsFromPlotIndices_() ;
        end
         
        function didAddDigitalInputChannel(self)
            self.IsDigitalChannelDisplayed_(1,end+1) = true ;
            self.PlotHeightFromDigitalChannelIndex_ = horzcat(self.PlotHeightFromDigitalChannelIndex_, 1) ;
            nRowsBefore = length(self.RowIndexFromAnalogChannelIndex_) + length(self.RowIndexFromDigitalChannelIndex_) ;
            self.RowIndexFromDigitalChannelIndex_ = horzcat(self.RowIndexFromDigitalChannelIndex_, nRowsBefore+1) ;
            self.updateMappingsFromPlotIndices_() ;
            %self.clearData_() ;
        end

        function didDeleteAnalogInputChannels(self, wasDeleted)
            wasKept = ~wasDeleted ;
            self.IsAnalogChannelDisplayed_ = self.IsAnalogChannelDisplayed_(wasKept) ;
            self.AreYLimitsLockedTightToDataForAnalogChannel_ = self.AreYLimitsLockedTightToDataForAnalogChannel_(wasKept) ;
            self.YLimitsPerAnalogChannel_ = self.YLimitsPerAnalogChannel_(:,wasKept) ;
            self.PlotHeightFromAnalogChannelIndex_ = self.PlotHeightFromAnalogChannelIndex_(wasKept) ;
            self.RowIndexFromAnalogChannelIndex_ = self.RowIndexFromAnalogChannelIndex_(wasKept) ;
            [self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_] = ...
                ws.Display.renormalizeRowIndices(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_) ;            
            self.updateMappingsFromPlotIndices_() ;
            %self.clearData_() ;
        end
        
        function didDeleteDigitalInputChannels(self, wasDeleted)            
            wasKept = ~wasDeleted ;
            self.IsDigitalChannelDisplayed_ = self.IsDigitalChannelDisplayed_(wasKept) ;
            self.PlotHeightFromDigitalChannelIndex_ = self.PlotHeightFromDigitalChannelIndex_(wasKept) ;
            self.RowIndexFromDigitalChannelIndex_ = self.RowIndexFromDigitalChannelIndex_(wasKept) ;
            [self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_] = ...
                ws.Display.renormalizeRowIndices(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_) ;            
            self.updateMappingsFromPlotIndices_() ;
            %self.clearData_() ;
        end
        
%         function didSetAnalogInputChannelName(self, didSucceed, oldValue, newValue) %#ok<INUSD>
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%             self.broadcast('Update') ;            
%         end
        
%         function didSetDigitalInputChannelName(self, didSucceed, oldValue, newValue) %#ok<INUSD>
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%             self.broadcast('Update') ;            
%         end
        
%         function didSetIsInputChannelActive(self) 
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%             self.broadcast('Update') ;            
%         end
        
%         function toggleIsGridOn_(self)
%             self.IsGridOn = ~(self.IsGridOn) ;
%         end

%         function toggleAreColorsNormal_(self)
%             self.AreColorsNormal = ~(self.AreColorsNormal) ;
%         end
        
        function set.IsGridOn(self,newValue)
            self.IsGridOn_ = newValue ;
        end
        
        function result = get.IsGridOn(self)
            result = self.IsGridOn_ ;
        end
            
        function set.AreColorsNormal(self, newValue)
            self.AreColorsNormal_ = newValue ;
        end
        
        function result = get.AreColorsNormal(self)
            result = self.AreColorsNormal_ ;
        end
            
        function set.DoShowZoomButtons(self, newValue)
            self.DoShowZoomButtons_ = logical(newValue) ;
        end
        
        function result = get.DoShowZoomButtons(self)
            result = self.DoShowZoomButtons_ ;
        end                    
        
        function set.DoColorTraces(self,newValue)
            self.DoColorTraces_ = newValue ;
        end
        
        function result = get.DoColorTraces(self)
            result = self.DoColorTraces_ ;
        end                    
        
        function channelIndex = scrollUp(self, plotIndex)  % works on analog channels only
            if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
                isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
                nPlots = length(isAnalogFromPlotIndex) ;
                if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
                    channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
                    yLimits = self.YLimitsPerAnalogChannel_(:,channelIndex) ;  % NB: a 2-el col vector
                    yMiddle=mean(yLimits);
                    ySpan=diff(yLimits);
                    yRadius=0.5*ySpan;
                    newYLimits=(yMiddle+0.1*ySpan)+yRadius*[-1 +1]' ;
                    self.YLimitsPerAnalogChannel_(:,channelIndex) = newYLimits ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            %self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to scrollUp() must be a valid AI channel index') ;
            end                            
        end  % function
        
        function channelIndex = scrollDown(self, plotIndex)  % works on analog channels only
            if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
                isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
                nPlots = length(isAnalogFromPlotIndex) ;
                if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
                    channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
                    yLimits = self.YLimitsPerAnalogChannel_(:,channelIndex) ;  % NB: a 2-el col vector
                    yMiddle=mean(yLimits);
                    ySpan=diff(yLimits);
                    yRadius=0.5*ySpan;
                    newYLimits=(yMiddle-0.1*ySpan)+yRadius*[-1 +1]' ;
                    self.YLimitsPerAnalogChannel_(:,channelIndex) = newYLimits ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            %self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to scrollDown() must be a valid AI channel index') ;
            end                           
        end  % function
                
        function channelIndex = zoomIn(self, plotIndex)  % works on analog channels only
            if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
                isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
                nPlots = length(isAnalogFromPlotIndex) ;
                if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
                    channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
                    yLimits = self.YLimitsPerAnalogChannel_(:,channelIndex) ;  % NB: a 2-el col vector
                    yMiddle=mean(yLimits);
                    yRadius=0.5*diff(yLimits);
                    newYLimits=yMiddle+0.5*yRadius*[-1 +1]' ;
                    self.YLimitsPerAnalogChannel_(:,channelIndex) = newYLimits ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to zoomOut() must be a valid AI plot index') ;
            end                                    
        end  % function
                
        function channelIndex = zoomOut(self, plotIndex)  % works on analog channels only
            if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
                isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
                nPlots = length(isAnalogFromPlotIndex) ;
                if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
                    channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
                    yLimits = self.YLimitsPerAnalogChannel_(:,channelIndex) ;  % NB: a 2-el col vector
                    yMiddle=mean(yLimits);
                    yRadius=0.5*diff(yLimits);
                    newYLimits=yMiddle+2*yRadius*[-1 +1]' ;
                    self.YLimitsPerAnalogChannel_(:,channelIndex) = newYLimits ;
                    isValid = true ;
                else
                    isValid = false ;
                end
            else
                isValid = false ;
            end
            %self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
            if ~isValid ,
                error('ws:invalidPropertyValue', ...
                      'Argument to zoomIn() must be a valid AI plot index') ;
            end                
        end  % function
                
%         function setYAxisLimitsTightToData(self, plotIndex)            
%             if isnumeric(plotIndex) && isscalar(plotIndex) && isreal(plotIndex) && (plotIndex==round(plotIndex)) && 1<=plotIndex,
%                 isAnalogFromPlotIndex = self.IsAnalogFromPlotIndex_ ;
%                 nPlots = length(isAnalogFromPlotIndex) ;
%                 if plotIndex <= nPlots && isAnalogFromPlotIndex(plotIndex),
%                     channelIndex = self.ChannelIndexWithinTypeFromPlotIndex_(plotIndex) ;
%                     self.setYAxisLimitsTightToData_(channelIndex) ;
%                     isValid = true ;
%                 else
%                     isValid = false ;
%                 end
%             else
%                 isValid = false ;
%             end
%             self.broadcast('DidSetYAxisLimits', plotIndex, channelIndex);
%             if ~isValid ,
%                 error('ws:invalidPropertyValue', ...
%                       'Argument to setYAxisLimitsTightToData() must be a valid AI channel index') ;
%             end                
%         end  % function

        function setAreYLimitsLockedTightToDataForSingleChannel_(self, channelIndex, newValue)            
            % this has an underscore b/c it doesn't do an update
            self.AreYLimitsLockedTightToDataForAnalogChannel_(channelIndex) = newValue ;
        end        
        
%         function didSetAnalogInputTerminalID(self)
%             % This should only be called by the parent, hence the
%             % underscore.
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%         end        
        
%         function didSetDigitalInputTerminalID_(self)
%             % This should only be called by the parent, hence the
%             % underscore.
%             %self.clearData_() ;
%             self.broadcast('ClearData') ;
%         end
        
        function setPlotHeightsAndOrder(self, isDisplayed, plotHeights, rowIndexFromChannelIndex)
            % Typically called by ws.PlotArrangementDialogController after OK
            % button is pressed.  Does no argument checking.
            nAIChannels = length(self.IsAnalogChannelDisplayed_) ;
            % Set properties
            
            % We'll need to decide whether to clear the displayed traces or
            % not.  We do this only if the height of one or more plots is
            % changing.  To determine whether this is the case, we need to
            % cache the original values of some things, and compare them to
            % the new values.
%             oldIsAnalogChannelDisplayed = self.IsAnalogChannelDisplayed_ ;
%             oldIsDigitalChannelDisplayed = self.IsDigitalChannelDisplayed_ ;
%             oldPlotHeightFromAnalogChannelIndex = self.PlotHeightFromAnalogChannelIndex_ ;
%             oldPlotHeightFromDigitalChannelIndex = self.PlotHeightFromDigitalChannelIndex_ ;
            newIsAnalogChannelDisplayed = isDisplayed(1:nAIChannels) ;
            newIsDigitalChannelDisplayed = isDisplayed(nAIChannels+1:end) ;
            newPlotHeightFromAnalogChannelIndex = plotHeights(1:nAIChannels) ;
            newPlotHeightFromDigitalChannelIndex = plotHeights(nAIChannels+1:end) ;
%             doNeedToClearDataCache = ...
%                 ~isequal(newIsAnalogChannelDisplayed, oldIsAnalogChannelDisplayed) || ...
%                 ~isequal(newIsDigitalChannelDisplayed, oldIsDigitalChannelDisplayed) || ...
%                 ~isequal(newPlotHeightFromAnalogChannelIndex, oldPlotHeightFromAnalogChannelIndex) || ...
%                 ~isequal(newPlotHeightFromDigitalChannelIndex, oldPlotHeightFromDigitalChannelIndex) ;
            
            % OK, now we can actually set instance variables                
            self.IsAnalogChannelDisplayed_ = newIsAnalogChannelDisplayed ;
            self.IsDigitalChannelDisplayed_ = newIsDigitalChannelDisplayed ;
            self.PlotHeightFromAnalogChannelIndex_ = newPlotHeightFromAnalogChannelIndex ;
            self.PlotHeightFromDigitalChannelIndex_ = newPlotHeightFromDigitalChannelIndex ;
            self.RowIndexFromAnalogChannelIndex_ = rowIndexFromChannelIndex(1:nAIChannels) ;
            self.RowIndexFromDigitalChannelIndex_ = rowIndexFromChannelIndex(nAIChannels+1:end) ;
            self.updateMappingsFromPlotIndices_() ;
%             if doNeedToClearData ,
%                 self.broadcast('ClearData') ;
%             end
%             self.broadcast('Update') ;            
        end
    end  % public methods block
    
    methods (Access=protected)        
%         function setYAxisLimitsTightToData_(self, aiChannelIndex)            
%             % this core function does no arg checking and doesn't call
%             % .broadcast.  It just mutates the state.
%             yMinAndMax=self.dataYMinAndMax_(aiChannelIndex);
%             if any(~isfinite(yMinAndMax)) ,
%                 return
%             end
%             yCenter=mean(yMinAndMax);
%             yRadius=0.5*diff(yMinAndMax);
%             if yRadius==0 ,
%                 yRadius=0.001;
%             end
%             newYLimits = yCenter + 1.05*yRadius*[-1 +1]' ;
%             self.YLimitsPerAnalogChannel_(:,aiChannelIndex) = newYLimits ;            
%         end
%         
%         function yMinAndMax=dataYMinAndMax_(self, aiChannelIndex)
%             % Min and max of the data, across all plotted channels.
%             % Returns a 1x2 array.
%             % If all channels are empty, returns [+inf -inf].
%             activeChannelIndexFromChannelIndex = self.Parent.Acquisition.ActiveChannelIndexFromChannelIndex ;
%             indexWithinData = activeChannelIndexFromChannelIndex(aiChannelIndex) ;
%             y = self.YData(:,indexWithinData) ;
%             yMinRaw=min(y);
%             yMin=ws.fif(isempty(yMinRaw),+inf,yMinRaw);
%             yMaxRaw=max(y);
%             yMax=ws.fif(isempty(yMaxRaw),-inf,yMaxRaw);            
%             yMinAndMax=double([yMin yMax]);
%         end
        
        function completingOrStoppingOrAbortingRun_(self)
            if ~isempty(self.CachedDisplayXSpan_)
                self.XSpan = self.CachedDisplayXSpan_;
            end
            self.CachedDisplayXSpan_ = [];
        end        
                
%         function clearData_(self)
%             self.XData_ = zeros(0,1) ;
%             acquisition = self.Parent.Acquisition ;
%             nActiveChannels = acquisition.NActiveAnalogChannels + acquisition.NActiveDigitalChannels ;
%             self.YData_ = zeros(0,nActiveChannels) ;
%         end
%         
%         function indicesOfAIChannelsNeedingYLimitUpdate = addData_(self, t, recentScaledAnalogData, recentRawDigitalData, sampleRate)
%             % t is a scalar, the time stamp of the scan *just after* the
%             % most recent scan.  (I.e. it is one dt==1/fs into the future.
%             % Queue Doctor Who music.)
% 
%             % Get the uint8/uint16/uint32 data out of recentRawDigitalData
%             % into a matrix of logical data, then convert it to doubles and
%             % concat it with the recentScaledAnalogData, storing the result
%             % in yRecent.
%             nActiveDigitalChannels = self.Parent.Acquisition.NActiveDigitalChannels ;
%             if nActiveDigitalChannels==0 ,
%                 yRecent = recentScaledAnalogData ;
%             else
%                 % Might need to write a mex function to quickly translate
%                 % recentRawDigitalData to recentDigitalData.
%                 nScans = size(recentRawDigitalData,1) ;                
%                 recentDigitalData = zeros(nScans,nActiveDigitalChannels) ;
%                 for j = 1:nActiveDigitalChannels ,
%                     recentDigitalData(:,j) = bitget(recentRawDigitalData,j) ;
%                 end
%                 % End of code that might need to mex-ify
%                 yRecent = horzcat(recentScaledAnalogData, recentDigitalData) ;
%             end
%             
%             % Compute a timeline for the new data            
%             nNewScans = size(yRecent, 1) ;
%             dt = 1/sampleRate ;  % s
%             t0 = t - dt*nNewScans ;  % timestamp of first scan in newData
%             xRecent = t0 + dt*(0:(nNewScans-1))' ;
%             
%             % Figure out the downsampling ratio
%             self.broadcast('ItWouldBeNiceToKnowXSpanInPixels') ;
%               % At this point, self.XSpanPixels_ should be set to the
%               % correct value, or the fallback value if there's no view
%             %xSpanInPixels=ws.ScopeFigure.getWidthInPixels(self.AxesGH_);
%             xSpanInPixels = self.XSpanInPixels_ ;
%             xSpan = self.XSpan ;
%             r = ws.ratioSubsampling(dt, xSpan, xSpanInPixels) ;
%             
%             % Downsample the new data
%             [xForPlottingNew, yForPlottingNew] = ws.minMaxDownsampleMex(xRecent, yRecent, r) ;            
%             
%             % deal with XData
%             xAllOriginal = self.XData ;  % these are already downsampled
%             yAllOriginal = self.YData ;            
%             
%             % Concatenate the old data that we're keeping with the new data
%             xAllProto = vertcat(xAllOriginal, xForPlottingNew) ;
%             yAllProto = vertcat(yAllOriginal, yForPlottingNew) ;
%             
%             % Trim off scans that would be off the screen anyway
%             doKeepScan = (self.XOffset_<=xAllProto) ;
%             xNew = xAllProto(doKeepScan) ;
%             yNew = yAllProto(doKeepScan,:) ;
% 
%             % Commit the data to self
%             self.XData_ = xNew ;
%             self.YData_ = yNew ;
%             
% %             % Update the x offset in the scope to match that in the Display
% %             % subsystem
% %             fprintf('xOffset: %20g     self.XOffset: %20g\n', xOffset, self.XOffset) ;
% %             if xOffset ~= self.XOffset , 
% %                 fprintf('About to change x offset\n') ;
% %                 self.XOffset = xOffset ;
% %             end
%             
%             % Change the y limits to match the data, if appropriate
%             indicesOfAIChannelsNeedingYLimitUpdate = self.setYAxisLimitsTightToDataIfAreYLimitsLockedTightToData_() ;
%         end        
    end
        
    methods    
        function startingSweep(self)
            self.ClearOnNextData_ = true;
        end
        
        function [doesNeedClear, doesNeedDidSetXOffset] = ...
                dataAvailable(self, isSweepBased, t, scaledAnalogData, rawAnalogData, rawDigitalData, timeSinceRunStartAtStartOfData, xSpan)  %#ok<INUSL>
            % t is a scalar, the time stamp of the scan *just after* the
            % most recent scan.  (I.e. it is one dt==1/fs into the future.
            % Queue Doctor Who music.)
            
            % Get relevant state
            doesNeedClear = self.ClearOnNextData_ ;
            originalXOffset = self.XOffset_ ;
            xAutoScroll = self.XAutoScroll_ ;

            % update the x offset
            if xAutoScroll ,                
                if doesNeedClear ,
                    doesNeedDidSetXOffset = true ;
                    newXOffset = 0 ;
                else                    
                    scale = min(1,xSpan) ;
                    tNudged = scale * ceil(100*t/scale)/100 ;  % Helps keep the axes aligned to tidy numbers
                    xOffsetNudged = tNudged - xSpan ;
                    if xOffsetNudged > originalXOffset ,
                        doesNeedDidSetXOffset = true ;
                        newXOffset = xOffsetNudged ;
                    else
                        doesNeedDidSetXOffset = false ;
                    end
                end
            else
                doesNeedDidSetXOffset = false ;
            end
            
            % Set state
            if doesNeedDidSetXOffset ,
                self.XOffset_ = newXOffset ;
            end
            self.ClearOnNextData_ = false ;
        end  % function
        
%         function result = getPlotIndexFromChannelIndex(self)
%             % The "channel index" here is is equal to the AI channel index
%             % for AI channels, and is equal to the DI channel index
%             % plus the number of AI channels for DI channels.  The plot
%             % index is set to nan for undisplayed channels.
%             isChannelDisplayed = horzcat(self.IsAnalogChannelDisplayed_, self.IsDigitalChannelDisplayed_) ;
%             rowIndexFromChannelIndex = horzcat(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_) ;
%             rowIndexFromChannelIndexAmongDisplayed = rowIndexFromChannelIndex(isChannelDisplayed) ;
%             plotIndexFromChannelIndexAmongDisplayed = ws.sortedOrder(rowIndexFromChannelIndexAmongDisplayed) ;
%             nChannels = length(isChannelDisplayed) ;
%             result = nan(1,nChannels) ;
%             result(isChannelDisplayed) = plotIndexFromChannelIndexAmongDisplayed ;            
%         end        
        
%         function [channelIndexWithinTypeFromPlotIndex, isAnalogFromPlotIndex] = getChannelIndexFromPlotIndexMapping(self)
%             isAnalogChannelDisplayed = self.IsAnalogChannelDisplayed_ ;
%             isDigitalChannelDisplayed = self.IsDigitalChannelDisplayed_ ;
%             nAnalogChannels = length(isAnalogChannelDisplayed) ;
%             nDigitalChannels = length(isDigitalChannelDisplayed) ;
%             %nChannels = nAnalogChannels + nDigitalChannels ;
%             isAnalogFromChannelIndex = horzcat( true(1,nAnalogChannels), false(1,nDigitalChannels) ) ;
%             isDisplayedFromChannelIndex = horzcat(isAnalogChannelDisplayed, isDigitalChannelDisplayed) ;
%             rowIndexFromChannelIndex = horzcat(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_) ;
%             channelIndexFromPlotIndex = ws.Display.computeChannelIndexFromPlotIndexMapping(rowIndexFromChannelIndex, isDisplayedFromChannelIndex) ;
%             isAnalogFromPlotIndex = isAnalogFromChannelIndex(channelIndexFromPlotIndex) ;
%             channelIndexWithinTypeFromPlotIndex = ...
%                 arrayfun(@(channelIndex)(ws.fif(channelIndex>nAnalogChannels,channelIndex-nAnalogChannels,el)), channelIndexFromPlotIndex) ;
%         end
        
%         function didSetAreSweepsFiniteDuration(self)
%             % Called by the parent to notify of a change to the acquisition
%             % duration            
%             self.broadcast('ClearData') ;
%             self.broadcast('DidSetXSpan') ;
%         end
        
%         function didSetSweepDurationIfFinite(self, isXSpanSlavedToAcquistionDuration)
%             % Called by the parent to notify of a change to the acquisition
%             % duration.            
%             if isXSpanSlavedToAcquistionDuration ,
%                 self.broadcast('ClearData') ;
%             end                
%             self.broadcast('DidSetXSpan') ;
%         end
        
%         function out = get.NPlots(self)
%             out = length(self.Scopes);
%         end
                
%         % Need to override the decodeProperties() method supplied by
%         % ws.Encodable() to get correct behavior when the number of
%         % scopes changes.
%         function decodeProperties(self, propSet)
%             % Sets the properties in self to the values encoded in propSet.
%             % Returns the _old_ property values from self in
%             % originalValues.
%             
%             assert(isstruct(propSet));
%             
%             % Need to clear the existing scopes first
%             self.removeScopes();
%             
%             % Now call the superclass method
%             %originalValues=self.decodeProperties@ws.Encodable(propSet);  % not _really_ the originalValues, but I don't think it matters...
%             self.decodeProperties@ws.Encodable(propSet);  % not _really_ the originalValues, but I don't think it matters...
% 
%             % Update the view
%             self.broadcast('NPlotsMayHaveChanged');
%         end  % function
        
%         function didSetScopeIsVisibleWhenDisplayEnabled(self)
%             self.broadcast('DidSetScopeIsVisibleWhenDisplayEnabled');
%         end
    end  % pulic methods block
    
%     methods (Access = protected)        
%         % Need to override the decodeUnwrappedEncodingCore_() method supplied
%         % by ws.Encodable() to get correct behavior when the number of
%         % scopes changes.
%         function decodeUnwrappedEncodingCore_(self, encoding)            
%             % Need to clear the existing scopes first
%             self.removeScopes_();
%             
%             % Now call the superclass method
%             self.decodeUnwrappedEncodingCore_@ws.Encodable(encoding);
% 
%             % Update the view
%             %self.broadcast('NPlotsMayHaveChanged');  % do I need this?
%         end  % function        
%     end  % protected methods block
    
    methods      
        % Allows access to protected and protected variables from ws.Encodable.
        function out = getPropertyValue_(self, name)            
            out = self.(name);
        end
        
        % Allows access to protected and protected variables from ws.Encodable.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end  % function        
    end  % protected methods
    
    methods
        function sanitizePersistedStateGivenChannelCounts_(self, nAIChannels, nDIChannels)
            %nAIChannels = self.Parent.Acquisition.NAnalogChannels ;
            self.IsAnalogChannelDisplayed_ = ws.sanitizeRowVectorLength(self.IsAnalogChannelDisplayed_, nAIChannels, true) ;
            self.AreYLimitsLockedTightToDataForAnalogChannel_ = ...
                ws.sanitizeRowVectorLength(self.AreYLimitsLockedTightToDataForAnalogChannel_, nAIChannels, false) ;
            self.YLimitsPerAnalogChannel_ = ...
                ws.Display.sanitizeYLimitsArrayLength(self.YLimitsPerAnalogChannel_, nAIChannels, [-10 +10]') ;
            self.PlotHeightFromAnalogChannelIndex_  = ...
                ws.sanitizeRowVectorLength(self.PlotHeightFromAnalogChannelIndex_, nAIChannels, 1) ;
            
            %nDIChannels = self.Parent.Acquisition.NDigitalChannels ;
            self.IsDigitalChannelDisplayed_ = ws.sanitizeRowVectorLength(self.IsDigitalChannelDisplayed_, nDIChannels, true) ;            
            self.PlotHeightFromDigitalChannelIndex_  = ...
                ws.sanitizeRowVectorLength(self.PlotHeightFromDigitalChannelIndex_, nDIChannels, 1) ;

            % The analog row indices have to be fixed using
            % knowledge of the digital row indices, and vice-versa
            [self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_] = ...
                ws.Display.sanitizeRowIndices(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_, nAIChannels, nDIChannels) ;
        end
        
        function synchronizeTransientStateToPersistedStateHelper(self)
            self.updateMappingsFromPlotIndices_() ;            
            %self.clearData_() ;  % This will ensure that the size of YData is appropriate
            %self.broadcast('ClearData') ;
        end                
    end
    
    methods
        function sanitizePersistedState_(self) %#ok<MANU>
            % This method should perform any sanity-checking that might be
            % advisable after loading the persistent state from disk.
            % This is often useful to provide backwards compatibility
            
%             nAIChannels = self.Parent.Acquisition.NAnalogChannels ;
%             self.IsAnalogChannelDisplayed_ = ws.sanitizeRowVectorLength(self.IsAnalogChannelDisplayed_, nAIChannels, true) ;
%             self.AreYLimitsLockedTightToDataForAnalogChannel_ = ...
%                 ws.sanitizeRowVectorLength(self.AreYLimitsLockedTightToDataForAnalogChannel_, nAIChannels, false) ;
%             self.YLimitsPerAnalogChannel_ = ...
%                 ws.Display.sanitizeYLimitsArrayLength(self.YLimitsPerAnalogChannel_, nAIChannels, [-10 +10]') ;
%             self.PlotHeightFromAnalogChannelIndex_  = ...
%                 ws.sanitizeRowVectorLength(self.PlotHeightFromAnalogChannelIndex_, nAIChannels, 1) ;
%             
%             nDIChannels = self.Parent.Acquisition.NDigitalChannels ;
%             self.IsDigitalChannelDisplayed_ = ws.sanitizeRowVectorLength(self.IsDigitalChannelDisplayed_, nDIChannels, true) ;            
%             self.PlotHeightFromDigitalChannelIndex_  = ...
%                 ws.sanitizeRowVectorLength(self.PlotHeightFromDigitalChannelIndex_, nDIChannels, 1) ;
% 
%             % The analog row indices have to be fixed using
%             % knowledge of the digital row indices, and vice-versa
%             [self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_] = ...
%                 ws.Display.sanitizeRowIndices(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_, nAIChannels, nDIChannels) ;
        end
        
        function synchronizeTransientStateToPersistedState_(self) %#ok<MANU>
%             self.updateMappingsFromPlotIndices_() ;            
%             %self.clearData_() ;  % This will ensure that the size of YData is appropriate
%             %self.broadcast('ClearData') ;
        end        
        
    end  % protected methods block
    
%     methods (Access=protected)    
%         function disableAllBroadcastsDammit_(self)
%             self.disableBroadcasts() ;
%         end
%         
%         function enableBroadcastsMaybeDammit_(self)
%             self.enableBroadcastsMaybe() ;
%         end
%     end  % protected methods block
    
    methods (Access=protected)    
%         function setIsEnabledImplementation_(self, newValue)
%             if isscalar(newValue) && (islogical(newValue) || (isnumeric(newValue) && (newValue==1 || newValue==0))) ,
%                 self.IsEnabled_ = logical(newValue) ;
%                 didSucceed = true ;
%             else
%                 didSucceed = false ;
%             end
%             if ~didSucceed ,
%                 error('ws:invalidPropertyValue', ...
%                       'IsEnabled must be a scalar, and must be logical, 0, or 1') ;
%             end
%         end
        
        function updateMappingsFromPlotIndices_(self)
            isAnalogChannelDisplayed = self.IsAnalogChannelDisplayed_ ;
            isDigitalChannelDisplayed = self.IsDigitalChannelDisplayed_ ;
            nAnalogChannels = length(isAnalogChannelDisplayed) ;
            nDigitalChannels = length(isDigitalChannelDisplayed) ;
            isAnalogFromChannelIndex = horzcat( true(1,nAnalogChannels), false(1,nDigitalChannels) ) ;
            isDisplayedFromChannelIndex = horzcat(isAnalogChannelDisplayed, isDigitalChannelDisplayed) ;
            rowIndexFromChannelIndex = horzcat(self.RowIndexFromAnalogChannelIndex_, self.RowIndexFromDigitalChannelIndex_) ;
            channelIndexFromPlotIndex = ws.Display.computeChannelIndexFromPlotIndexMapping(rowIndexFromChannelIndex, isDisplayedFromChannelIndex) ;
            isAnalogFromPlotIndex = isAnalogFromChannelIndex(channelIndexFromPlotIndex) ;
            channelIndexWithinTypeFromPlotIndex = ...
                arrayfun(@(channelIndex)(ws.fif(channelIndex>nAnalogChannels,channelIndex-nAnalogChannels,channelIndex)), channelIndexFromPlotIndex) ;
            % Determine the channel index -> plot index mapping
            plotIndexFromChannelIndex = ws.sortedOrderLeavingNansInPlace(rowIndexFromChannelIndex) ;
            % Finally, set the state variables that we need to set
            self.IsAnalogFromPlotIndex_ = isAnalogFromPlotIndex ;
            self.ChannelIndexWithinTypeFromPlotIndex_ = channelIndexWithinTypeFromPlotIndex ;            
            self.ChannelIndexFromPlotIndex_ = channelIndexFromPlotIndex ;
            self.PlotIndexFromChannelIndex_ = plotIndexFromChannelIndex ;
        end  % function
        
    end  % protected methods block    
    
    methods (Static=true)
        function y = sanitizeYLimitsArrayLength(x, targetLength, defaultValue)
            % If x is 2xtargetLength, with all(x(1,:)<x(2,:)), return x.
            % Otherwise, massage x in various ways to make the result
            % 2xtargetLength, with all(y(1,:)<y(2,:)).
            [nRowsOriginal,nColsOriginal] = size(x) ;
            
            % As a first step, fix the shape, so that yProto is 2xn, with
            % all(yProto(1,:)<yProto(2,:)).
            if nRowsOriginal==2 && nColsOriginal==2 ,
                if all(x(1,:)<x(2,:)) ,
                    % Everything looks good.
                    yProto = x ;
                elseif all(x(:,1)<x(:,2)),
                    % if each row has the first el less than
                    % the second, can fix things by transposing.
                    yProto = x' ;                    
                else
                    % WTF?  Force elements to be in right order
                    yProto = ws.prewashYLimitsArray(x) ;
                end                    
            elseif nRowsOriginal==2 ,
                if all(x(1,:)<x(2,:)) ,
                    % Everything looks good.
                    % This should be the common case, want it to be fast, and hopefully
                    % involve no copying...                    
                    yProto = x ;
                else
                    yProto = ws.prewashYLimitsArray(x) ;
                end
            elseif nColsOriginal==2 ,
                yProto = ws.prewashYLimitsArray(x') ;
            else
                % Just use the default value, repeated
                yProto = repmat(defaultValue,[1 targetLength]) ;
            end
            
            % At this point yProto is 2xn, for some n, with 
            % all( yProto(1,:)<yProto(2,:) ).
            nCols = size(yProto,2) ;
            if nCols>targetLength ,
                y = yProto(:,targetLength) ;
            elseif nCols<targetLength ,
                nNewCols = targetLength-nCols ;
                y = horzcat(yProto, repmat(defaultValue,[1 nNewCols])) ;
            else
                % yProto has the right number of cols
                y = yProto ;
            end               
        end        
        
        function y = prewashYLimitsArray(x)
            % x must be 2xn.  Makes sure each col has the first element
            % strictly less than the second.
            n = size(x,2) ;
            y = zeros(2,n) ;
            for i = 1:n ,
                if x(1,i)<x(2,i) ,
                    y(:,i) = x(:,i) ;
                elseif x(1,i)>x(2,i) ,
                    y(:,i) = flipud(x(:,i)) ;
                else
                    % both els equal, so just add/subtract one to make a range
                    y(:,i) = x(1,i) + [-1 +1]' ;
                end
            end
        end
    
        function [newRowIndexFromAnalogChannelIndex, newRowIndexFromDigitalChannelIndex] = ...
                renormalizeRowIndices(rowIndexFromAnalogChannelIndex, rowIndexFromDigitalChannelIndex)
            % Used, e.g. after channel deletion, to maintain the ordering
            % of the row indices, but eliminate any gaps, so that they go
            % from 1 to the number of channels.
            nAIChannels = length(rowIndexFromAnalogChannelIndex) ;
            %nDIChannels = length(rowIndexFromDigitalChannelIndex) ;
            %nChannels = nAIChannels + nDIChannels ;
            rowIndexFromChannelIndex = horzcat(rowIndexFromAnalogChannelIndex, rowIndexFromDigitalChannelIndex) ;  % this may have gaps in the ordering
            newRowIndexFromChannelIndex = ws.sortedOrder(rowIndexFromChannelIndex) ;
            %[~,channelIndexFromRowIndex] = sort(rowIndexFromChannelIndex) ;
            %newRowIndexFromChannelIndex(channelIndexFromRowIndex) = 1:nChannels ;
            newRowIndexFromAnalogChannelIndex = newRowIndexFromChannelIndex(1:nAIChannels) ;
            newRowIndexFromDigitalChannelIndex = newRowIndexFromChannelIndex(nAIChannels+1:end) ;
        end
        
        function [newRowIndexFromAnalogChannelIndex, newRowIndexFromDigitalChannelIndex] = ...
                sanitizeRowIndices(rowIndexFromAnalogChannelIndex, rowIndexFromDigitalChannelIndex, nAIChannels, nDIChannels)
            protoNewRowIndexFromAnalogChannelIndex = ws.sanitizeRowVectorLength(rowIndexFromAnalogChannelIndex, nAIChannels, +inf) ;
            protoNewRowIndexFromDigitalChannelIndex = ws.sanitizeRowVectorLength(rowIndexFromDigitalChannelIndex, nDIChannels, +inf) ;
            [newRowIndexFromAnalogChannelIndex, newRowIndexFromDigitalChannelIndex] = ...
                ws.Display.renormalizeRowIndices(protoNewRowIndexFromAnalogChannelIndex, protoNewRowIndexFromDigitalChannelIndex) ;
        end

        function result = computeChannelIndexFromPlotIndexMapping(rowIndexFromChannelIndex, isDisplayedFromChannelIndex)
            % Computes the mapping from plot index to channel index, given
            % the relevant inputs.  "rows" here refers to rows in the
            % dialog box that the user uses to order the channels for
            % display.  If isDisplayedFromChannelIndex is all true, then
            % the result is simply the inverse permutation of
            % rowIndexFromChannelIndex.
            channelIndexFromRowIndex = ws.invertPermutation(rowIndexFromChannelIndex) ;
            isDisplayedFromRowIndex(rowIndexFromChannelIndex) = isDisplayedFromChannelIndex ;
            result = channelIndexFromRowIndex(isDisplayedFromRowIndex) ;
        end  % function
    end  % static methods block

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
    
    methods
        function result = get.IsEnabled(self)
            result = self.IsEnabled_ ;
        end
        
        function set.IsEnabled(self, value)
            self.IsEnabled_ = value ;
        end
        
        function addData(self, t, recentScaledAnalogData, recentRawDigitalData, nActiveDIChannels, sampleRate, xSpan)
            % t is a scalar, the time stamp of the scan *just after* the
            % most recent scan.  (I.e. it is one dt==1/fs into the future.
            % Queue Doctor Who music.)

            % Get the uint8/uint16/uint32 data out of recentRawDigitalData
            % into a matrix of logical data, then convert it to doubles and
            % concat it with the recentScaledAnalogData, storing the result
            % in yRecent.
            %wsModel = self.Model_ ;
            %nActiveDIChannels = wsModel.getNActiveDIChannels() ;
            if nActiveDIChannels==0 ,
                yRecent = recentScaledAnalogData ;
            else
                % Might need to write a mex function to quickly translate
                % recentRawDigitalData to recentDigitalData.
                nScans = size(recentRawDigitalData,1) ;                
                recentDigitalData = zeros(nScans,nActiveDIChannels) ;
                for j = 1:nActiveDIChannels ,
                    recentDigitalData(:,j) = bitget(recentRawDigitalData,j) ;
                end
                % End of code that might need to mex-ify
                yRecent = horzcat(recentScaledAnalogData, recentDigitalData) ;
            end
            
            % Compute a timeline for the new data            
            nNewScans = size(yRecent, 1) ;
            %sampleRate = wsModel.AcquisitionSampleRate ;
            dt = 1/sampleRate ;  % s
            t0 = t - dt*nNewScans ;  % timestamp of first scan in newData
            xRecent = t0 + dt*(0:(nNewScans-1))' ;
            
            % Figure out the downsampling ratio
            xSpanInPixels = self.XSpanInPixels_ ;
            r = ws.ratioSubsampling(dt, xSpan, xSpanInPixels) ;
            
            % Downsample the new data
            [xForPlottingNew, yForPlottingNew] = ws.minMaxDownsampleMex(xRecent, yRecent, r) ;            
            
            % deal with XData
            xAllOriginal = self.XData_ ;  % these are already downsampled
            yAllOriginal = self.YData_ ;            
            
            % Concatenate the old data that we're keeping with the new data
            xAllProto = vertcat(xAllOriginal, xForPlottingNew) ;
            yAllProto = vertcat(yAllOriginal, yForPlottingNew) ;
            
            % Trim off scans that would be off the screen anyway
            doKeepScan = (self.XOffset_<=xAllProto) ;
            xNew = xAllProto(doKeepScan) ;
            yNew = yAllProto(doKeepScan,:) ;

            % Commit the data to self
            self.XData_ = xNew ;
            self.YData_ = yNew ;            
        end  % function        
        
        function updateTraces(self, scaledAnalogData, digitalDataAsUint, cachedDigitalSignalCount, t, xSpan, sampleRate)
            % t is a scalar, the time stamp of the scan *just after* the
            % most recent scan.  (I.e. it is one dt==1/fs into the future.
            % Queue Doctor Who music.)

            %wsModel = self.Model_ ;
            %scaledAnalogData = wsModel.getAIDataFromCache() ;
            %[digitalDataAsUint, cachedDigitalSignalCount] = wsModel.getDIDataFromCache() ;
            %t = wsModel.getTimestampsForDataInCache() ;
            
            % Get the uint8/uint16/uint32 data out of recentRawDigitalData
            % into a matrix of logical data, then convert it to doubles and
            % concat it with the recentScaledAnalogData, storing the result
            % in yRecent.
            %display = wsModel.Display ;
            digitalDataAsLogical = ws.logicalColumnsFromUintColumn(digitalDataAsUint, cachedDigitalSignalCount) ;
            y = horzcat(scaledAnalogData, digitalDataAsLogical) ;  % horzcat will convert logical to double
            
            % Figure out the downsampling ratio
            xSpanInPixels = self.XSpanInPixels_ ;
            dt = 1/sampleRate ;
            r = ws.ratioSubsampling(dt, xSpan, xSpanInPixels) ;
            
            % Downsample the new data
            [xForPlotting, yForPlotting] = ws.minMaxDownsampleMex(t, y, r) ;            
            
            % Trim off scans that would be off the screen anyway
            doKeepScan = (self.XOffset_ <= xForPlotting) ;
            xNew = xForPlotting(doKeepScan) ;
            yNew = yForPlotting(doKeepScan,:) ;

            % Commit the data to self
            self.XData_ = xNew ;
            self.YData_ = yNew ;
        end  % function        
        
        function result = get.XSpanInPixels(self)
            result = self.XSpanInPixels_ ;
        end        
        
        function set.XSpanInPixels(self, newValue)
            self.XSpanInPixels_ = newValue ;            
        end
        
        function result = get.YData(self)
            result = self.YData_ ;
        end        
        
        function result = get.XData(self)
            result = self.XData_ ;
        end
    end  % public methods block

end
