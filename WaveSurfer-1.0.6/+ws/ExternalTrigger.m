classdef ExternalTrigger < ws.Model %& ws.HasPFIIDAndEdge  % & matlab.mixin.Heterogeneous  (was second in list)
    % A class that represents a trigger destination, i.e. a digital input
    % to the daq board that could potentially be used to trigger
    % acquisition, etc.  A trigger destination has a name, a device ID, 
    % a PFI identifier (the zero-based PFI index),
    % and an edge type (rising or falling).  ALT, 2014-05-24
    
    properties (Dependent=true)
        Name
        DeviceName  % the NI device ID string, e.g. 'Dev1'
        PFIID  
            % If the destination is user-defined, this indicates the PFI
            % line to be used as an input. If the destination is automatic,
            % this indicates the PFI line to be used as an _output_.
        Edge
            % Whether rising edges or falling edges constitute a trigger
            % event.
        IsMarkedForDeletion
    end
    
    properties (Access=protected)
        Name_
        DeviceName_
        PFIID_
        Edge_
    end
    
    properties (Access=protected, Transient=true)
        IsMarkedForDeletion_
    end
    
    methods
        function self=ExternalTrigger()
            %self@ws.Model() ;  % ignore parent arg
            self.Name_ = 'Destination';
            self.DeviceName_ = '' ;
            self.PFIID_ = 0;
            self.Edge_ = 'rising';  
            self.IsMarkedForDeletion_ = false ;
        end
        
        function value=get.Name(self)
            value=self.Name_;
        end
        
        function value=get.DeviceName(self)
            value=self.DeviceName_;
        end
        
        function set.DeviceName(self, deviceName)  % should only be called by triggering subsystem
            self.DeviceName_ = deviceName ;
        end        
        
        function value=get.PFIID(self)
            value=self.PFIID_;
        end
        
        function value=get.Edge(self)
            value=self.Edge_;
        end
        
        function set.Name(self, value)
            if ws.isString(value) && ~isempty(value) ,
                self.Name_ = value ;
            else
                error('ws:invalidPropertyValue', ...
                      'Name must be a nonempty string');                  
            end                    
        end
        
        function set.PFIID(self, value)
            if isnumeric(value) && isscalar(value) && isreal(value) && value==round(value) && value>=0 ,
                value = double(value) ;
                self.PFIID_ = value ;
            else
                error('ws:invalidPropertyValue', ...
                      'PFIID must be a (scalar) nonnegative integer');                  
            end                    
        end
        
        function set.Edge(self, value)
            if ws.isAnEdgeType(value) ,
                self.Edge_ = value;
            else
                error('ws:invalidPropertyValue', ...
                      'Edge must be ''rising'' or ''falling''');                  
            end                                        
        end  % function 
        
        function value = get.IsMarkedForDeletion(self)
            value = self.IsMarkedForDeletion_ ;
        end

        function set.IsMarkedForDeletion(self, value)
            if (islogical(value) || isnumeric(value)) && isscalar(value) ,
                self.IsMarkedForDeletion_ = logical(value) ;
            else
                error('ws:invalidPropertyValue', ...
                      'IsMarkedForDeletion must be a truthy scalar');                  
            end                    
        end
        
    end  % methods
    
    methods 
        function out = getPropertyValue_(self, name)
            out = self.(name);
        end  % function
        
        % Allows access to protected and protected variables from ws.Encodable.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end  % function
    end    
    
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
