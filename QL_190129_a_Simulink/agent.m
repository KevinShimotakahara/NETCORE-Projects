
%% Class
classdef agent < handle
    
    properties
        nodeID
        oldest_crtTS
        cqiBuffer
        schedInfo
        scheduler
        action
    end
    
    methods
        
        %% UD Constructor
        function obj = agent
            obj.nodeID = 0;
            obj.oldest_crtTS = 0;
            obj.cqiBuffer = 9;
            
            obj.schedInfo.prAction = 1;
            obj.schedInfo.prState = 1;
            
            obj.scheduler = [];
            obj.action = [];
        end
        
    end
    
end
    
    