
% *********************************************************************** %
% QLearning scheduler:
% Performs the resource allocation for DL / UL transmission. It can be used
% by the eNB, or the SBS.
% *********************************************************************** %

classdef Fixed < handle
    
    %% QLearning scheduler properties
    properties
        agentNode   % node that runs this scheduler instance
        
        ueQueue
    end
    
    methods
        %% Constructor
        function obj = Fixed
            obj.agentNode = [];
            obj.ueQueue = [];
        end
        % *************************************************************** %
        
        %% Scheduler initialization
        function initScheduler_fn(obj)
            
        end
        % *************************************************************** %
        
        %% Add UD
        function addUser_fn(obj, udTemp)
            
            switch(udTemp.nodeType)
                case 'UE'
                    existing = find(udTemp == obj.ueQueue);
                    
                    if(isempty(existing))
                        obj.ueQueue = [obj.ueQueue, udTemp];
                    end
            end
        end
        % *************************************************************** %
        
        %% Remove user from the queue
        function removeUser_fn(obj, udTemp)
            
            switch(udTemp.nodeType)
                case 'UE'
                    userID = find(obj.ueQueue == udTemp);
                    
                    if(~isempty(userID))
                        obj.ueQueue(userID) = [];
                    end
            end
        end
        % *************************************************************** %
        
        %% Scheduling algorithm
        % usersSchd: users scheduled (ids)
        % ulGrid: the resources allocated to every users
        function [usersSchd, ulGrid] = schedule_fn(obj)
            
            global simParameters
            
            if(~isempty(obj.ueQueue))
                udsTemp = obj.ueQueue;

                nUDs = length(udsTemp);
                ulGrid = zeros(simParameters.nRBs, 1);
                nRBsPerUD = floor(simParameters.nRBs / nUDs);

                for u = 1:nUDs
                    ulGrid((u-1)*nRBsPerUD+1 : u*nRBsPerUD) = udsTemp(u).nodeID;
                end

                usersSchd = [udsTemp.nodeID];
            else
                usersSchd = [];
                ulGrid = zeros(simParameters.nRBs, 1);
            end            
            
        end
        % *************************************************************** %
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end