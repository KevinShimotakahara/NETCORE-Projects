
% *********************************************************************** %
% QLearning scheduler:
% Performs the resource allocation for DL / UL transmission. It can be used
% by the eNB, or the SBS.
% *********************************************************************** %

classdef Random < handle
    
    %% QLearning scheduler properties
    properties
        agentNode   % node that runs this scheduler instance
        
        % QL parameters
        alpha           
        gamma
        epsilon
        qTable          % Q-table
        ndelayStates
        actionSet
        cqiReward
        delayReward
        reward
        nsfExplore      % number of subframes used for exploration
        nsfExploreCounter
        
        % Algorithm parameters
        beta
        targetDelay     % Target utility (Delay)
        cqiMin
        cqiIdeal
        delayIdeal
        delayMax
        nrbsRBG         % number of RBs per RBG
        explTime
    end
    
    methods
        %% Constructor
        function obj = Random
            obj.agentNode = [];
            
            obj.alpha = 0.5;
            obj.gamma = 0.9;
            obj.epsilon = 0.9;           
            obj.qTable = [];
            obj.ndelayStates = 3;
            obj.cqiReward = [];
            obj.delayReward = [];
            obj.reward = [];
            obj.actionSet = [];
            obj.nsfExplore = 50;
            obj.nsfExploreCounter = 0;
            
            obj.beta = 0.5;                 % 1: delay, 0: cqi, 0.5: Both
            obj.targetDelay = [];   
            obj.nrbsRBG = 10;
            obj.cqiMin = 3;
            obj.cqiIdeal = 9;
            obj.delayIdeal = 10;              % msec
            obj.delayMax = 100;             % msec
            obj.explTime = 3500;
        end
        % *************************************************************** %
        
        %% Scheduler initialization
        function initScheduler_fn(obj, ~, ~, delayIdeal, delayMax, cqiMin, cqiIdeal, node, beta)
            obj.agentNode = node;
            
            % init the qTable and action-set
            obj.actionSet = obj.actionSet_fn();
            nActions = length(obj.actionSet);
            obj.qTable = zeros(2*obj.ndelayStates, nActions);
            
            % init the targetDelay
            obj.delayIdeal = delayIdeal;
            obj.delayMax = delayMax;
            
            % init the cqi target
            obj.cqiMin = cqiMin;
            obj.cqiIdeal = cqiIdeal;
            
            % init beta
            obj.beta = beta;
        end
        % *************************************************************** %
        
        %% Generate the action set
        function [actionSet] = actionSet_fn(obj)
            
            subframes1 = '0000000011111111111111111111111111111111';
            subframes2 = '0000000000000000000000001111111111111111';
            subframes3 = '0000000011111111111111110000000000000000';

            % RBG = 5 RBs
            actionSet = {};            
            for i = 1:3
                actionSet = [actionSet; {i, subframes1}];
            end
            
            for i = 1:3
                actionSet = [actionSet; {i, subframes2}];
            end
            
            for i = 1:3
                actionSet = [actionSet; {i, subframes3}];
            end
            
            actionSet = [actionSet; {0, ''}];
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
        
        %% Scheduling algorithm (Q-Learning)
        function [users, grid, timePattern] = schedule_fn(obj)   
            global gClk
            global logData
            
            nActions = length(obj.actionSet);

            actIdx = randi([1, nActions], 1, 1);
            action = obj.actionSet(actIdx, :);
            
            % Save the action
            obj.agentNode.sltxConfig.action = action;
            
            % update the outputs (users, grid)
            users = obj.agentNode.nodeID;        % ID of the user scheduled.
            grid = obj.maptoGrid_fn(action(1));
            timePattern = action(2);
            
            % Save convergence results
            convInfo.nodeType = obj.agentNode.nodeType;
            convInfo.nodeID = obj.agentNode.nodeID;
            convInfo.nSF = gClk;
            convInfo.reward = 0;
            convInfo.cqi = obj.agentNode.sltxConfig.cqiBuffer(1);
            convInfo.delay = obj.agentNode.sltxConfig.recentDelay;
            convInfo.remDelay = 0;
            convInfo.prState = 0;
            convInfo.newState = 0;
            convInfo.action = actIdx;
            convInfo.qValue = 0;
            
            logData.saveconvLog_fn(convInfo);
        end
        % *************************************************************** %
        
        %% Map action to grid
        function [grid] = maptoGrid_fn(obj, action)
            grid = zeros(100, 1);
            action = cell2mat(action);
            
            if(action ~= 0)
                for i = 1:length(action)
                   gridPos = ((action(i)-1)*obj.nrbsRBG+1) : ((action(i))*obj.nrbsRBG);
                   grid(gridPos) = obj.agentNode.nodeID;
                end
            end
        end
        % *************************************************************** %

        %% state table
        function [state] = observeState_fn(obj, cqi, delay)
            if(cqi < obj.cqiMin)
                index1 = 1;
            else
                index1 = 2;
            end
            
            interval = ((obj.delayMax - obj.delayIdeal)/(obj.ndelayStates-2));
            delayStates = obj.delayIdeal:interval:obj.delayMax;
            
            if(isempty(find(delayStates>=delay, 1)))
                index2 = 5;
            else
                index2 = find(delayStates>=delay, 1);
            end
            
            state = (index1-1)*obj.ndelayStates + index2;
        end
        % *************************************************************** %
    end
end






