
% *********************************************************************** %
% QLearning scheduler:
% Performs the resource allocation for DL / UL transmission. It can be used
% by the eNB, or the SBS.
% *********************************************************************** %

classdef QLearning < handle
    
    %% QLearning scheduler properties
    properties
        agentNode   % node that runs this scheduler instance
        
        % QL parameters
        alpha           
        gamma
        epsilon
        epsilon0
        zeta            % exploitation parameter
        qTable          % Q-table
        stateVisits
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
        schedType
    end
    
    methods
        %% Constructor
        function obj = QLearning
            obj.agentNode = [];
            obj.schedType = 'QLearning';
            obj.alpha = 0.5;
            obj.gamma = 0.9;
            obj.epsilon = 0.1; 
            obj.epsilon0 = 0.95;
            obj.zeta = 3;
            obj.qTable = [];
            obj.stateVisits = [];
            obj.ndelayStates = 3;
            obj.cqiReward = [];
            obj.delayReward = [];
            obj.reward = [];
            obj.actionSet = [];
            obj.nsfExplore = 50;
            obj.nsfExploreCounter = 0;
            
            obj.beta = 0;                 % 1: delay, 0: cqi, 0.5: Both
            obj.targetDelay = [];   
            obj.nrbsRBG = 4;
            obj.cqiMin = 3;
            obj.cqiIdeal = 9;
            obj.delayIdeal = 10;              % msec
            obj.delayMax = 100;             % msec
            obj.explTime = 1.999 * 1e3;       
        end
        % *************************************************************** %
        
        %% Scheduler initialization
        function initScheduler_fn(obj, ~, ~, delayIdeal, delayMax, cqiMin, cqiIdeal, node, beta)
            obj.agentNode = node;
            
            % init the qTable and action-set
            obj.actionSet = obj.actionSet_fn();
            nActions = length(obj.actionSet);
            obj.qTable = zeros(2*obj.ndelayStates, nActions);
            obj.stateVisits = zeros(2*obj.ndelayStates, 1);
%             obj.epsilon = zeros(2*obj.ndelayStates, 1);
            
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
            
%             subframes1 = '0000000011111111111111111111111111111111';
            subframes2 = '0000000000000000000000001111111111111111';
            subframes3 = '0000000011111111111111110000000000000000';

            % RBG = 10 RBs
            actionSet = {};            
            
            for i = 1:6
                actionSet = [actionSet; {i, subframes2}];
            end
            
            for i = 1:6
                actionSet = [actionSet; {i, subframes3}];
            end
            
%             actionSet = [actionSet; {0, ''}];
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
            
            % Algorithm steps:
            % 1. Compute the recent packet delay and update the reward.
            % 2. Observe the new state.
            % 3. Update the q-value and the q-table.
            % 4. Update the state.
            % 5. Action selection (epsilon-greedy).
            
            nodeTemp = obj.agentNode;
            prAction = nodeTemp.schedInfo.prAction;  % idx of present action
            prState = nodeTemp.schedInfo.prState;    % idx of present state
            nActions = length(obj.actionSet);
            
            if(isempty(nodeTemp.RLCtxEntity.oldest_crtTS))
                peakQDelay = rand;
            else
                peakQDelay = gClk - nodeTemp.RLCtxEntity.oldest_crtTS;
            end
            
            % 1. Update the reward
            % Delay reward
            if(peakQDelay <= obj.delayIdeal)
                obj.delayReward = 1;
            elseif(peakQDelay > obj.delayIdeal && peakQDelay <= obj.delayMax)
                obj.delayReward = (obj.delayIdeal/peakQDelay);
            else
                obj.delayReward = -1;
            end
            
            % CQI reward
            if ~isempty(nodeTemp.sltxConfig.cqiBuffer)
                %Only advance algorithm if data was actually sent last
                %period
                if(nodeTemp.sltxConfig.cqiBuffer(1) >= obj.cqiIdeal)
                    obj.cqiReward = 1;
    %             elseif(nodeTemp.sltxConfig.cqiBuffer(1) > obj.cqiMin && nodeTemp.sltxConfig.cqiBuffer(1) <= obj.cqiIdeal)
    %                 obj.cqiReward = (nodeTemp.sltxConfig.cqiBuffer(1)/obj.cqiIdeal);
                else
                    obj.cqiReward = -1;
                end

                % Total reward
                obj.reward = obj.beta * obj.delayReward + (1 - obj.beta) * obj.cqiReward;

                % 2. Observe the new state.
                newState = obj.observeState_fn(nodeTemp.sltxConfig.cqiBuffer(1), peakQDelay);

                % 3. Update the q-value and the q-table.
                oldqValue = obj.qTable(prState, prAction);

                newqValue = ((1 - obj.alpha) * oldqValue) +...
                                obj.alpha * (obj.reward + obj.gamma * (max(obj.qTable(newState, :))));

                obj.qTable(prState, prAction) = newqValue;
                obj.stateVisits(prState) = obj.stateVisits(prState)+1;

                % 4. Update the state.
                nodeTemp.schedInfo.prState = newState;

                % 5. Action selection (Decaying epsilon-greedy).
                % update epsilon
    %             obj.epsilon(newState) = (1 - obj.epsilon0)^(obj.stateVisits(newState)/(obj.zeta*nActions)) * obj.epsilon0;

                r = rand; % get 1 uniform random number
                epsGreedy = (r <= obj.epsilon);

                if(gClk < obj.explTime)         
                    if(epsGreedy >= 1)   % explore 
                        actIdx = randi([1, nActions], 1, 1);
                        action = obj.actionSet(actIdx, :);
                        nodeTemp.schedInfo.prAction = actIdx;
                    else                 % exploit 
                        [ActVal, actIdx] = max(obj.qTable(newState, :));
                        [idx] = find(obj.qTable(newState, :) == ActVal);
                        if(length(idx) > 1)
                            actIdx = datasample(idx, 1);
                        end
                        action = obj.actionSet(actIdx, :);
                        nodeTemp.schedInfo.prAction = actIdx;
                    end
                else
                    [ActVal, actIdx] = max(obj.qTable(newState, :));
                    [idx] = find(obj.qTable(newState, :) == ActVal);
                    if(length(idx) > 1)
                        actIdx = datasample(idx, 1);
                    end
                    action = obj.actionSet(actIdx, :);
                    nodeTemp.schedInfo.prAction = actIdx;
                end                
                % Save the action
                obj.agentNode.sltxConfig.action = action;

                % 6. update the outputs (users, grid)
                users = nodeTemp.nodeID;        % ID of the user scheduled.
                grid = obj.maptoGrid_fn(action(1));
                timePattern = action(2);

                % Save convergence results
                convInfo.nodeType = obj.agentNode.nodeType;
                convInfo.nodeID = obj.agentNode.nodeID;
                convInfo.nSF = gClk;
                convInfo.reward = obj.reward;
                convInfo.cqi = nodeTemp.sltxConfig.cqiBuffer(1);
                convInfo.delay = nodeTemp.sltxConfig.recentDelay;
                convInfo.remDelay = peakQDelay;
                convInfo.prState = prState;
                convInfo.newState = newState;
                convInfo.action = actIdx;
                convInfo.qValue = newqValue;
                convInfo.epsilon = obj.epsilon;

                logData.saveconvLog_fn(convInfo);
            else
                % return same scheduling decision as last time
                action = obj.agentNode.sltxConfig.action;
                users = nodeTemp.nodeID;        % ID of the user scheduled.
                grid = obj.maptoGrid_fn(action(1));
                timePattern = action(2);
            end
                
        end
        % *************************************************************** %
        
        %% Map action to grid
        function [grid] = maptoGrid_fn(obj, action)
            grid = zeros(25, 1);
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
                index2 = obj.ndelayStates;
            else
                index2 = find(delayStates>=delay, 1);
            end
            
            state = (index1-1)*obj.ndelayStates + index2;
        end
        % *************************************************************** %
    end
end






