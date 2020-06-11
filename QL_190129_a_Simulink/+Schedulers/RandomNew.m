
% *********************************************************************** %
% QLearning scheduler:
% Performs the resource allocation for DL / UL transmission. It can be used
% by the eNB, or the SBS.
% *********************************************************************** %

classdef RandomNew < handle
    
    %% QLearning scheduler properties
    properties
        agentNode   % node that runs this scheduler instance
        nrbsRBG         % number of RBs per RBG
        schedType
        itrp
    end
    
    methods
        %% Constructor
        function obj = RandomNew
            obj.agentNode = [];        
            obj.nrbsRBG = 5;
            obj.schedType = 'RandomNew';
        end
        % *************************************************************** %
        
        %% Scheduler initialization
        function initScheduler_fn(obj,~, ~,  ~, ~, ~, ~, node,~)
            obj.agentNode = node;
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
        
        %% Scheduling algorithm (TM-2 LTE R12)
        function [users, nprb, timePattern,itrp] = schedule_fn(obj)   
            global gClk
        
            users = obj.agentNode.nodeID;
            timePattern = '0000000011111111111111111111111111111111';
            
            %for now, set CQI to 9 always
            obj.agentNode.sltxConfig.cqiBuffer = 9;
            
            %Get imcs based on CQI feedback (if this is annoying, then just
            %keep constant CQI
            imcsTable = [-1 0 0 2 4 6 8 11 13 16 18 21 23 25 27 27];
            if(isempty(obj.agentNode.sltxConfig.cqiBuffer))
                imcs = imcsTable(4);
            elseif(obj.agentNode.sltxConfig.cqiBuffer == 0)
                imcs = imcsTable(obj.agentNode.sltxConfig.cqiBuffer+2);
            else
                imcs = imcsTable(obj.agentNode.sltxConfig.cqiBuffer+1);
            end
            
            %Get itbs based on imcs
            [itbs, ~, ~] = lteMCS(imcs);
            
            %Calculate desired TBS based on buffer status
            RLCbufferBits = 0;
            for k = 1:length(obj.agentNode.RLCtxEntity.TXbuff)
                RLCbufferBits = RLCbufferBits + length(obj.agentNode.RLCtxEntity.TXbuff{k});
            end
            
            %See what nprb is for 1-8 TBs
            i = 1;
            timeFreqSpace = {};
            for nTB = [1 2 4 8]
                desiredPDUsize = ceil(RLCbufferBits/nTB);
                if desiredPDUsize == 0
                    %don't need to schedule anything
                    obj.nrbsRBG = 0;
                    break
                    %grid = zeros(25,1);
                else
                    MACheaderSize = 32;

                    TBSdesired = desiredPDUsize + MACheaderSize;

                    %Loop through tbs table for the row corresponding to itbs
                    %found, find first nprb that can result in tbs large enough; if
                    %buffer is too large, must settle with nprb = entire bandwidth
                    %of PSSCH resource pool
                    SLPRBs = 25;
                    for k = 1:SLPRBs
                        if lteTBS(k, itbs) >= TBSdesired 
                            timeFreqSpace{i} = [nTB, k];
                            i = i+1;
                            break
                        end
                        if k == SLPRBs && nTB == 8
                            timeFreqSpace{i} = [nTB, k];
                            i = i+1;
                            break
                        end
                    end
                end
            end
            
            if gClk < 40
                obj.nrbsRBG = 4;
                itrp = 106;
                nprb = obj.nrbsRBG;
            elseif ~isempty(timeFreqSpace)
               lenTFspace = length(timeFreqSpace);
               switch lenTFspace
                   case 1
                       %only nTB = 8 is an option
                       timeFreqMagnitudeChosen = timeFreqSpace{1};
                   case 2
                       %nTB = 8 or nTB = 4 are options; must bias
                       %selection of timeFreqMagnitude to reflect that
                       %there are more event elements for nTB = 4                        roll = randi([36 106]);
                       roll = randi([36 106]);
                       if roll == 106
                           timeFreqMagnitudeChosen = timeFreqSpace{2};
                       else
                           timeFreqMagnitudeChosen = timeFreqSpace{1};
                       end
                   case 3
                       %nTB = 8 or nTB = 4 or nTB = 2 are options; must bias
                       %selection of timeFreqMagnitude to reflect that
                       %there are more event elements for nTB = 4                        roll = randi([36 106]);
                       roll = randi([8 106]);
                       if roll == 106
                           timeFreqMagnitudeChosen = timeFreqSpace{3};
                       elseif roll >= 36
                           timeFreqMagnitudeChosen = timeFreqSpace{2};
                       else
                           timeFreqMagnitudeChosen = timeFreqSpace{1};
                       end
                   case 4
                       %nTB = 8 or nTB = 4 or nTB = 2 are options; must bias
                       %selection of timeFreqMagnitude to reflect that
                       %there are more event elements for nTB = 4                        roll = randi([36 106]);
                       roll = randi([0 106]);
                       if roll == 106
                           timeFreqMagnitudeChosen = timeFreqSpace{4};
                       elseif roll >= 36
                           timeFreqMagnitudeChosen = timeFreqSpace{3};
                       elseif roll >= 8
                           timeFreqMagnitudeChosen = timeFreqSpace{2};
                       else
                           timeFreqMagnitudeChosen = timeFreqSpace{1};
                       end
                   otherwise
                       roll = 106;
               end
                
               obj.nrbsRBG = timeFreqMagnitudeChosen(2);
               itrp = roll;
               nprb = obj.nrbsRBG;          
            else
                %no packets to tx
                nprb = 0;
                itrp = 0;
            end
            obj.itrp = itrp;
        end
        % *************************************************************** %
%         %% Map action to grid (may not need if I'm going to do frequency hopping)
%         function [grid] = maptoGrid_fn(obj, idx)
%             grid = zeros(25, 1);
%             grid(idx:idx+obj.nrbsRBG-1) = obj.agentNode.nodeID;
%         end
%         % *************************************************************** %
    end
end






