
% *********************************************************************** %
% Dscription: Generic Node class
% *********************************************************************** %

classdef UD < handle
    
    %% UD Properties
    properties
        % Node identification parameters
        nodeID                  % Node ID
        absID                   % absolute ID for every type of device
        nodeType
        nodeDummy
        commDir         % communication direction: tx, or rx
        cellID                  % eNodeB ID
        cellType
        cellAttached            % cell to which UD is attached (struct)
        txUsers                 % transmitters that this user receives from
        rxUsers             % receivers connected to this transmitter
        
        % Geographical data
        disttoCell
        disttoTxer
        disttoRxer
        nodePosition
        udsInterfering
        
        % Transmission parameters
        txPower         % Tx Power in dB
        ud
        ulConfig        % uplink configuration
        dlConfig        % downlink configuration
        sltxConfig      % sidelink configuration
        slrxConfig
        RLCtxEntity
        RLCrxEntity
        
        cqiBuffer       % Buffer for the CQI
        sinrBuffer      % sinr buffer
        
        % traffic data and control
        ultxPkt
        sltxPkt
        
        % layer
        macPadding
        expectedNSAID
        pscchPeriod
        MACheaderSL
        
        % Packet
        avgPktSize
        txBuffer
        rxBuffer                % Rxer buffer for upper layer
        pktNumber                % Packet number
        pktEvents
        arrivTimes
        recentDelay
        
        % Scheduling parameters
        d2dScheduler
        schedInfo
        nsrUE          % variable to keep track of 4 subframes after reciving SR
        nsrD2D
        
        % other
        nTBs            % Variable to keep track of 4 TBs repetition
        txDone
    end
    
    methods
        
        %% UD Constructor
        function obj = UD
            obj.nodeID = [];
            obj.absID = 0;
            obj.nodeType = [];
            obj.nodeDummy = 0;
            obj.commDir = [];
            obj.cellID = [];
            obj.cellType = [];
            obj.cellAttached = [];
            obj.txUsers = [];
            obj.rxUsers = [];
            
            obj.disttoCell = [];
            obj.disttoTxer = [];
            obj.disttoRxer = [];
            obj.nodePosition = [];
            obj.udsInterfering = [];
            
            obj.txPower = [];
            obj.ud = [];
            obj.ulConfig = [];
            obj.dlConfig = [];
            obj.sltxConfig = [];
            obj.slrxConfig = [];
            obj.RLCtxEntity = [];
            obj.RLCrxEntity = [];
            
            obj.cqiBuffer = [];
            obj.sinrBuffer = [];
            
            obj.ultxPkt = [];
            obj.sltxPkt = [];
            
            obj.macPadding = [];
            obj.expectedNSAID = [];
            obj.pscchPeriod = [];
            obj.MACheaderSL = [];
            
            obj.avgPktSize = [];          % 25 bytes
            obj.txBuffer = [];
            obj.rxBuffer = {};
            obj.pktNumber = 0;
            obj.pktEvents = [];
            obj.arrivTimes = [];
            obj.recentDelay = 0;
            
            obj.d2dScheduler = [];
            
            obj.schedInfo = struct;
            obj.schedInfo.prAction = 1;     % idx of present action
            obj.schedInfo.prState = 1;      % idx of present state
            
            obj.nsrUE = 0;
            obj.nsrD2D = 0;
            
            obj.nTBs = 0;
            obj.txDone = 0;
        end
        % *************************************************************** %
        
        %% Scheduler configuration
        function setScheduler_fn(obj, SchType)
            % Select the scheduler
            switch(SchType)
                case 'RoundRobin'
                    obj.d2dScheduler = Schedulers.RoundRobin;
                    
                case 'MaxThroughput'
                    obj.d2dScheduler = Schedulers.MaxThroughput;
                    
                case 'PropFairness'
                    obj.d2dScheduler = Schedulers.PropFairness;
                    
                case 'QLearning'
                    obj.d2dScheduler = Schedulers.QLearning;
                    
                case 'Fixed'
                    obj.d2dScheduler = Schedulers.Fixed;
                    
                case 'Random'
                    obj.d2dScheduler = Schedulers.Random;
                    
                otherwise   % Round Robin is the default scheduler
                    obj.d2dScheduler = Schedulers.RoundRobin;
            end
        end
        % *************************************************************** %
        
        %% Generate all the traffic beforehand
        function allTraffic_fn(obj, lambda)
            
            %global RunTime
            RunTime = 20000;
            
            switch(obj.nodeType)
                % Poisson traffic for UEs
                case 'UE'
                    TTI = 1e-3;
                    T = RunTime*1e-3; % simulation time in second (10 hours)
                    delta = TTI;          % simulation step size in second
                    N = (T/delta);        % number of simulation steps
                    event = zeros(N,1);     % array recording at each step if a "packet" arrived.
                    R = rand(size(event));
                    event(R < (lambda*delta)) = 1;
                    inds = find(event == 1);
                    intArrivtime = diff(inds) * delta;
                    
                    AbsArrivT = fix(cumsum(intArrivtime)*1000);
                    AbsArrivT(AbsArrivT >= RunTime) = [];
                    
                    event = zeros(N, 1);
                    Idx = ceil(AbsArrivT) + 1;
                    Idx = tabulate(Idx);
                    event(Idx(:, 1)) = Idx(:, 2);
                    padding = RunTime - length(event);
                    event = [event; zeros(padding, 1)];
                    
                    obj.pktEvents = event;
                    obj.arrivTimes = AbsArrivT;
                % ------------------------------------------------------- %
                
                % case 'D2D'
                    % if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))
                        % TTI = 1e-3;
                        % T = RunTime*1e-3; % simulation time in second (10 hours)
                        % delta = TTI;          % simulation step size in second
                        % N = (T/delta);        % number of simulation steps
                        % event = zeros(N,1);     % array recording at each step if a "packet" arrived.
                        % R = rand(size(event));
                        % event(R < (lambda*delta)) = 1;
                        % inds = find(event == 1);
                        % intArrivtime = diff(inds) * delta;
                        
                        % AbsArrivT = fix(cumsum(intArrivtime)*1000);
                        % AbsArrivT(AbsArrivT >= RunTime) = [];
                        
                        % event = zeros(N, 1);
                        % Idx = ceil(AbsArrivT)+1;
                        % Idx = tabulate(Idx);
                        % event(Idx(:, 1)) = Idx(:, 2);
                        % padding = RunTime - length(event);
                        % event = [event; zeros(padding, 1)];
                        
                        % obj.pktEvents = event;
                        % obj.arrivTimes = AbsArrivT;
                    % end
                % ------------------------------------------------------- %
                
                case 'DD2D'
                    if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))
                        TTI = 1e-3;
                        T = RunTime*1e-3; % simulation time in second (10 hours)
                        delta = TTI;          % simulation step size in second
                        N = (T/delta);        % number of simulation steps
                        event = zeros(N,1);     % array recording at each step if a "packet" arrived.
                        R = rand(size(event));
                        event(R < (lambda*delta)) = 1;
                        inds = find(event == 1);
                        intArrivtime = diff(inds) * delta;
                        
                        AbsArrivT = fix(cumsum(intArrivtime)*1000);
                        AbsArrivT(AbsArrivT >= RunTime) = [];
                        
                        event = zeros(N, 1);
                        Idx = ceil(AbsArrivT)+1;
                        Idx = tabulate(Idx);
                        event(Idx(:, 1)) = Idx(:, 2);
                        padding = RunTime - length(event);
                        event = [event; zeros(padding, 1)];
                        
                        obj.pktEvents = event;
                        obj.arrivTimes = AbsArrivT;
                    end
                % ------------------------------------------------------- %
            end
        end
        % *************************************************************** %

        %% update interferers list to rx D2D users
        function udsInterfer_fn(obj)
            global UDs
            
            cellRange = obj.cellAttached.cellRadius;
            
            if(strcmp(obj.commDir, "txrx") || strcmp(obj.commDir, "rx"))
                for u = 1:length(UDs)
                    udTemp = UDs(u);
                    
                    distance = sqrt(sum((obj.nodePosition - udTemp.nodePosition).^2));
                    
                    % check interferers only if the user is D2D or UE
                    % transmitter or it is outside the cell.
                    if(udTemp.nodeID == obj.nodeID || strcmp(udTemp.commDir, "rx") || distance > cellRange)
                        continue;
                    else
                        obj.udsInterfering = [obj.udsInterfering, udTemp];
                    end
                end
            end
        end
        % *************************************************************** %
        
        %% Schedule Traffic data generation
        function trafficGen_fn(obj)
            
            global gClk
            
            % Generate new data for UE only (D2D traffic is generated from
            % simulink).
            switch(obj.nodeType)
                case 'UE'
                    txFlag = obj.pktEvents(gClk+1);
                
                    while(txFlag > 0)
                        pktSize = ceil(exprnd(obj.avgPktSize,[1, 1]));  % Packet size
                        obj.pktNumber = obj.pktNumber + 1;
                        remainSize = pktSize;                   % Remaining number of bits To Send

                        % txBuffer
                        ipPkt.pktSize = pktSize;                % Total packet size
                        ipPkt.pktNumber = obj.pktNumber;        % Packet Number
                        ipPkt.crtTS = obj.arrivTimes(obj.pktNumber);  % Packet creation Time Stamp
                        ipPkt.remSize = remainSize;

                        obj.txBuffer    = [ipPkt; obj.txBuffer];

                        txFlag = txFlag - 1;
                    end
                % ------------------------------------------------------- %
                
                % D2D
                % case 'D2D'
                    % if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))
                        % txFlag = obj.pktEvents(gClk+1);

                        % while(txFlag > 0)

                            % obj.pktNumber = obj.pktNumber + 1;

                            % % txBuffer
                            % ipPkt.pktSize = ceil(exprnd(obj.avgPktSize,[1, 1]));    % total packet size  
                            % if(ipPkt.pktSize > 1800)
                                % ipPkt.pktSize = 1800;
                            % end
                            % ipPkt.pktNumber = obj.pktNumber;        % Packet Number
                            % ipPkt.crtTS = obj.arrivTimes(obj.pktNumber);  % Packet creation Time Stamp
                            % ipPkt.remSize = ipPkt.pktSize;
                            % ipPkt.data = randi([0, 1], 1, ipPkt.pktSize);

                            % % create the header
                            % devID = fliplr(de2bi(obj.nodeID, 5));
                            % nPkt = fliplr(de2bi(obj.pktNumber, 24));
                            % crtTS = fliplr(de2bi(gClk, 24));
                            % txTS = fliplr(de2bi(gClk, 24));
                            % ipHeader = [devID, nPkt, crtTS, txTS];
                            % ipData = [ipHeader, ipPkt.data]; 

                            % obj.RLCtxEntity.TXbuff = [obj.RLCtxEntity.TXbuff, {ipData}];

                            % txFlag = txFlag - 1;
                        % end 
                    % end
                % ------------------------------------------------------- %
                
                % Dummy D2D
                case 'DD2D'
                    if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))
                        txFlag = obj.pktEvents(gClk+1);

                        while(txFlag > 0)

                            obj.pktNumber = obj.pktNumber + 1;

                            % txBuffer
                            ipPkt.pktSize = ceil(exprnd(obj.avgPktSize,[1, 1]));    % total packet size  
                            if(ipPkt.pktSize > 1800)
                                ipPkt.pktSize = 1800;
                            end
                            ipPkt.pktNumber = obj.pktNumber;        % Packet Number
                            ipPkt.crtTS = obj.arrivTimes(obj.pktNumber);  % Packet creation Time Stamp
                            ipPkt.remSize = ipPkt.pktSize;
                            ipPkt.data = randi([0, 1], 1, ipPkt.pktSize);

                            % create the header
                            devID = fliplr(de2bi(obj.nodeID, 5));
                            nPkt = fliplr(de2bi(obj.pktNumber, 24));
                            crtTS = fliplr(de2bi(gClk, 24));
                            txTS = fliplr(de2bi(gClk, 24));
                            ipHeader = [devID, nPkt, crtTS, txTS];
                            ipData = [ipHeader, ipPkt.data]; 

                            obj.RLCtxEntity.TXbuff = [obj.RLCtxEntity.TXbuff, {ipData}];

                            txFlag = txFlag - 1;
                        end 
                    end
                % ------------------------------------------------------- %
            end
        end
        % *************************************************************** %
        
        %% This function performs one of the following:
        % ul or sl transmission
        % Send SR on the uplink from UE or from D2D
        function tx_fn(obj)     
            global gClk
            
            switch(obj.nodeType)
                % UE
                case 'UE'
                    if(~isempty(obj.txBuffer))
                        % Check if eNB has done ul RB allocation (Scheduling falg = 1)
                        if(obj.ulConfig.schFlag == 1)
                            % Perform ul tx
                            obj.ultx_fn();
                            
                            obj.ulConfig.schFlag = 0;
                            
                            if(isempty(obj.txBuffer))
                                obj.ulConfig.SR = 0;
                            end
                        else    % Send SR request
                            obj.ulConfig.SR = 1;
                        end
                    else % clear all tx parameters
                        obj.ulConfig.SR = 0;
                        obj.ulConfig.schFlag = 0;
                    end
                % ------------------------------------------------------- %  
                
                % D2D
                case 'D2D'     
                    if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))  
                        if(mod(gClk, 40)==0)
                            % sl scheduling (RB allocation)
                            if(isempty(obj.sltxConfig.cqiBuffer))
                                obj.sltxConfig.cqiBuffer = 1;
                            end
                            if(isempty(obj.sltxConfig.recentDelay))
                                obj.sltxConfig.recentDelay = 10*rand();
                            end
                                    
                            [~, slGrid] = obj.d2dScheduler.schedule_fn();
                            obj.sltxConfig.NPRB = length(find(slGrid == obj.nodeID));

                            % sci message 
                            [sciMessage, TBS, ue] = obj.sci_fn();
                            obj.sltxConfig.sciMessage = sciMessage;
                            obj.sltxConfig.TBS = TBS;
                            obj.sltxConfig.ue = ue;
                            [obj.sltxConfig.pscchsubframes, ~] = obj.pscchPeriod.getPSCCHResources(sciMessage);
                            [obj.sltxConfig.psschsubframes] = obj.pscchPeriod.getPSSCHResources(sciMessage);
                            obj.sltxConfig.txWaveform = [];
                            
                            % other 
                            obj.nTBs = 0;
                        end
                    
                        subframeNo = mod(gClk, 40);
                        
                        if(any(subframeNo == obj.sltxConfig.pscchsubframes))
                            obj.slPhypscch_fn();
                            obj.sltxConfig.phyTB = [];
                        elseif(any(subframeNo == obj.sltxConfig.psschsubframes))
                            if(obj.nTBs == 4)
                                obj.nTBs = 0; 
                                obj.sltxConfig.phyTB = [];
                                obj.sltxConfig.txWaveform = [];
                                obj.sltxConfig.txChInfo = [];
                                obj.sltxConfig.txInfo = [];
                            end
                            
                            if(~isempty(obj.sltxConfig.phyTB))
                                obj.slPhypssch_fn();
                                obj.nTBs = obj.nTBs + 1;
                            elseif(~isempty(obj.RLCtxEntity.TXbuff))
                                if(obj.nTBs == 0)
                                    obj.slRLCMAC_fn();
                                    obj.slPhypssch_fn();
                                    obj.nTBs = obj.nTBs + 1;
                                end
                            end

                        else
                            obj.sltxConfig.phyTB = [];
                            obj.sltxConfig.txWaveform = [];
                            obj.sltxConfig.txChInfo = [];
                            obj.sltxConfig.txInfo = [];
                        end
                    end
                % ------------------------------------------------------- %    
                
                % D2D
                case 'DD2D'     
                    if(strcmp(obj.commDir, 'tx') || strcmp(obj.commDir, 'txrx'))  
                        if(mod(gClk, 40)==0)
                            % sl scheduling (RB allocation)
                            if(isempty(obj.sltxConfig.cqiBuffer))
                                obj.sltxConfig.cqiBuffer = 1;
                            end
                            if(isempty(obj.sltxConfig.recentDelay))
                                obj.sltxConfig.recentDelay = 10*rand();
                            end
                                    
                            [~, slGrid] = obj.d2dScheduler.schedule_fn();
                            obj.sltxConfig.NPRB = length(find(slGrid == obj.nodeID));

                            % sci message 
                            [sciMessage, TBS, ue] = obj.sci_fn();
                            obj.sltxConfig.sciMessage = sciMessage;
                            obj.sltxConfig.TBS = TBS;
                            obj.sltxConfig.ue = ue;
                            [obj.sltxConfig.pscchsubframes, ~] = obj.pscchPeriod.getPSCCHResources(sciMessage);
                            [obj.sltxConfig.psschsubframes] = obj.pscchPeriod.getPSSCHResources(sciMessage);
                            obj.sltxConfig.txWaveform = [];
                            
                            % other 
                            obj.nTBs = 0;
                        end
                    
                        subframeNo = mod(gClk, 40);
                        
                        if(any(subframeNo == obj.sltxConfig.pscchsubframes))
                            obj.slPhypscch_fn();
                            obj.sltxConfig.phyTB = [];
                        elseif(any(subframeNo == obj.sltxConfig.psschsubframes))
                            if(obj.nTBs == 4)
                                obj.nTBs = 0; 
                                obj.sltxConfig.phyTB = [];
                                obj.sltxConfig.txWaveform = [];
                                obj.sltxConfig.txChInfo = [];
                                obj.sltxConfig.txInfo = [];
                            end
                            
                            if(~isempty(obj.sltxConfig.phyTB))
                                obj.slPhypssch_fn();
                                obj.nTBs = obj.nTBs + 1;
                            elseif(~isempty(obj.RLCtxEntity.TXbuff))
                                if(obj.nTBs == 0)
                                    obj.slRLCMAC_fn();
                                    obj.slPhypssch_fn();
                                    obj.nTBs = obj.nTBs + 1;
                                end
                            end

                        else
                            obj.sltxConfig.phyTB = [];
                            obj.sltxConfig.txWaveform = [];
                            obj.sltxConfig.txChInfo = [];
                            obj.sltxConfig.txInfo = [];
                        end
                    end
                % ------------------------------------------------------- %    
            end
        end
        % *************************************************************** %
        
        %% ul/sl tx (RLC, MAC, PHY)
        function ultx_fn(obj)
            
            global gClk
            
            % Steps:
            % 1. Update the ul configuration modified by the eNB (PRBSet, MCS,
            % TBS) for use in RLC, MAC, and PHY layers processing.
            % 2. RLC: Create the SDU from the IP packet in the txBuffer,
            % perform segmentation if needed.
            % 3. MAC: For UE, only HARQ scheduling is needed since the eNB
            % is responsible for RB allocation.
            % 4. PHY: Perform PHY processing (TB, Coding, Modulation,
            % waveform, etc).
            % 5. Apply fading channel to the waveform.
            
            % 1. Update the ul configuration modified by the eNB (PRBSet, MCS,
            % TBS) for use in RLC, MAC, and PHY layers processing.
            
            subframeNo = gClk;
            obj.ulConfig.NSubframe = subframeNo;
%             harqProcessSequence = obj.ulConfig.HARQ.harqProcessSequence;
            harqProcesses = obj.ulConfig.HARQ.harqProcesses;
            
            %             harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);
            harqID = mod(obj.ulConfig.HARQ.harqProcIdx, 8) + 1;
            obj.ulConfig.HARQ.harqProcIdx = harqID;
            obj.ulConfig.HARQ.harqID = harqID;
            
            % Update current HARQ process
            harqProcesses(harqID) = hPUSCHHARQScheduling(harqProcesses(harqID));
            obj.ulConfig.PUSCH.RV = harqProcesses(harqID).rvSeq(harqProcesses(harqID).rvIdx);
            obj.ulConfig.PUSCH.RVSeq = harqProcesses(harqID).rvSeq(harqProcesses(harqID).rvIdx);
            
            % HARQ: check CRC from previous transmission, i.e. is a
            % retransmission required?
            if(harqProcesses(harqID).crc)
                newTrBlkRequired = false; % Retransmit old data
                if (harqProcesses(harqID).rvIdx==1)
                    newTrBlkRequired = true;
                end
            else % Get new data
                newTrBlkRequired = true;
            end
            
            if(newTrBlkRequired)
                fprintf('HARQ number = %d \t New transmission.\n', harqID);
                % 2. RLC: Create the SDU from the IP packet in the txBuffer,
                % perform segmentation if needed.
                macPDU = obj.ulRLC_fn(obj.ulConfig.TBS);
                
                % 3. MAC: For UE, only HARQ scheduling is needed since the eNB
                % is responsible for RB allocation. (Not complete yet).
                phyTB = macPDU;
                harqProcesses(harqID).ulschTransportBlk = phyTB;
                
                % update harq
                harqProcesses(harqID).rvIdx = 1;
                harqProcesses(harqID).decState = [];
                harqProcesses(harqID).crc = 0;
                harqProcesses(harqID).trBlkSize = obj.ulConfig.TBS;
                
                % 4. PHY: Perform PHY processing (TB, Coding, Modulation,
                % waveform, etc).
                [txWaveform, txGrid, codedTrBlock, puschIndices, txInfo,...
                    txChInfo, obj.ulConfig] = ulPHY_fn(obj.ulConfig, harqProcesses(harqID).ulschTransportBlk);
                
                harqProcesses(harqID).codedTrBlkSize = obj.ulConfig.PUSCH.CodedTrBlkSizes(1);
                
                % Get the HARQ ID sequence from 'enbOut' for HARQ processing
                harqProcessSequence = obj.ulConfig.PUSCH.HARQProcessSequence;
                obj.ulConfig.HARQ.harqProcessSequence = harqProcessSequence;
                obj.ulConfig.HARQ.harqProcesses = harqProcesses;
                
                % transmitted packet (ultxPkt) for the base station
                obj.ultxPkt.txWaveform = txWaveform;
                obj.ultxPkt.txGrid = txGrid;
                obj.ultxPkt.txInfo = txInfo;
                obj.ultxPkt.txChInfo = txChInfo;
                obj.ultxPkt.macPDU = macPDU;
                obj.ultxPkt.phyTB = phyTB;
                obj.ultxPkt.puschIndices = puschIndices;
                obj.ultxPkt.codedTrBlock = codedTrBlock;
                
                obj.ultxPkt.pktType = obj.nodeType;
                obj.ultxPkt.TBS = obj.ulConfig.TBS;
                obj.ultxPkt.txTS = gClk;
                
                % store information of the harq
                obj.ulConfig.HARQ.harqInfo(harqID).macPDU = macPDU;
                obj.ulConfig.HARQ.harqInfo(harqID).phyTB = phyTB;
                obj.ulConfig.HARQ.harqInfo(harqID).pktNumber = obj.ultxPkt.pktNumber;
                obj.ulConfig.HARQ.harqInfo(harqID).pktSize = obj.ultxPkt.pktSize;
                obj.ulConfig.HARQ.harqInfo(harqID).crtTS = obj.ultxPkt.crtTS;
            else
                % re-transmission of data stored in harqProcesses
                fprintf('HARQ number = %d \t Re-transmission.\n', harqID);
                
                % 2. RLC: re-transmit the same mac pdu
                macPDU = obj.ulConfig.HARQ.harqInfo(harqID).macPDU;
                
                % 3. MAC: re-transmit the same TB
                phyTB = obj.ulConfig.HARQ.harqInfo(harqID).phyTB;
                harqProcesses(harqID).ulschTransportBlk = phyTB;
                
                % 4. PHY: Perform PHY processing (TB, Coding, Modulation,
                % waveform, etc).
                [txWaveform, txGrid, codedTrBlock, puschIndices, txInfo,...
                    txChInfo, obj.ulConfig] = ulPHY_fn(obj.ulConfig, harqProcesses(harqID).ulschTransportBlk);
                
                harqProcesses(harqID).codedTrBlkSize = obj.ulConfig.PUSCH.CodedTrBlkSizes(1);
                
                % Get the HARQ ID sequence from 'enbOut' for HARQ processing
                harqProcessSequence = obj.ulConfig.PUSCH.HARQProcessSequence;
                obj.ulConfig.HARQ.harqProcessSequence = harqProcessSequence;
                obj.ulConfig.HARQ.harqProcesses = harqProcesses;
                
                % transmitted packet (ultxPkt) for the base station
                obj.ultxPkt.txWaveform = txWaveform;
                obj.ultxPkt.txGrid = txGrid;
                obj.ultxPkt.txInfo = txInfo;
                obj.ultxPkt.txChInfo = txChInfo;
                obj.ultxPkt.macPDU = macPDU;
                obj.ultxPkt.phyTB = phyTB;
                obj.ultxPkt.puschIndices = puschIndices;
                obj.ultxPkt.codedTrBlock = codedTrBlock;
                
                % update packet information
                obj.ultxPkt.pktNumber = obj.ulConfig.HARQ.harqInfo(harqID).pktNumber;
                obj.ultxPkt.pktSize = obj.ulConfig.HARQ.harqInfo(harqID).pktSize;
                obj.ultxPkt.crtTS = obj.ulConfig.HARQ.harqInfo(harqID).crtTS;
                obj.ultxPkt.pktType = obj.nodeType;
                obj.ultxPkt.TBS = obj.ulConfig.TBS;
                obj.ultxPkt.txTS = gClk;
            end
            
            
        end
        % *************************************************************** %
        
        %% ul RLC
        function [macPDU] = ulRLC_fn(obj, TBS)
            
            ipPkt = obj.txBuffer(end);
            obj.ultxPkt.crtTS = ipPkt.crtTS;
            obj.ultxPkt.pktNumber = ipPkt.pktNumber;
            ipPktSize = ipPkt.pktSize;
            
            if(ipPktSize > TBS)
                obj.ultxPkt.pktSize = TBS;
                obj.ultxPkt.remSize = ipPktSize - TBS;
                obj.txBuffer(end).remSize = ipPktSize - TBS;
            else
                obj.ultxPkt.pktSize = ipPktSize;
                obj.ultxPkt.remSize = 0;
                obj.txBuffer(end) = [];
            end
            
            macPDU = randi([0, 1], TBS, 1);
        end
        % *************************************************************** %
        
        %% perform sl rlc and mac functions (every pscch period)
        function slRLCMAC_fn(obj)
            % 1. RLC: Generate the macPDU                    
            [macPDU, pduLength] = obj.sltxRLC_fn(obj.sltxConfig.TBS);

            % 2. MAC: Generate the phyTB
            [phyTB, slmacPadding, slmacHeader] = obj.sltxMAC_fn(macPDU, pduLength, obj.sltxConfig.TBS);

            obj.MACheaderSL = slmacHeader;
            obj.macPadding = slmacPadding;
            obj.sltxConfig.phyTB = phyTB;
            obj.sltxConfig.macPDU = macPDU;       
        end
        % *************************************************************** %
        
        %% sidelink transmission function
        function slPhypssch_fn(obj)            
            % perform Phy processing and apply fading channel
            phyTB = {obj.sltxConfig.phyTB};
            sciMessage = obj.sltxConfig.sciMessage;
            ue = obj.sltxConfig.ue;
            
            [txWaveform, txChInfo, txInfo] = obj.sltxPHY_fn(phyTB, obj.sltxConfig.NPRB, sciMessage, ue);
            
            % transmission information
            obj.sltxConfig.txWaveform = txWaveform;
            obj.sltxConfig.txChInfo = txChInfo;
            obj.sltxConfig.txInfo = txInfo;
        end
        % *************************************************************** %
        
        %% sidelink pscch transmission function
        function slPhypscch_fn(obj)
            global gClk
            global ulChannel
            global perfectChannel
            
            % perform Phy processing and apply fading channel
            sciMessage = obj.sltxConfig.sciMessage;
            ue = obj.sltxConfig.ue;
            
            sfNumber = mod(gClk, 40);

            [RIVsetTemp,range] = obj.pscchPeriod.getAllowedRIV(sciMessage);
            RIVset = RIVsetTemp(range(:,1) == obj.sltxConfig.NPRB);
            
            % Update PSCCH obj.pscchPeriod number
            obj.pscchPeriod.Config.NPSCCHPeriod = 0; % p = 1;
            
            % Choose a random PSSCH resource
            if(length(RIVset) < 2)
                sciMessage.Allocation.RIV = RIVset(1);
            else
                sciMessage.Allocation.RIV = RIVset(2); %randi(length(RIVset)));
            end
            
            % save sci message 
            obj.sltxConfig.sciMessage = sciMessage;
            
            % % Generate PSCCH obj.pscchPeriod waveform
            txWaveform = pscch_fn(obj.nodeID, obj.pscchPeriod, sciMessage, sfNumber);
%             txWaveform = obj.pscchPeriod.genWaveform(sciMessage, phyTB, sfNumber);
            
            % Set the sampling rate and init time for frequency-selective fading
            ue.CyclicPrefixSL = obj.pscchPeriod.Config.sc_CP_Len_r12;
            
            % SC-FDMA transmission information
            txInfo = lteSLSCFDMAInfo(ue);
            
            % Apply channel
            ulChannel.slchcfg.InitTime = gClk/1000;
            ulChannel.slchcfg.SamplingRate  = txInfo.SamplingRate;
            
            if(~perfectChannel)
                [txWaveform, txChInfo] = ulChannel.lteCh_fn(txWaveform, 'sideLink');
            else
                txChInfo = [];
            end
            
            % transmission information
            obj.sltxConfig.txWaveform = txWaveform;
            obj.sltxConfig.txChInfo = txChInfo;
            obj.sltxConfig.txInfo = txInfo;
        end
        % *************************************************************** %
        
        %% Create the sci message
        function [sciMessage, TBS, ue] = sci_fn(obj)
            global simParameters
            
            % Resource Pool Configuration
            % Configure the PSCCH obj.pscchPeriod parameters for reference pool #1-FDD, 5MHz
            obj.pscchPeriod.Config = PSCCHPeriod.defaultConfig(1, strcat(num2str(simParameters.bandwidth), 'MHz'));
            
            % Create Sidelink Control Information (SCI) Message
            % TBS, modscheme and rv
            TBS = obj.updateTBS_fn();
            
            % Create UE configuration and SCI message with frequency hopping enabled
            ue.NSLRB = obj.pscchPeriod.Config.NSLRB;
            sciMessage.FreqHopping = 0;
            sciMessage.SCIFormat = 'Format0';
            sciMessage = lteSCI(ue,sciMessage);
            
            % Set the hopping bits for Type 2 hopping and the time resource pattern for
            % RMC CC.3 FDD, TS 36.101 Table A.6.4-1
            sciMessage.Allocation.HoppingBits = 1;
            sciMessage.TimeResourcePattern = 106;   % TS 136.213 (page 301) Table 14.1.1.1.1-1: Time Resource pattern Index mapping for NTRP = 8
            
            % Set the PSSCH MCS according to the TBS
            ITBSs = lteMCS(0:28,'PUSCH');
            sciMessage.ModCoding = find(lteTBS(obj.sltxConfig.NPRB,ITBSs)==TBS,1,'first') - 1;
            
            % Set NSAID value, would be assigned by higher layers
            sciMessage.NSAID = 117;
            
            % Store NSAID value in 'expectedNSAID' so that the receiver won't have to
            % access 'sciMessage' at all
            obj.expectedNSAID = sciMessage.NSAID;
            
            % RBs used for sci message 
            sciMessage.PSCCHResource = randi([0 obj.pscchPeriod.NumPSCCHResource-1]);
        end
        % *************************************************************** %
        
        %% update tbs, modulation scheme, and rv
        function [tbs, modScheme, rv] = updateTBS_fn(obj)
            nprb = obj.sltxConfig.NPRB;
            
            imcsTable = [-1 0 0 2 4 6 8 11 13 16 18 21 23 25 27 27];
            if(isempty(obj.sltxConfig.cqiBuffer))
                imcs = imcsTable(4);
            else
                imcs = imcsTable(obj.sltxConfig.cqiBuffer+1);
            end
            
            [itbs, modScheme, rv] = lteMCS(imcs);
            tbs = lteTBS(nprb, itbs);
        end
        % *************************************************************** %
        
        %% sidelink (sl) RLC processing
        function [macPDU, pduSize] = sltxRLC_fn(obj, TBS)
                     
            ipheaderlength = 5+24*3;
            macheaderlength = 32;
            pduSize = TBS - macheaderlength - ipheaderlength; 
            
            [macPDU, obj.RLCtxEntity] = obj.RLCtxEntity.txPDU(pduSize);
            pduSize = length(macPDU); % 772;
                        
            obj.RLCtxEntity.PDUsize = pduSize;
        end
        % *************************************************************** %
        
        %% sidelink MAC processing
        function [phyTB, padding, slmacHeader] = sltxMAC_fn(obj, macPDU, pduLength, TBS)
            
            srcID = obj.nodeID;
            destID = 2;
            LCID = [0 0 0 0 0];
            CE = false;
            
            % Generate the MAC header
            slmacHeader = buildMACheaderSL(srcID, destID, pduLength, LCID, CE);
            
            % Combine header and payload, add padding if needed
            HeaderPayload = [slmacHeader macPDU];
            
            if length(HeaderPayload) < TBS
                padding = zeros(1,TBS - length(HeaderPayload));
                phyTB = [HeaderPayload padding];
            elseif length(HeaderPayload) > TBS
                fprintf("error: MAC PDU is larger than transport block size.\n")
            else
                padding = 0;
                phyTB = HeaderPayload;
            end
        end
        % *************************************************************** %
        
        %% sidelink PHY processing
        function [txWaveform, txChInfo, txInfo] = sltxPHY_fn(obj, phyTB, NPRB, sciMessage, ue)
            
            global ulChannel
            global gClk
            global perfectChannel
            
            sfNumber = mod(gClk, 40);

            [RIVset,range] = obj.pscchPeriod.getAllowedRIV(sciMessage);
            RIVset = RIVset(range(:,1) == NPRB);
            
            p = 1;
            
            % Update PSCCH obj.pscchPeriod number
            obj.pscchPeriod.Config.NPSCCHPeriod = p - 1;
            
            % Choose a random PSSCH resource
            sciMessage.Allocation.RIV = RIVset(2);%randi(length(RIVset)));
            
            % % Generate PSCCH obj.pscchPeriod waveform
            txWaveform = obj.pscchPeriod.genWaveform(obj.nodeID, sciMessage, phyTB, sfNumber);
            
            % Set the sampling rate and init time for frequency-selective fading
            ue.CyclicPrefixSL = obj.pscchPeriod.Config.sc_CP_Len_r12;
            
            % SC-FDMA transmission information
            txInfo = lteSLSCFDMAInfo(ue);
            
            % Apply fading channel
            subframeNo = gClk;
            ulChannel.slchcfg.InitTime = subframeNo/1000;
            ulChannel.slchcfg.SamplingRate  = txInfo.SamplingRate;
            
            if(~perfectChannel)
                [txWaveform, txChInfo] = ulChannel.lteCh_fn(txWaveform, 'sideLink');
            else
                txChInfo = [];
            end
        end
        % *************************************************************** %
        
        %% sidelink receiving function
        function slrx_fn(obj)
            global gClk
            
            % Call the physical layer to receive TBs. 
            % If there is a received TB, call the MAC and RLC layers (This
            % is performed inside the slrxPHY_fn function)
            if(strcmp(obj.commDir, 'rx') || strcmp(obj.commDir, 'txrx'))  
                % Clear rx parameters at the pscch period boundary
                if(mod(gClk, 40) == 0)        
                    % reset tx parameters at the PSCCH period
                    % boundary
                    for i = 1:length(obj.slrxConfig)
                        obj.slrxConfig{i}.sciDecoded = 0;
                        obj.slrxConfig{i}.sciMessageRx = [];
                    end
                end
                
                obj.slrxPHY_fn();    
            end
        end
        % *************************************************************** %
        
        %% sidelink receiving function
        function slrxPHY_fn(obj)
            global logData
            global gClk
            global ulChannel
            global perfectChannel
            
            % Steps:
            % 1. Apply large-scale fading, calculate sinr and cqi
            % 2. Perform sl rx processing (phy, mac, rlc) to compute the throughput
            % 3. Perform sl rx MAC and RLC processing
            % 4. Store the results of ul tx/rx (BER, BLER, Delay,
            % Throughput, etc)
            
            for u = 1:length(obj.txUsers)
                udTemp = obj.txUsers(u);
                udnodeID = udTemp.absID;
                
                if(~isempty(udTemp.sltxConfig.txWaveform))
                    
                    % 1. Apply large-scale fading, calculate sinr and cqi
                    if(perfectChannel)
                        rxWaveform = udTemp.sltxConfig.txWaveform;
                        sinrdB = 50;
                        cqi = 9;
                    else
                        [rxWaveform, sinrdB, cqi] = obj.largeFading_fn(udTemp);
                    end
                    
                    udTemp.sltxConfig.sinrBuffer = sinrdB;
                    udTemp.sltxConfig.cqiBuffer = cqi;
                    
                    % 2. Perform sl rx phy processing
                    [rxTB, slschsDecoded] = obj.slrxProc_fn(rxWaveform, udTemp);
                    
                    % 3. Perform sl rx MAC and RLC processing
                    if(~isempty(rxTB))
                        rxpktSize = 0;
                        for i = 1:length(rxTB)
                            if slschsDecoded(i)
                                rxrlcPDU = rxTB{i}(length(udTemp.MACheaderSL)+1:udTemp.sltxConfig.TBS-length(udTemp.macPadding));
                                obj.RLCrxEntity{udnodeID} = obj.RLCrxEntity{udnodeID}.rxPDU(rxrlcPDU');
                                rxpktSize = rxpktSize + length(rxrlcPDU);
                            end
                        end
                        
                        % Put SDUs in a receiver buffer
                        if(~isempty(obj.RLCrxEntity{udnodeID}.SDUs))
                            sdusTemp = obj.RLCrxEntity{udnodeID}.SDUs;
                            
                            for sduIdx = 1:length(sdusTemp)
                                sdu = {sdusTemp{sduIdx}};
                                sduHeader = readMetadata(sdu);
                                obj.rxBuffer = [obj.rxBuffer; sdusTemp{sduIdx}];
                                pktSize = length(sdusTemp{sduIdx});
                                
                                pktInfo.txType = udTemp.nodeType;
                                pktInfo.txID = sduHeader.nodeID;
                                pktInfo.rxType = obj.nodeType;
                                pktInfo.rxID = obj.nodeID;
                                pktInfo.pktNumber = sduHeader.packetNum;
                                pktInfo.pktSize = pktSize;
                                pktInfo.crtTS = sduHeader.TOC;
                                pktInfo.txTS = sduHeader.TOTX;
                                pktInfo.rxTS = gClk + ulChannel.propDelay_fn(obj.disttoNode_fn(udTemp));
                                pktInfo.delay = pktInfo.rxTS - pktInfo.crtTS;
                                pktInfo.bitTp = pktSize / ((pktInfo.rxTS - pktInfo.crtTS)*1e-3);     
                                pktInfo.pktErr = 0;

                                logData.saveLog_fn(pktInfo);
                            end
                            
                            udTemp.sltxConfig.recentDelay = pktInfo.rxTS - pktInfo.crtTS;
                            obj.RLCrxEntity{udnodeID}.SDUs = [];  
                            
                        elseif(obj.RLCrxEntity{udnodeID}.SDUsDropped > 0)
                            pktInfo.txType = udTemp.nodeType;
                            pktInfo.txID = udTemp.nodeID;
                            pktInfo.rxType = obj.nodeType;
                            pktInfo.rxID = obj.nodeID;
                            pktInfo.pktNumber = 0;
                            pktInfo.pktSize = 0;
                            pktInfo.crtTS = 0;
                            pktInfo.txTS = 0;
                            pktInfo.rxTS = 0;
                            pktInfo.delay = 0;
                            pktInfo.bitTp = 0;     
                            pktInfo.pktErr = obj.RLCrxEntity{udnodeID}.SDUsDropped;

                            logData.saveLog_fn(pktInfo);
                        end           
                    end
                end
                
            end
            
        end
        % *************************************************************** %
        
        %% ApplyChannel: Apply pathloss and shadowing to both main link and interferers
        function [waveform, sinrdB, cqi] = largeFading_fn(obj, udTemp)
            
            global ulChannel
            global simParameters
            
            % Recieved waveform
            txWaveform = udTemp.sltxConfig.txWaveform;
            
            [rxWaveform] = ulChannel.pathLoss_fn(txWaveform, udTemp.disttoRxer(obj.absID));
%             rxWaveform = txWaveform;
            
            % Interferers
            udsInterferers = obj.udsInterfering;
            intrxWaveform = zeros(size(rxWaveform));
            
            for i = 1 : length(udsInterferers)
                inttxWaveform = udsInterferers(i).sltxConfig.txWaveform;
                intDistance = disttoNode_fn(obj, udsInterferers(i));
                
                [waveformTemp]  = ulChannel.pathLoss_fn(inttxWaveform, intDistance);
                
                if(~isempty(waveformTemp))
                    intrxWaveform = intrxWaveform + waveformTemp;
                end
            end
            
            % White noise
            No_dB = -174 + 10 * log10(simParameters.bandwidth * 1e6) - 30;
            No = 10^(No_dB/10);
            noise = No * complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
            
            waveform = rxWaveform + intrxWaveform + noise;
            
            % SINR estimation
            sinrdB = snr(rxWaveform, intrxWaveform+noise);
            
            % CQI: map sinr wideband to cqi value (wideband)
            cqi = snrTocqi(sinrdB);
        end
        % *************************************************************** %
        
        %% sl rx processing (phy, mac, rlc)
        function [trbRxs, slschsDecoded] = slrxProc_fn(obj, rxWaveform, udTemp)
            global gClk
            
            sfNumber = mod(gClk, 40);
                        
            % init
            trbRxs = [];
            slschsDecoded = [];
            udnodeID = udTemp.absID;
            
            % Create PSCCH receiver configuration
            pscch.SidelinkMode = 'D2D';
            pscch.NSLRB = udTemp.pscchPeriod.Config.NSLRB;
            pscch.CyclicPrefixSL = udTemp.pscchPeriod.Config.sc_CP_Len_r12;

            % Get the number of current frame if it is PSCCH frame
            nsfpscch = udTemp.pscchPeriod.PSCCHSubframePool(sfNumber == udTemp.pscchPeriod.PSCCHSubframePool);
            
            % Decode SCI or traffic depending on the current subframe
            if(~isempty(nsfpscch) && ~obj.slrxConfig{udnodeID}.sciDecoded)
                % If this is a PSCCH subframe, decode the control message.
                [obj.slrxConfig{udnodeID}.sciMessageRx, obj.slrxConfig{udnodeID}.sciDecoded] =...
                    decodePSCCHsf(rxWaveform, udTemp, pscch, sfNumber);
                
                fprintf("User: %d \t Receiving PSCCH subframe = %d \t Decoding error = %d\n",...
                    obj.nodeID, sfNumber+1, ~obj.slrxConfig{udnodeID}.sciDecoded);
                
            elseif(~isempty(obj.slrxConfig{udnodeID}.sciMessageRx))
                [sf, prb] = udTemp.pscchPeriod.getPSSCHResources(obj.slrxConfig{udnodeID}.sciMessageRx);

                if(any(sfNumber == sf)) 

                    nsfpssch = find(sfNumber == sf);

                    if obj.slrxConfig{udnodeID}.sciDecoded
                        % Create PSSCH receiver configuration
                        pssch.SidelinkMode = 'D2D';
                        pssch.NSLRB = udTemp.pscchPeriod.Config.NSLRB;
                        pssch.CyclicPrefixSL = udTemp.pscchPeriod.Config.data_CP_Len_r12;
                        pssch.NTurboDecIts = 5;

                        % Remove any resources that overlap with synchronization
                        % subframes. PSSCH will not be transmitted in these subframes
                        [~,syncIndex] = intersect(sf, udTemp.pscchPeriod.SyncSubframes);
                        prb(:,syncIndex) = [];

                        % Configure the PSSCH receiver NSAID value and modulation from
                        % the decoded SCI message, and establish the SL-SCH transport
                        % block size (TBS)
                        pssch.NSAID = obj.slrxConfig{udnodeID}.sciMessageRx.NSAID;
                        [ITBS,modulation] = lteMCS(obj.slrxConfig{udnodeID}.sciMessageRx.ModCoding,'PUSCH');
                        if (strcmpi(modulation,'64QAM'))
                            modulation = '16QAM';
                        end
                        pssch.Modulation = modulation;
                        NPRB = max(size(prb,1),1);
                        TBS = lteTBS(NPRB,ITBS);

                        % Repeat for each PSSCH transmission instance until the SL-SCH
                        % is decoded
                        rvsequence = [0 2 3 1];
                        trbRxs = {};

                        % Configure the PSSCH receiver for the PSSCH subframe
                        % number, PRB allocation and redundancy version (RV). If
                        % the RV index is zero (corresponding to the first PSSCH
                        % transmission instance in a block of 4 transmissions),
                        % reset the receiver buffer
                        pssch.NSubframePSSCH = nsfpssch-1;
                        pssch.PRBSet = prb(:, nsfpssch);       
                        pssch.RV = rvsequence(mod(nsfpssch-1, 4)+1);                            

                        % Perform timing synchronization, extract the appropriate
                        % subframe of the received waveform, and perform SC-FDMA
                        % demodulation
                        offset = lteSLFrameOffsetPSSCH(pssch, rxWaveform);

                        % Assume perfect synchronization
                        offset = 0;
                        subframeWaveform = rxWaveform(offset+1:end, :);
                        subframe = lteSLSCFDMADemodulate(pssch, subframeWaveform);

                        % Perform channel estimation, extract the received PSSCH
                        % symbols and the corresponding channel estimate, and
                        % perform equalization
                        cec.PilotAverage = 'UserDefined';
                        cec.TimeWindow = 15;
                        cec.FreqWindow = 23;
                        cec.InterpType = 'linear';

                        [hest, nest] = lteSLChannelEstimatePSSCH(pssch, cec, subframe);
                        psschIndices = ltePSSCHIndices(pssch);
                        [psschRx, psschHest] = lteExtractResources(psschIndices, subframe, hest);
                        psschSymbols = lteEqualizeMMSE(psschRx, psschHest, nest);

                        % Demodulate the PSSCH
                        codedSlschBits = ltePSSCHDecode(pssch, psschSymbols);

                        % Decode the SL-SCH including soft combining into the
                        % receiver buffer and check the CRC
                        [trbRx, slschCRC, ~] = lteSLSCHDecode(pssch, TBS, codedSlschBits,[]);
                        
                        fprintf("User: %d \t Receiving PSSCH subframe = %d \t Decoding error = %d\n", obj.nodeID, sfNumber+1, slschCRC);
                        if (slschCRC==0)
                            trbRxs = [trbRxs; trbRx];
                            slschsDecoded = true;
                        else
                            slschsDecoded = false;
                        end
                    else
                        trbRxs = [];
                        slschsDecoded = 0;
                    end
                end
            end
        end
        % *************************************************************** %
        
        %% sl rx processing (phy, mac, rlc)
        %         function [trbRxs, delayRXs, slschsDecoded] = slrxProc_fn2(obj, rxWaveform, udTemp)
        %
        %             TBs = udTemp.sltxPkt.phyTB;
        %
        %             % Create PSCCH receiver configuration
        %             pscch.SidelinkMode = 'D2D';
        %             pscch.NSLRB = udTemp.pscchPeriod.Config.NSLRB;
        %             pscch.CyclicPrefixSL = udTemp.pscchPeriod.Config.sc_CP_Len_r12;
        %
        %             % Establish the number of time-domain samples 'Nt' per subframe and
        %             % number of SC-FDMA symbols 'L' per subframe
        %             pscchGridInfo = lteSLSCFDMAInfo(pscch);
        %             Nt = pscchGridInfo.SamplingRate * 1e-3;
        %             L = numel(pscchGridInfo.CyclicPrefixLengths);
        %
        %             % If SCI BLER measurement is configured:
        %             sciDecoded = false;
        %
        %             % Synchronize and extract each subframe of the PSCCH subframe
        %             % pool from the received waveform and perform SC-FDMA
        %             % demodulation. The timing offset used is the offset
        %             % corresponding to the strongest correlation peak when
        %             % correlating with the PSCCH DRS for each PRB in the PSCCH
        %             % subframe pool
        %             nSubframes = udTemp.pscchPeriod.Config.sc_Period_r12;
        %             pscchPoolGrid = repmat(lteSLResourceGrid(pscch),1,nSubframes,size(rxWaveform,2));
        %             for l = udTemp.pscchPeriod.PSCCHSubframePool
        %                 subframeWaveform = rxWaveform(l*Nt + (1:Nt),:);
        %                 bestCorr = 0.0;
        %                 for prb = udTemp.pscchPeriod.PSCCHResourceBlockPool.'
        %                     pscch.PRBSet = prb;
        %                     [thisoffset,thiscorr] = lteSLFrameOffsetPSCCH(pscch,subframeWaveform);
        %                     % Find the best correlation across all antennas and
        %                     % store the corresponding offset
        %                     maxCorr = max(thiscorr(:));
        %                     if maxCorr > bestCorr
        %                         bestCorr = maxCorr;
        %                         offset = thisoffset;
        %                     end
        %                 end
        %                 subframe = lteSLSCFDMADemodulate(pscch,rxWaveform(l*Nt + offset + (1:Nt),:));
        %                 pscchPoolGrid(:,l*L + (1:L),:) = subframe;
        %             end
        %
        %             % Repeat for each PSCCH resource until the SCI is decoded
        %             pscchResource = 0;
        %             while (pscchResource < udTemp.pscchPeriod.NumPSCCHResource && ~sciDecoded)
        %
        %                 % Get the PSCCH subframes 'sf' and PRB allocations 'prb'
        %                 % for the current PSCCH resource
        %                 [sf,prb] = udTemp.pscchPeriod.getPSCCHResources(pscchResource);
        %
        %                 % Remove any resources that overlap with synchronization
        %                 % subframes. PSCCH will not be transmitted in these
        %                 % subframes
        %                 [~,syncIndex] = intersect(sf,udTemp.pscchPeriod.SyncSubframes);
        %                 sf(syncIndex) = [];
        %                 prb(syncIndex) = [];
        %
        %                 % Repeat for each PSCCH transmission instance until the SCI
        %                 % is decoded
        %                 tx = 1;
        %                 while (tx <= min(length(sf),2) && ~sciDecoded)
        %
        %                     % Configure the PSCCH receiver for the PRB allocation
        %                     pscch.PRBSet = prb(:,tx);
        %
        %                     % Select the appropriate subframe resource grid from
        %                     % the PSCCH subframe pool
        %                     subframe = pscchPoolGrid(:,sf(tx)*L + (1:L),:);
        %
        %                     % Perform channel estimation, extract the received
        %                     % PSCCH symbols and the corresponding channel estimate,
        %                     % and perform equalization
        %                     cec.PilotAverage = 'UserDefined';
        %                     cec.TimeWindow = 15;
        %                     cec.FreqWindow = 23;
        %                     cec.InterpType = 'linear';
        %
        %                     [hest,nest] = lteSLChannelEstimatePSCCH(pscch,cec,subframe);
        %                     [pscchIndices,pscchIndicesInfo] = ltePSCCHIndices(pscch);
        %                     [pscchRx,pscchHest] = lteExtractResources(pscchIndices,subframe,hest);
        %                     pscchSymbols = lteEqualizeMMSE(pscchRx,pscchHest,nest);
        %
        %                     % If we are receiving the first PSCCH transmission
        %                     % instance, reset the receiver buffer
        %                     if (tx==1)
        %                         codedSciBits = zeros(pscchIndicesInfo.G,1);
        %                     end
        %
        %                     % Demodulate the PSCCH and add the result into the
        %                     % receiver buffer
        %                     codedSciBits = codedSciBits + ltePSCCHDecode(pscchSymbols);
        %
        %                     % Decode the SCI message. If successful (CRC=0), check
        %                     % the NSAID field of the decoded message against the
        %                     % transmitted message
        %                     sciInfo = lteSCIInfo(pscch);
        %                     [sciBits,sciCRC] = lteSCIDecode(sciInfo.Format0,codedSciBits);
        %                     if (sciCRC==0)
        %                         sciMessageRx = lteSCI(pscch,sciBits);
        %                         if (sciMessageRx.NSAID == udTemp.expectedNSAID)
        %                             sciDecoded = true;
        %                         end
        %                     end
        %
        %                     % Increment the PSCCH transmission instance index
        %                     tx = tx + 1;
        %
        %                 end
        %
        %                 % Increment the PSCCH resource number
        %                 pscchResource = pscchResource + 1;
        %
        %             end
        %
        %             % If SCI decoding failed and the SCI is assumed for SL-SCH
        %             % decoding, set the decoded SCI message equal to the transmitted
        %             % SCI message
        %             if (~sciDecoded)
        %                 fprintf('   SCI decoding failed\n');
        %             end
        %
        %             % If SL-SCH BLER measurement is configured and if the SCI was
        %             % successfully decoded or if SCI decoding success is assumed,
        %             % perform PSSCH reception and SL-SCH decoding
        %             slschDecodedFinal = false;
        %             if sciDecoded
        %
        %                 % Create PSSCH receiver configuration
        %                 pssch.SidelinkMode = 'D2D';
        %                 pssch.NSLRB = udTemp.pscchPeriod.Config.NSLRB;
        %                 pssch.CyclicPrefixSL = udTemp.pscchPeriod.Config.data_CP_Len_r12;
        %                 pssch.NTurboDecIts = 5;
        %
        %                 % Establish the number of time-domain samples 'Nt' per subframe
        %                 psschGridInfo = lteSLSCFDMAInfo(pssch);
        %                 Nt = psschGridInfo.SamplingRate * 1e-3;
        %
        %                 % Get the PSSCH subframes 'sf', PRB allocations 'prb' and PSSCH
        %                 % subframe numbers 'nsf' for the configuration given by the
        %                 % decoded SCI message. 'sf' are the subframe numbers within the
        %                 % SC udTemp.pscchPeriod, 'nsf' are the subframe numbers within the PSSCH
        %                 % subframe pool
        %                 [sf,prb,nsf] = udTemp.pscchPeriod.getPSSCHResources(sciMessageRx);
        %
        %                 % Create the redundancy version index sequence (IRV)
        %                 % corresponding to the PSSCH transmission instances. This must
        %                 % be computed "up front" prior to the removal of any
        %                 % synchronization subframes in the next step, otherwise the RV
        %                 % sequence is lost
        %                 IRV = 0:length(sf)-1;
        %
        %                 % Remove any resources that overlap with synchronization
        %                 % subframes. PSSCH will not be transmitted in these subframes
        %                 [~,syncIndex] = intersect(sf,udTemp.pscchPeriod.SyncSubframes);
        %                 sf(syncIndex) = [];
        %                 prb(:,syncIndex) = [];
        %                 nsf(syncIndex) = [];
        %                 IRV(syncIndex) = [];
        %
        %                 % Configure the PSSCH receiver NSAID value and modulation from
        %                 % the decoded SCI message, and establish the SL-SCH transport
        %                 % block size (TBS)
        %                 pssch.NSAID = sciMessageRx.NSAID;
        %                 [ITBS,modulation] = lteMCS(sciMessageRx.ModCoding,'PUSCH');
        %                 if (strcmpi(modulation,'64QAM'))
        %                     modulation = '16QAM';
        %                 end
        %                 pssch.Modulation = modulation;
        %                 NPRB = max(size(prb,1),1);
        %                 TBS = lteTBS(NPRB,ITBS);
        %
        %                 % Repeat for each PSSCH transmission instance until the SL-SCH
        %                 % is decoded
        %                 tx = 1;
        %                 rvsequence = [0 2 3 1];
        %                 trbRxs = cell(1,length(TBs));
        %                 delayRXs = zeros(1,length(TBs));
        %                 slschsDecoded = zeros(1,length(TBs));
        %                 k = 1;
        %                 while (tx <= length(TBs)*4 && ~slschDecodedFinal)
        %
        %                     % Configure the PSSCH receiver for the PSSCH subframe
        %                     % number, PRB allocation and redundancy version (RV). If
        %                     % the RV index is zero (corresponding to the first PSSCH
        %                     % transmission instance in a block of 4 transmissions),
        %                     % reset the receiver buffer
        %                     pssch.NSubframePSSCH = nsf(tx);
        %                     pssch.PRBSet = prb(:,tx);
        %                     pssch.RV = rvsequence(mod(IRV(tx),4)+1);
        %
        %                     % Perform timing synchronization, extract the appropriate
        %                     % subframe of the received waveform, and perform SC-FDMA
        %                     % demodulation
        %                     offset = lteSLFrameOffsetPSSCH(pssch,rxWaveform(sf(tx)*Nt + (1:Nt),:));
        %                     subframeWaveform = rxWaveform(sf(tx)*Nt + offset + (1:Nt),:);
        %                     subframe = lteSLSCFDMADemodulate(pssch,subframeWaveform);
        %
        %                     % Perform channel estimation, extract the received PSSCH
        %                     % symbols and the corresponding channel estimate, and
        %                     % perform equalization
        %                     [hest,nest] = lteSLChannelEstimatePSSCH(pssch,cec,subframe);
        %                     psschIndices = ltePSSCHIndices(pssch);
        %                     [psschRx,psschHest] = lteExtractResources(psschIndices,subframe,hest);
        %                     psschSymbols = lteEqualizeMMSE(psschRx,psschHest,nest);
        %
        %                     % Demodulate the PSSCH
        %                     codedSlschBits = ltePSSCHDecode(pssch,psschSymbols);
        %
        %                     % Decode the SL-SCH including soft combining into the
        %                     % receiver buffer and check the CRC
        %                     [trbRx,slschCRC,~] = lteSLSCHDecode(pssch,TBS,codedSlschBits,[]);
        %                     if (slschCRC==0)
        %                         %Save TB k in kth index of trbRxs
        %                         trbRxs{k} = trbRx;
        %                         %sf(tx) is number of subframes into PSCCH udTemp.pscchPeriod that TB k
        %                         %was recieved
        %                         delayRXs(k) = sf(tx); %may need to refine this delay estimate; should be a good guess for now though
        %                         %kth TB was successfully decoded
        %                         slschsDecoded(k) = true;
        %                         k = k+1;
        %                         if tx/4 > k-1
        %                             %decoded the last TB in the PSCCH Period
        %                             slschDecodedFinal = true;
        %                         else
        %                             %set tx to index of first RV of next TB
        %                             tx = tx + 5 - mod(tx,4);
        %                         end
        %                     else
        %                         % Increment the PSSCH transmission instance index
        %                         tx = tx + 1;
        %                         if mod(tx,4) == 1
        %                             %since this is after incrementation of tx, this will
        %                             %only be true if the current TB failed to decode for
        %                             %this TTI, and the next subframe is the start of a new
        %                             %TB RV sequence
        %                             trbRxs{k} = [];
        %                             slschsDecoded(k) = false;
        %                             k = k+1;
        %                         end
        %                     end
        %                 end
        %             else
        %                 trbRxs = [];
        %                 delayRXs = 0;
        %                 slschsDecoded = 0;
        %             end
        %         end
        % *************************************************************** %
        
        %% Distance between the object and a specified node
        function distance = disttoNode_fn(udTemp, Node)
            distance = sqrt(sum((udTemp.nodePosition - Node.nodePosition).^2));
        end
        % *************************************************************** %
        
        %% D2D distance to transmitter
        function disttoTxer_fn(udTemp, txTemp)
            
            for u = 1:length(txTemp)
                distTemp = txTemp(u).nodePosition;
                udTemp.disttoTxer(u) = sqrt(sum((udTemp.nodePosition - distTemp).^2));
            end
        end
        % *************************************************************** %

        %% D2D distance to receiver
        function disttoRxer_fn(udTemp, txTemp)
            udTemp.disttoRxer = zeros(1, length(txTemp));
            for u = 1:length(txTemp)
                distTemp = txTemp(u).nodePosition;
                udTemp.disttoRxer(u) = sqrt(sum((udTemp.nodePosition - distTemp).^2));
            end
        end
        % *************************************************************** %
    end
end






