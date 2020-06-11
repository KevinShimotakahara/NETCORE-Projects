
% *********************************************************************** %
% Dscription: Generic Node class
% *********************************************************************** %

classdef eNodeB < handle
    properties
        % Node identification
        nodeID              % Node ID
        nodeType
        cellRadius
        nodePosition        % (x, y) points on the grid
        
        udsAttached
        nudsAttached
        udsInterfering
        
        enbScheduler
        sendPools
        ulPool
        slPool
    end
    
    methods
        %% Constructor
        function obj = eNodeB
            obj.nodeID = [];
            obj.nodeType = [];
            obj.cellRadius = [];
            obj.nodePosition = [];
            
            obj.udsAttached = [];
            obj.nudsAttached = [];
            obj.udsInterfering  = [];
            
            obj.enbScheduler = [];
            obj.sendPools = 0;
            obj.ulPool = [];
            obj.slPool = [];
        end
        % *************************************************************** %
        
        %% Scheduler configuration
        function setScheduler_fn(obj, SchType)
            % Select the scheduler
            switch(SchType)
                case 'RoundRobin'
                    obj.enbScheduler = Schedulers.RoundRobin;
                    
                case 'MaxThroughput'
                    obj.enbScheduler = Schedulers.MaxThroughput;
                    
                case 'PropFairness'
                    obj.enbScheduler = Schedulers.PropFairness;
                    
                case 'QLearning'
                    obj.enbScheduler = Schedulers.QLearning;
                    
                case 'Fixed'
                    obj.enbScheduler = Schedulers.Fixed;
                    
                otherwise   % Round Robin is the default scheduler
                    obj.enbScheduler = Schedulers.RoundRobin;
            end
        end
        % *************************************************************** %
        
        %% Attach users to eNB
        function udsAttach_fn(obj, udTemp)
            if(~ismember(udTemp, obj.udsAttached))
                obj.udsAttached = [obj.udsAttached, udTemp];
                obj.nudsAttached = obj.nudsAttached + 1;     % Number of attached users.
            end
        end
        % *************************************************************** %
        
        %% Positions Allocation: Allocate position of the eNB
        function enbPosalloc_fn(obj)
            
            global simParameters
            
            gridSize  = simParameters.GridSize * 1000;
            gridLimit = gridSize - obj.cellRadius;
            flag = 1;
            
            while(flag == 1)
                xPos = uniformRV(obj.cellRadius, gridLimit);    % unifrom distribution
                yPos = uniformRV(obj.cellRadius, gridLimit);
                
                obj.nodePosition = [xPos, yPos];
                
                if(xPos > gridLimit || yPos > gridLimit || xPos < obj.cellRadius...
                        || yPos < obj.cellRadius)
                    flag = 1;
                else
                    flag = 0;
                end
            end
        end
        % *************************************************************** %
        
        %% Allocate positions to the attached users
        function udsPosalloc_fn(obj)
            for u = 1:length(obj.udsAttached)
                udTemp = obj.udsAttached(u);
                enbRadius = obj.cellRadius;
                
%                 if(strcmp(udTemp.commDir, 'tx') || strcmp(udTemp.commDir, 'txrx'))
%                     udTemp.nodePosition = posAlloc(enbRadius/4, obj.nodePosition);
%                 else
                    udTemp.nodePosition = posAlloc(enbRadius, obj.nodePosition);
%                 end
                
                udTemp.disttoCell = sqrt((obj.nodePosition(1) - udTemp.nodePosition(1))^2 +...
                    (obj.nodePosition(2) - udTemp.nodePosition(2))^2);
            end
        end
        % *************************************************************** %
        
        %% Interfering users (UE of other eNBs, and D2D transmitters
        function udsInterfer_fn(obj)
            global UDs
            
            for u = 1:length(UDs)
                udTemp = UDs(u);
                
                if(ismember(udTemp, obj.udsAttached) || strcmp(udTemp, "CTRLER"))
                    continue;
                else
                    obj.udsInterfering = [obj.udsInterfering, udTemp];
                end
            end
        end
        % *************************************************************** %
        
        %% eNB ul rx from UEs
        function ulrx_fn(obj)
            global logData
            global gClk
            
            % Steps:
            % A loop to check all the attached devices to the eNB.
            % 1. Apply large-scale fading, calculate sinr and cqi.
            % 2. Perform ul rx processing to compute the throughput.
            % 3. Store the results of ul tx/rx (BER, BLER, Delay,
            % Throughput, etc)
            
            for u = 1:length(obj.udsAttached)
                udTemp = obj.udsAttached(u);
                
                switch(udTemp.nodeType)
                    
                    % UE: Check for ul tx
                    case 'UE'
                        if(~isempty(udTemp.ultxPkt))
                            % 1. Apply large-scale fading, calculate sinr
                            % and cqi and update configuration parameters
                            [rxWaveform, ~, sinrdB, cqi] = obj.largeFading_fn(udTemp);
%                             rxWaveform = udTemp.ultxPkt.txWaveform;
                            
                            udTemp.ulConfig.sinrBuffer = sinrdB;
                            udTemp.ulConfig.cqiBuffer = cqi;
                            
                            % 2. Perform ul rx processing to compute the throughput
                            [udTemp.ulConfig, bitTp, pktErr] = obj.ulrxProc_fn(rxWaveform, udTemp.ulConfig);

                            % 3. Store the results of ul tx/rx (BER, BLER, Delay,
                            % Throughput, etc)
                            pktInfo.txType = udTemp.nodeType;
                            pktInfo.txID = sduHeader.nodeID;
                            pktInfo.rxType = obj.nodeType;
                            pktInfo.rxID = obj.nodeID;
                            pktInfo.pktNumber = udTemp.ultxPkt.pktNumber;
                            pktInfo.pktSize = udTemp.ultxPkt.pktSize;
                            pktInfo.crtTS = udTemp.ultxPkt.crtTS;
                            pktInfo.txTS = udTemp.ultxPkt.txTS;
                            pktInfo.rxTS = gClk;
                            pktInfo.bitTp = bitTp;
                            pktInfo.pktErr = pktErr;

                            logData.saveLog_fn(pktInfo);
                            
                            % clear the ultxPkt
                            udTemp.ultxPkt = [];
                            
                            % Print results of current harq
                            fprintf('HARQ number = %d \t blkerr = %d\n', udTemp.ulConfig.HARQ.harqID,...
                                udTemp.ulConfig.HARQ.harqProcesses(udTemp.ulConfig.HARQ.harqID).decState.BLKCRC);
                        end
                end
            end
        end
        % *************************************************************** %
        
        %% ApplyChannel: Apply pathloss and shadowing to both main link and interferers
        function [waveform, rxGrid, sinrdB, cqi] = largeFading_fn(obj, udTemp)
            
            global ulChannel
            global simParameters
            
            % Recieved waveform
            txWaveform = udTemp.ultxPkt.txWaveform;
            txGrid = udTemp.ultxPkt.txGrid;
            
            [rxWaveform, rxGrid] = ulChannel.pathLoss_fn(txWaveform, txGrid, udTemp.disttoCell);
            
            % Interferers
            udsInterferers = obj.udsInterfering;
            intrxWaveform = zeros(size(rxWaveform));
            intrxGrid = zeros(size(txGrid));
            
            for i = 1 : length(udsInterferers)
                switch(udsInterferers(i).nodeType)
                    case 'UE'
                        inttxWaveform = udsInterferers(i).ultxPkt.txWaveform;
                        inttxGrid = udsInterferers(i).ultxPkt.txGrid;
                        
                    case 'D2D'
                        inttxWaveform = udsInterferers(i).sltxPkt.txWaveform;
                        inttxGrid = udsInterferers(i).sltxPkt.txGrid;
                end
                
                intDistance = disttoNode_fn(obj, udsInterferers(i));
                
                [waveformTemp, gridTemp]  = ulChannel.pathLoss_fn(inttxWaveform, inttxGrid, intDistance);
                
                if(~isempty(waveformTemp))
                    intrxWaveform = intrxWaveform + waveformTemp;
                    intrxGrid = intrxGrid + gridTemp;
                end
            end
            
            % White noise
            No_dB = -174 + 10 * log10(simParameters.bandwidth * 1e6) - 30;
            No = 10^(No_dB/10);
            noise = No * complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
            
            waveform = rxWaveform + intrxWaveform + noise;
            
            % SINR estimation
            [sinrWB, sinrSB] = ComputeSINR(rxWaveform, intrxWaveform, 5, ...
                simParameters.bandwidth*1e6, udTemp.ultxPkt.txInfo);    % sinr wideband and subband
            
            sinrdB = [sinrWB; sinrSB];
            
            % CQI: map sinr wideband to cqi value
            cqi = snrTocqi(sinrWB);
        end
        % *************************************************************** %
        
        %% ul processing: Phy, MAC, RLC
        function [ulConfig, bitTp, blkError] = ulrxProc_fn(obj, rxWaveform, ulConfig)
            
            % channel estimation configuration
            cec.PilotAverage = 'UserDefined'; % Type of pilot averaging 
            cec.FreqWindow = 13;              % Frequency averaging windows in REs
            cec.TimeWindow = 1;               % Time averaging windows in REs
            cec.InterpType = 'cubic';         % Interpolation type
            cec.Reference = 'Antennas';       % Reference for channel estimation
            cec.Window = 'Left';
            
            % Calculate synchronization offset
            offsetused = lteULFrameOffset(ulConfig, ulConfig.PUSCH, rxWaveform);
%             if (offset < 25)
%                 offsetused = offset;
%             end
            
            % SC-FDMA demodulation
            rxSubframe = lteSCFDMADemodulate(ulConfig, rxWaveform(1+offsetused:end, :));
            
            if(~isempty(rxSubframe))
                % Channel and noise power spectral density estimation
                [estChannelGrid, noiseEst] = lteULChannelEstimate(ulConfig, ulConfig.PUSCH, cec, rxSubframe);

                % Extract REs corresponding to the PUSCH from the given subframe
                % across all receive antennas and channel estimates
                puschIndices = ltePUSCHIndices(ulConfig, ulConfig.PUSCH);
                [puschRx, puschEstCh] = lteExtractResources(puschIndices, rxSubframe, estChannelGrid);

                % MMSE equalization
                rxSymbols = lteEqualizeMMSE(puschRx, puschEstCh, noiseEst);

                % Update frc.PUSCH to carry complete information of the UL-SCH
                % coding configuration
                harqProcesses = ulConfig.HARQ.harqProcesses;
                harqID = ulConfig.HARQ.harqID;
                ulConfig.PUSCH = lteULSCHInfo(ulConfig, ulConfig.PUSCH, harqProcesses(harqID).trBlkSize, 'chsconcat');

                % Decode the PUSCH
                rxEncodedBits = ltePUSCHDecode(ulConfig, ulConfig.PUSCH, rxSymbols);
                
                % Decode the UL-SCH channel and store the block CRC error for given
                % HARQ process
                [rxDecodedBits, harqProcesses(harqID).crc, harqProcesses(harqID).decState] = lteULSCHDecode(...
                    ulConfig, ulConfig.PUSCH, ulConfig.TBS, ...
                    rxEncodedBits, harqProcesses(harqID).decState);

                ulConfig.HARQ.harqProcesses = harqProcesses;
                
                % bit throughput calculation
                bitTp = harqProcesses(harqID).trBlkSize.*(1-harqProcesses(harqID).crc);
                blkError = harqProcesses(harqID).decState.BLKCRC; 
            else               
                harqProcesses = ulConfig.HARQ.harqProcesses;
                harqID = ulConfig.HARQ.harqID;
                
                harqProcesses(harqID).decState.BLKCRC = 1;
                harqProcesses(harqID).crc = 1;
                
                ulConfig.HARQ.harqProcesses = harqProcesses;
                bitTp = 0;
                blkError = 1;
            end
            
        end
        % *************************************************************** %
        
        %% eNB transmits either pools configuration to (UEs/D2D2) or ul grants to UEs. (Update)
        function dltx_fn(obj)
            global UDs
            
            % UEs ul grant: Perform ul scheduling (We use equal portions
            % for testing only)
            [users, ulGrid] = obj.enbScheduler.schedule_fn();
                        
            for u = users
                udIdx = ([UDs.nodeID] == u);
                udTemp = UDs(udIdx);
                udTemp.ulConfig.PUSCH.PRBSet = find(ulGrid == u) - 1;
                udTemp.ulConfig.SR = 0;
                udTemp.ulConfig.schFlag = 1;
                
                % Other configurations (MCS, HARQ, TBS, etc)
                IMCSTable = [-1 0 0 2 4 6 8 11 13 16 18 21 23 25 27 27];

                if(isempty(udTemp.cqiBuffer))
                    IMCS = IMCSTable(4);
                else
                    IMCS = IMCSTable(udTemp.cqiBuffer+1);
                end
            
                [ITBS, modScheme] = lteMCS(IMCS);

                nRBs = length(udTemp.ulConfig.PUSCH.PRBSet);
                TBS  = lteTBS(nRBs, ITBS);
                udTemp.ulConfig.TBS = TBS;
                udTemp.ulConfig.Modulation = modScheme;
                udTemp.ulConfig.PUSCH.Modulation = modScheme;
                udTemp.ulConfig.PUSCH.TrBlkSizes = double(TBS) * ones(1, 10);
                
                udTemp.ulConfig.HARQ.harqProcesses(udTemp.ulConfig.HARQ.harqID).trBlkSize = double(TBS); 
            end
        end
        % *************************************************************** %
        
        %% Check SR from users and add/remove them to scheduler
        % Note that: Pools configuration is performed automatically every
        % PSCCH period (e.g., 40 msec).
        function checkSR_fn(obj)
            global gClk
            global simParameters
            
            for u_ = 1:length(obj.udsAttached)
                udTemp = obj.udsAttached(u_);
                
                switch(udTemp.nodeType)
                    % UE: check if 4 subframes passed and add to ul scheduler
                    case 'UE'
                        if(udTemp.ulConfig.SR == 1)
                            if(udTemp.nsrUE >= 3)
                                obj.enbScheduler.addUser_fn(udTemp);
                                
                                udTemp.nsrUE = 0;
                                udTemp.ulConfig.SR = 0;
                            else
                                udTemp.nsrUE = udTemp.nsrUE + 1;
                            end
                        else
                            obj.enbScheduler.removeUser_fn(udTemp);
                            udTemp.nsrUE = 0;
                            udTemp.ulConfig.SR = 0;
                        end
                end
            end
            
            % Pools configuration
            if(mod(gClk, 40) == 0)
               obj.ulPool = 1:(simParameters.nRBs/2);
               obj.slPool = 1:(simParameters.nRBs/2);
            end
        end
        % *************************************************************** %
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        %% eNB DL Tx to SBSs
        function eNB_DLTx_fn(obj)
            
            global MyGV
            
            obj.eNB_lteDLConfig.enb.NSubframe = mod(MyGV.Clock, 10);
            
            % Perform Scheduling to get the DL/UL resource allocation
            DLRG = zeros(obj.eNB_lteDLConfig.enb.NDLRB, 1);
            ULRG = zeros(obj.eNB_lteULConfig.ud.NULRB, 1);
            
            [obj.eNB_Pkt.DLRG, obj.eNB_Pkt.ULRG, UD_id, ~, ~, obj.Next_SBS] =...
                obj.enbScheduler.Schedule(DLRG, ULRG, 0, 0, obj.Next_SBS);
            
            % Compute the Imcs, MCS, and TBS, update each UD buffer and
            % construct the data for each UD.
            pdschSymbols        = [];
            pdschIndices        = [];
            dlschTransportBlk   = [];
            
            if(~isempty(UD_id))
                for u_ = 1 : length(UD_id)
                    
                    % Get the current user
                    UD_temp                       = obj.udsAttached(UD_id(u_));
                    UD_temp.UD_Pkt.DLRG           = find(obj.eNB_Pkt.DLRG == UD_temp.nodeID) - 1;
                    UD_temp.UD_Pkt.ULRG           = find(obj.eNB_Pkt.ULRG == UD_temp.nodeID) - 1;
                    UD_temp.UD_Pkt.NodeScheduled  = 1;
                    
                    obj.eNB_lteDLConfig.enb.PDSCH.PRBSet = UD_temp.UD_Pkt.DLRG;
                    
                    % Update the TBS based on the LTE standard tables
                    IMCSTable   = [-1 0 0 2 4 6 8 11 13 16 18 21 23 25 27 27];
                    if(~isempty(UD_temp.UD_Pkt.cqi))
                        CQI     = UD_temp.UD_Pkt.cqi(1);                % Picking the wideband cqi only
                    else
                        CQI     = [];
                    end
                    
                    if(isempty(CQI))
                        IMCS    = IMCSTable(4);
                    else
                        IMCS    = IMCSTable(CQI+1);
                    end
                    
                    [ITBS, modscheme]       = lteMCS(IMCS);
                    TBS                     = lteTBS(length(UD_temp.UD_Pkt.DLRG), ITBS);
                    
                    switch(modscheme)
                        case 'QPSK'
                            obj.eNB_lteDLConfig.M     = 2;
                        case '16QAM'
                            obj.eNB_lteDLConfig.M     = 4;
                        case '64QAM'
                            obj.eNB_lteDLConfig.M     = 6;
                        case '256QAM'
                            obj.eNB_lteDLConfig.M     = 8;
                    end
                    
                    obj.eNB_lteDLConfig.enb.Modulation          = modscheme;
                    obj.eNB_lteDLConfig.enb.PDSCH.Modulation    = modscheme;
                    
                    dlschTransportBlk_temp  = randi([0, 1], TBS, 1);
                    
                    [pdschIndices_temp, pdschInfo] = ...
                        ltePDSCHIndices(obj.eNB_lteDLConfig.enb,...
                        obj.eNB_lteDLConfig.enb.PDSCH, obj.eNB_lteDLConfig.enb.PDSCH.PRBSet, {'1based'});
                    
                    codedTrBlkSize  = pdschInfo.Gd * obj.eNB_lteDLConfig.M;
                    
                    codedTrBlock = lteDLSCH(obj.eNB_lteDLConfig.enb,...
                        obj.eNB_lteDLConfig.enb.PDSCH, codedTrBlkSize, dlschTransportBlk_temp);
                    
                    pdschSymbols_temp = ltePDSCH(obj.eNB_lteDLConfig.enb,...
                        obj.eNB_lteDLConfig.enb.PDSCH, codedTrBlock);
                    
                    pdschIndices        = [pdschIndices; pdschIndices_temp];
                    pdschSymbols        = [pdschSymbols; pdschSymbols_temp];
                    dlschTransportBlk   = [dlschTransportBlk; dlschTransportBlk_temp];
                    
                    UD_temp.UD_Pkt.DLIndices = pdschIndices_temp;
                    UD_temp.UD_Pkt.TxBits    = dlschTransportBlk_temp;
                    UD_temp.UD_Pkt.DL_TBS    = TBS;
                    UD_temp.UD_Pkt.UL_TBS    = TBS;
                    UD_temp.UD_Pkt.DL_QamMod = modscheme;
                    UD_temp.UD_Pkt.UL_QamMod = modscheme;
                end
                
                % Construct the DL Tx Waveform (Grid, and OFDM modulation)
                [DLTxWaveform, TxGrid, DLTxInfo]    = lteDLPhyTx(obj.eNB_lteDLConfig,...
                    pdschIndices, pdschSymbols);
                
                % Packet to send to the UDs, before applying the channel.
                obj.eNB_Pkt.TxGrid      = TxGrid;
                obj.eNB_Pkt.TxWaveform  = DLTxWaveform;
                obj.eNB_Pkt.TxBits      = dlschTransportBlk;
                obj.eNB_Pkt.TxInfo      = DLTxInfo;
                obj.eNB_Pkt.TxPktNumber = 0;
                
                % empty other attached users to this eNB
                DevicesIDs = [];
                
                for u_ = 1:length(obj.udsAttached)
                    if(~ismember(obj.udsAttached(u_).nodeID, UD_id))
                        DevicesIDs =  [DevicesIDs, obj.udsAttached(u_).nodeID];
                    end
                end
                
                if(~isempty(DevicesIDs))
                    for u_ = DevicesIDs
                        obj.udsAttached(u_).UD_Pkt.DLIndices   = [];
                        obj.udsAttached(u_).UD_Pkt.TxBits      = [];
                        obj.udsAttached(u_).UD_Pkt.DL_TBS      = [];
                        obj.udsAttached(u_).UD_Pkt.UL_TBS      = [];
                        obj.udsAttached(u_).UD_Pkt.DL_QamMod   = [];
                        obj.udsAttached(u_).UD_Pkt.UL_QamMod   = [];
                        obj.udsAttached(u_).UD_Pkt.ULRG        = [];
                        obj.udsAttached(u_).UD_Pkt.DLRG        = [];
                    end
                end
            end
        end
        % *************************************************************** %
        
        
        
        
        
        %% Distance between the object and a specified node
        function distance = disttoNode_fn(obj, Node)
            distance = sqrt(sum((obj.nodePosition - Node.nodePosition).^2));
        end
        % *************************************************************** %
        
        
        
        %% ApplyChannel: Apply pathloss and shadowing to both main link and interferers
        %         function [TotalWaveform, RxGrid, sinr_dB, cqi] = ApplyULChannel_fn(obj, SBS_temp)
        %
        %             global ulChannel
        %             global simParameters
        %
        %             % Recieved waveform
        %             TxWaveform  = SBS_temp.SBS_eNB_Pkt.TxWaveform;
        %             TxGrid      = SBS_temp.SBS_eNB_Pkt.TxGrid;
        %
        %             [RxWaveform, RxGrid]  = ulChannel.ApplyPathLoss_fn(...
        %                 TxWaveform, TxGrid, SBS_temp.disttoCell);
        %
        %             % Interferers
        %             Interfering_SBSs    = obj.InterferingSBSs;
        %
        %             IntRxWaveform   = zeros(size(RxWaveform));
        %             IntRxGrid       = zeros(size(TxGrid));
        %
        %             for i = 1 : length(Interfering_SBSs)
        %                 IntTxWaveform  = Interfering_SBSs(i).SBS_eNB_Pkt.TxWaveform;
        %                 IntTxGrid      = Interfering_SBSs(i).SBS_eNB_Pkt.TxGrid;
        %                 Int_Distance    = disttoNode_fn(obj, Interfering_UDs(i));
        %
        %                 [IntRxWaveform_temp, IntRxGrid_temp]  = ulChannel.ApplyPathLoss_fn(...
        %                     IntTxWaveform, IntTxGrid, Int_Distance);
        %
        %                 if(~isempty(IntRxWaveform_temp))
        %                     IntRxWaveform   = IntRxWaveform + IntRxWaveform_temp;
        %                     IntRxGrid       = IntRxGrid + IntRxGrid_temp;
        %                 end
        %             end
        %
        %             % White noise
        %             No_dB   = -174 + 10 * log10(simParameters.bandwidth * 1e6) - 30;
        %             No      = 10^(No_dB/10);
        %             noise   = No * complex(randn(size(RxWaveform)), randn(size(RxWaveform)));
        %
        %             TotalWaveform = RxWaveform + IntRxWaveform + noise;
        %
        %             % SNR estimation
        %             [sinr_Wb, sinr_Sb] = ComputeSINR(RxWaveform, IntRxWaveform, 5, ...
        %                 simParameters.bandwidth*1e6, SBS_temp.SBS_eNB_Pkt.TxInfo);    % sinr wideband and subband
        %
        %             sinr_dB = [sinr_Wb; sinr_Sb];
        %
        %             % CQI
        %             cqi         = snr_cqi(sinr_Wb);
        %         end
        %         % *************************************************************** %
        
        
        
        %% Interfering UDs
        function Init_InterferingSBSs_fn(obj)
            global SBSs
            
            obj.InterferingSBSs     = [];
            SBS_attched             = obj.SBSs_Attached;
            
            for sbs_ = 1:length(SBSs)
                SBS_temp    = SBSs(sbs_);
                Distance    = sqrt(sum((obj.nodePosition - SBS_temp.nodePosition).^2));
                
                if((~ismember(SBS_temp, SBS_attched)) && (Distance < (obj.cellRadius)))
                    obj.InterferingSBSs = [obj.InterferingSBSs, SBS_temp];
                end
            end
        end
        % *************************************************************** %
        
        %% BER calculation
        function [BitErrors, BER] = BER_fn(obj)
            BitErrors           = sum(obj.eNB_Pkt.RxBits ~= obj.eNB_Pkt.TxBits);
            BER                 = BitErrors / length(obj.eNB_Pkt.TxBits);
            obj.eNB_Pkt.BER     = BER;
        end
        % *************************************************************** %
        
        %% Distance between two nodes
        function Distance = Euclide_Distance_fn(obj, Node1, Node2)
            Node1Pos = [Node1.nodePosition]';
            Node2Pos = reshape([Node2.nodePosition], 2, []);
            Distance = sqrt(sum((Node1Pos - Node2Pos).^2, 1));
            Distance(Distance == 0) = [];
        end
        % *************************************************************** %
    end
end











