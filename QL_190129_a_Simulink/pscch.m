

classdef pscch < matlab.mixin.Copyable
    
    properties (SetAccess = private)
              
        NSubframeBegin;         % First subframe of PSCCH period ('jbegin'). Read-only
        PeriodLength;           % Number of subframes in pool period. Read-only
        TxMode;                 % PSSCH transmission mode ('Mode1','Mode2'). Read-only
        
        PSCCHSubframePool;      % PSCCH subframe pool indices. Read-only
        PSCCHResourceBlockPool; % PSCCH resource block pool indices. Read-only
        NumPSCCHResource;       % Number of nPSCCH values (0... NumPSCCHResource). Read-only
               
        PSSCHSubframePool;      % PSSCH subframe pool indices. Read-only
        PSSCHResourceBlockPool; % PSSCH resource block pool indices. Read-only
        AllowedITRP;            % Allowed ITRP values. Read-only
        
        SyncSubframes;          % PSSS/SSSS/PSBCH synchronization subframe indices. Read-only
        
    end
    
    properties (Access = public, SetObservable = true)
        Config;
    end
    
    properties (Access=private)
        
        % PSCCH resource pool related intermediate variables
        L_PSCCH;            % Number of subframes in PSCCH subframe pool
        M_RB_PSCCH_RP;      % Number of RB in PSCCH resource block pool

    end
    
    methods (Access = public, Static=true)
        
        %% Default configuration
        function p = defaultConfig(cid,bw)
            
            % Default NSLRB value (number of sidelink resource blocks)
            nslrb = 25;
            
            if nargin > 0
                % Validate the configuration ID input (if present)
                cidset = [1 2 5];
                if ~any(cid == cidset)
                    error('The default configuration identifier is not one of the set [%s\b\b]',sprintf('%d, ',cidset));
                end 
                if nargin>1
                    % Validate and convert the bandwidth input parameter (if present) into the associated NSLRB
                    bwset = {'1.4MHz','3MHz','5MHz','10MHz','15MHz','20MHz'};
                    bwrb = [6 15 25 50 75 100];
                    nslrb = bwrb(strcmpi(bw,bwset));
                    if isempty(nslrb)
                        error('The default configuration bandwidth is not one of the set [%s]',strjoin(bwset,', '));
                    end
                end
            else
                cid = 5;    % Default configuration ID
            end
         
            % TS 36.101 Section A.7.2.1 FDD,
            % Table A.7.2.1-1: ProSe Direct Communication configuration for E-UTRA FDD (Configuration #1-FDD)
            p = struct();
            
            % General settings (following LTE System Toolbox parameter naming where applicable)
            p.NSLRB = nslrb;
            p.DuplexMode = 'FDD';
            p.TDDConfig = 0;
            p.UESelected = 'On';
            p.SyncEnable = 'On';
            p.NPSCCHPeriod = 0;
            
            % SL-CommResourcePool part            
            p.sc_CP_Len_r12 = 'Normal';        % 'Normal', 'Extended'
            p.sc_Period_r12 = 40;              % 40,60,70,80,120,140,160,240,280,320 subframes

            % Control transport format 
            tf = struct();
            tf.prb_Num_r12 = ceil(p.NSLRB/2);
            tf.prb_Start_r12 = 0;
            tf.prb_End_r12 = p.NSLRB-1;
            tf.offsetIndicator_r12 = 0;
            tf.subframeBitmap_r12 = '0001100000000000000000000000000000000000';      % Bitmap of length 4,8,12,16,30,40,42 bits
            p.sc_TF_ResourceConfig_r12 = tf;
                        
            p.data_CP_Len_r12 = 'Normal';

            dh = struct(); 
            dh.hoppingParameter_r12 = 504;
            dh.numSubbands_r12 = 2;             % 1,2,4 subbands
            dh.rb_Offset_r12 = 0;          
            p.dataHoppingConfig_r12 = dh;

            % Data transport format (mode 2) - UE selected part, or PSSCH resource configuration (36.213)
            tf = struct();
            tf.prb_Num_r12 = ceil(p.NSLRB/2);
            tf.prb_Start_r12 = 0;
            tf.prb_End_r12 = p.NSLRB-1;
            tf.offsetIndicator_r12 = 0;
            tf.subframeBitmap_r12 = '0000000011111111111111110000000000000000'; % Bitmap of length 4,8,12,16,30,40,42 bits
            p.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12 = tf;
            p.ue_SelectedResourceConfig_r12.trpt_Subset_r12 = '010';
       
            sc = struct();                          % Instance of SL-SyncConfig... would be referenced by index
            sc.syncCP_Len_r12 = 'Normal';           % 'Normal','Extended'
            sc.syncOffsetIndicator_r12 = 1;         % 0...39
            sc.slssid_r12 = 0;
            p.syncConfig = sc;
 
            % TS 36.101 Section A.7.2.1 FDD,            
            % Table A.7.2.1-2: ProSe Direct Communication configuration for E-UTRA FDD (Configuration #2-FDD)           
            if cid==2
                p.SyncEnable = 'Off';
                p.sc_TF_ResourceConfig_r12.subframeBitmap_r12 = '0011110000000000000000000000000000000000';
                p.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.subframeBitmap_r12 = '0000000011111111111111110000000000000000';
            % 
            % Table A.7.2.1-5: ProSe Direct Communication configuration for E-UTRA FDD (Configuration #5-FDD)           
            elseif cid==5
                p.SyncEnable = 'Off';
                p.sc_TF_ResourceConfig_r12.subframeBitmap_r12 = '0001100000000000000000000000000000000000';
                p.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.subframeBitmap_r12 = '0000000011111111111111111111111111111111'; 
                p.ue_SelectedResourceConfig_r12.trpt_Subset_r12 = '001';
            end
                   
        end
              
    end
                       
    methods
        
        % Class constructor
        function obj = pscch(params)
            
            % Default parameter set, if no input provided
            if nargin == 0
                params = PSCCHPeriod.defaultConfig();
%                 params = PSCCHPeriod;
            end

            % Set the Config property to the input parameter set
            obj.Config = params;
            
        end    
        
        %% Setter associated with the Config parameter
        function set.Config(obj,params)
            
            % Cache the current Config state in the event that we have 
            % to roll back due to any errors in the new parameter set
            current = obj.Config;
            obj.Config = params;
            try
                obj.initialize();
            catch e
                obj.Config = current;     % Roll back to last known good
                obj.initialize();         % And re-initialize
                rethrow(e);
            end
          
        end
        
        %% generatewaveform
        function [waveform,grid] = genWaveform(obj, nodeID, dci, TBs, sfNumber)
           
            % update subframe number
            sfNumber = mod(sfNumber, obj.Config.sc_Period_r12)+1;
            
            % Preallocate output for PSCCH period resource grid and baseband waveform
            ue.NSLRB = obj.Config.NSLRB;
            ue.CyclicPrefixSL = obj.Config.syncConfig.syncCP_Len_r12;
            sfdims = lteSLResourceGridSize(ue);
            grid = zeros(sfdims(1), sfdims(2));
            waveinfo = lteSLSCFDMAInfo(ue);
            nsamplespersf = waveinfo.SamplingRate/1000;
            waveform = zeros(nsamplespersf, 1);
            
            % The SC-FDMA modulation will not use windowing due to the variable CP
            % and the piece-wise nature of the waveform construction
            windowing = 0;
                                                          
            % Create basic set of common parameters used by the low-level PHY functions
            ue = struct('NSLRB', obj.Config.NSLRB, 'CyclicPrefixSL', obj.Config.sc_CP_Len_r12,...
                'DuplexMode', obj.Config.DuplexMode, 'TDDConfig', obj.Config.TDDConfig);
            
            ue.SidelinkMode = 'D2D';
            
            % Get the resources (PRB/subframes) that will be used to carry
            % the pair of PSCCH which in turn carry the coded SCI
            [pscchsubframes, pscchprb] = obj.getPSCCHResources(dci);
            [pscchsubframes, Idx] = sort(pscchsubframes);
            pscchprb = pscchprb(Idx);
            pscchsubframes = pscchsubframes + 1;  % Make subframe indices 1-based for container indexing
            
            % Calculate the transport/physical channel content that is sent in both PSCCH subframes
            % Since the content doesn't change (and the DRS are independent of the subframe)
            % the only difference between the PSCCH is the PRB used for transmission 
            dci.SCIFormat = 'Format0';
            [sci, scibits] = lteSCI(ue, dci);
            cw = lteSCIEncode(ue, scibits);
            symbols = ltePSCCH(cw);
            drssymbols = ltePSCCHDRS();
            
            % Create an empty grid for PSCCH and PSSCH symbols
            lgrid = lteSLResourceGrid(ue); 
                    
            % PSCCH Control instances
            
            % If subframe carries sync then skip further creation
            if sfNumber ~= (obj.SyncSubframes+1)
                nsfpscch = find(sfNumber == pscchsubframes);
                
                % create the pscch control
                if(~isempty(nsfpscch))
                    fprintf("User: %d \t Txed PSCCH subframe = %d\n", nodeID, sfNumber);
                    
                    % Update the PRB allocation parameter for this PSCCH
                    ue.PRBSet = pscchprb(:, nsfpscch);      

                    % Write PSCCH
                    pscchindices = ltePSCCHIndices(ue);
                    lgrid(pscchindices) = symbols;
                    
                    % Write PSCCH DRS 
                    pscchdrsindices = ltePSCCHDRSIndices(ue);
                    lgrid(pscchdrsindices) = drssymbols;
                    
                    % Write subframe into overall period grid
                    % Modulate subframe and write the waveform section into the overall waveform vector  
                    waveform = lteSLSCFDMAModulate(ue, lgrid, windowing);   
                end
            end
            
            % Create PSSCH subframe sequence - sent on antenna port 1000      
                        
            % Get the resources (PRB/subframes) that will be used to carry
            % the sequence of PSSCH which in turn carry a sequence of encoded transport blocks
            [psschsubframes, psschprb, poolindices] = obj.getPSSCHResources(sci);
            psschsubframes = psschsubframes + 1;  % Make subframe indices 1-based for container indexing purposes                                        
                
            % Set up processing parameters that are specific to this PSSCH sequence      
            ue.NSAID = sci.NSAID;
            ue.CyclicPrefixSL = obj.Config.data_CP_Len_r12;     
            % Configure PSSCH modulation parameter from DCI/SCI MCS index
            [itbs, ue.Modulation] = lteMCS(sci.ModCoding, 'PUSCH');
            if (strcmpi(ue.Modulation, '64QAM'))
                ue.Modulation = '16QAM';
            end
            
            % And calculate transport block size from MCS too
            nprb = max(size(psschprb, 1), 1);    % It's possible that no PSSCH resources can be allocated (for example, to satisfy the multiple of 4 rule)
            tbs = lteTBS(nprb,itbs);
        
            % Each transport block created and sent with the period is 
            % transmitted four times according to a fixed RV sequence 
            rvsequence = [0 2 3 1];
            
            % The processing for each PSSCH instance is dependent of the RV value
            % and NSubframePSSCH number (subframe index relative to the PSSCH subframe pool)
            % therefore there are no intermediate vectors to pre-calculate and cache
            if sfNumber ~= (obj.SyncSubframes+1)
                nsfpssch = find(sfNumber == psschsubframes);
                
                if(~isempty(nsfpssch)) 
                    fprintf("User: %d \t Txed PSSCH subframe = %d\n", nodeID, sfNumber);
                    % Calculate the RV value associated with this PSSCH and 
                    % whether a new transport block is to be created
                    % Always perform this step even if the subframe is configured for sync 
                    ue.RV = rvsequence(mod(nsfpssch-1, 4)+1); 

                    nTBs = length(TBs);
                    if nTBs == 1
                        trblkinfobits = cell2mat(TBs);
                    else
                        trblkinfobits = TBs{mod(nsfpssch-1, nTBs)+1}; 
                    end

                    % Update the PRB allocation for this PSSCH
                    ue.PRBSet = psschprb(:, nsfpssch);
                    
                    % Update the PSSCH subframe number (wrt the PSSCH subframe pool)
                    ue.NSubframePSSCH = nsfpssch-1;      
                    
                    % Write PSSCH (encode transport block with new RV)
                    [psschindices, info] = ltePSSCHIndices(ue); 
                    cw = lteSLSCH(ue, info.G, trblkinfobits);
                    symbols = ltePSSCH(ue,cw);
                    lgrid(psschindices) = symbols;
                    
                    % Write PSSCH DRS
                    drssymbols = ltePSSCHDRS(ue);
                    psschdrsindices = ltePSSCHDRSIndices(ue);
                    lgrid(psschdrsindices) = drssymbols;

                    % Write subframe into overall period grid
                    % Modulate subframe and write the waveform section into the overall waveform vector  
                    waveform = lteSLSCFDMAModulate(ue, lgrid, windowing);  
                end 
            end

            % Synchronization subframes
            if sfNumber == (obj.SyncSubframes+1)
                nsfsync = find(sfNumber == obj.SyncSubframes);
                
                if ~isempty(nsfsync)

                    % Set up processing parameters that are specific to the sync sequence         
                    ue.CyclicPrefixSL = obj.Config.syncConfig.syncCP_Len_r12;
                    ue.NSLID = obj.Config.syncConfig.slssid_r12;   

                    syncsubframes = obj.SyncSubframes + 1;  % Make subframe indices 1-based for container indexing    

                    % Calculate the transport/physical channel content that is sent in all sync subframes
                    % The content doesn't change across these subframes except for the PSBCH since 
                    % the NSubframe/NFrame fields of the SL-MIB was change according to the subframe use
                    % Therefore we can pre-create a complete sync subframe and write the PSBCH as we go  

                    % PSSS symbols and indices
                    pssssymbols = sqrt(72/62)*ltePSSS(ue);
                    psssindices = ltePSSSIndices(ue);

                    % SSSS symbols and indices
                    sssssymbols = sqrt(72/62)*lteSSSS(ue);
                    ssssindices = lteSSSSIndices(ue);

                    % PSBCH indices and PSBCH DRS symbols and indices
                    psbchindices = ltePSBCHIndices(ue); 
                    psbchdrssymbols = ltePSBCHDRS(ue);
                    psbchdrsindices = ltePSBCHDRSIndices(ue);

                    % Write all fixed signals into an empty subframe grid
                    lgrid = lteSLResourceGrid(ue);
                    lgrid(psssindices) = pssssymbols;
                    lgrid(ssssindices) = sssssymbols;
                    lgrid(psbchdrsindices) = psbchdrssymbols;

                    % Processing each sync subframe instance in the period by updating the subframe/frame number
                    % in the MIB-SL message and recoding through the SL-BCH and onto the PSBCH
                    for i=1:length(syncsubframes)

                        % Create MIB-SL (letting unset fields default to 0)
                        nsf = obj.NSubframeBegin + syncsubframes(i)-1; % Current absolute subframe number
                        ue.NSubframe = mod(nsf,10);         % 0...9
                        ue.NFrame = floor(nsf/10);          % Will be 0...1024-1 in the message since it's 10 bits in the SL-MIB
                        ue.InCoverage = 1;                  % Indicates UE is in E-UTRAN coverage
                        mibsl = lteSLMIB(ue);               % Turn parameters into the 40 bit MIB-SL

                        % Encode MIB-SL on the SL-BCH and write PSBCH into subframe grid
                        codeword = lteSLBCH(ue,mibsl);
                        psbchsymbols = ltePSBCH(ue,codeword);
                        lgrid(psbchindices) = psbchsymbols;

                        % Write subframe into overall period grid
                        % Modulate subframe and write the waveform section into the overall waveform vector  
                        waveform = lteSLSCFDMAModulate(ue, lgrid, windowing);  
                    end
                end
            end
        end
        % *************************************************************** %

        %% generateWaveform
        function [waveform,grid] = generateWaveform(obj, dci, TBs, sfNumber)
            % Input argument check 
            narginchk(4,4);
        
            % The waveform generation is performed by creating the sequences of PSCCH, 
            % PSSCH and synchronization subframes separately and writing the individual
            % subframes into the output grid and waveform vector
                      
            % The cyclic prefix of the subframes in the PSCCH period can vary
            % depending on the channel type, so to help preallocate variables 
            % create a list of the number of SC-FDMA symbols in each of the 
            % subframes in the PSCCH period 
            nsubframes = obj.Config.sc_Period_r12;
            cplengths = 14*ones(1,nsubframes);

            % update subframe number
            sfNumber = mod(sfNumber, obj.Config.sc_Period_r12);
            
            % Mark the PSCCH 
            ue.NSLRB = obj.Config.NSLRB;
            ue.CyclicPrefixSL = obj.Config.sc_CP_Len_r12;
            sfdims = lteSLResourceGridSize(ue);
            cplengths(obj.PSCCHSubframePool+1) = sfdims(2);

            % Mark the PSSCH 
            ue.CyclicPrefixSL = obj.Config.data_CP_Len_r12;
            sfdims = lteSLResourceGridSize(ue);
            cplengths(obj.PSSCHSubframePool+1) = sfdims(2);

            % Overwrite the synchonization subframes
            ue.CyclicPrefixSL = obj.Config.syncConfig.syncCP_Len_r12;
            sfdims = lteSLResourceGridSize(ue);
            cplengths(obj.SyncSubframes+1) = sfdims(2);

            % Calculate the total number of SC-FDMA symbols across all 
            % the subframes in the PSCCH period
            sfcppos = [0 cumsum(cplengths)];
            nsymbols = sfcppos(end);
            sfcppos = sfcppos + 1;
            
            % Preallocate output for PSCCH period resource grid and baseband waveform
            grid = zeros(sfdims(1),nsymbols);
            waveinfo = lteSLSCFDMAInfo(ue);
            nsamplespersf = waveinfo.SamplingRate/1000;
            waveform = zeros(nsubframes*nsamplespersf,1);
            
            % Mark all the sync frames (if any) in a vector
            % This will be used to skip these subframes during PSCCH/PSSCH processing 
            sfmap = zeros(1,nsymbols);
            sfmap(obj.SyncSubframes+1) = 1;
            
            % The SC-FDMA modulation will not use windowing due to the variable CP
            % and the piece-wise nature of the waveform construction
            windowing = 0;
                                              
            % Create PSCCH subframe sequence - sent on antenna port 1000                          
            
            % Create basic set of common parameters used by the low-level PHY functions
            ue = struct('NSLRB',obj.Config.NSLRB,'CyclicPrefixSL',obj.Config.sc_CP_Len_r12,'DuplexMode',obj.Config.DuplexMode,'TDDConfig',obj.Config.TDDConfig);
            ue.SidelinkMode = 'D2D';
            
            % Get the resources (PRB/subframes) that will be used to carry
            % the pair of PSCCH which in turn carry the coded SCI
            [pscchsubframes,pscchprb] = obj.getPSCCHResources(dci);
            pscchsubframes = pscchsubframes + 1;  % Make subframe indices 1-based for container indexing
            
            % Calculate the transport/physical channel content that is sent in both PSCCH subframes
            % Since the content doesn't change (and the DRS are independent of the subframe)
            % the only difference between the PSCCH is the PRB used for transmission 
            dci.SCIFormat = 'Format0';
            [sci,scibits] = lteSCI(ue,dci);
            cw = lteSCIEncode(ue,scibits);
            symbols = ltePSCCH(cw);
            drssymbols = ltePSCCHDRS();
   
            % PSCCH Control instances...
            for i=1:length(pscchsubframes)               
                
                % If subframe carries sync then skip further creation
                if sfmap(pscchsubframes(i))
                    continue;
                end
               
                % Update the PRB allocation parameter for this PSCCH
                ue.PRBSet = pscchprb(:,i);
                
                % Write PSCCH and PSCCHDRS symbols into an empty subframe grid
                lgrid = lteSLResourceGrid(ue); 
                % Write PSCCH
                pscchindices = ltePSCCHIndices(ue);
                lgrid(pscchindices) = symbols;
                % Write PSCCH DRS 
                pscchdrsindices = ltePSCCHDRSIndices(ue);
                lgrid(pscchdrsindices) = drssymbols;

                % Write subframe into overall period grid
                % Modulate subframe and write the waveform section into the overall waveform vector  
                grid(:,sfcppos(pscchsubframes(i)):sfcppos(pscchsubframes(i)+1)-1) = lgrid;
                waveform((pscchsubframes(i)-1)*nsamplespersf+1:pscchsubframes(i)*nsamplespersf) = lteSLSCFDMAModulate(ue,lgrid,windowing);          
            end
           
            % Create PSSCH subframe sequence - sent on antenna port 1000      
                        
            % Get the resources (PRB/subframes) that will be used to carry
            % the sequence of PSSCH which in turn carry a sequence of encoded transport blocks
            [psschsubframes,psschprb,poolindices] = obj.getPSSCHResources(sci);
            psschsubframes = psschsubframes + 1;  % Make subframe indices 1-based for container indexing purposes                                        
                
            % Set up processing parameters that are specific to this PSSCH sequence      
            ue.NSAID = sci.NSAID;
            ue.CyclicPrefixSL = obj.Config.data_CP_Len_r12;     
            % Configure PSSCH modulation parameter from DCI/SCI MCS index
            [itbs,ue.Modulation] = lteMCS(sci.ModCoding,'PUSCH');
            if (strcmpi(ue.Modulation,'64QAM'))
                ue.Modulation = '16QAM';
            end
            % And calculate transport block size from MCS too
            nprb = max(size(psschprb,1),1);    % It's possible that no PSSCH resources can be allocated (for example, to satisfy the multiple of 4 rule)
            tbs = lteTBS(nprb,itbs);
        
            % Each transport block created and sent with the period is 
            % transmitted four times according to a fixed RV sequence 
            rvsequence = [0 2 3 1];
            
            % The processing for each PSSCH instance is dependent of the RV value
            % and NSubframePSSCH number (subframe index relative to the PSSCH subframe pool)
            % therefore there are no intermediate vectors to pre-calculate and cache
            k = 1;
            for i=1:length(psschsubframes)
                
                % Calculate the RV value associated with this PSSCH and 
                % whether a new transport block is to be created
                % Always perform this step even if the subframe is configured for sync 
                ue.RV = rvsequence(mod(i-1,4)+1);
                if ue.RV==0
                    if k > length(TBs)
                        break;
                    else
                        trblkinfobits = TBs{k};%randi([0 1],1,tbs);
                        k = k+1;
                    end
                end

                % If subframe carries sync then skip further creation
                if sfmap(psschsubframes(i))
                    continue;
                end
          
                % Update the PRB allocation for this PSSCH
                ue.PRBSet = psschprb(:,i);
                % Update the PSSCH subframe number (wrt the PSSCH subframe pool)
                ue.NSubframePSSCH = poolindices(i);     
                               
                % Write PSSCH and PSSCHDRS symbols into an empty subframe grid
                lgrid = lteSLResourceGrid(ue);                   
                % Write PSSCH (encode transport block with new RV)
                [psschindices,info] = ltePSSCHIndices(ue); 
                cw = lteSLSCH(ue,info.G,trblkinfobits);
                symbols = ltePSSCH(ue,cw);
                lgrid(psschindices) = symbols;
                % Write PSSCH DRS
                drssymbols = ltePSSCHDRS(ue);
                psschdrsindices = ltePSSCHDRSIndices(ue);
                lgrid(psschdrsindices) = drssymbols;
                
                % Write subframe into overall period grid
                % Modulate subframe and write the waveform section into the overall waveform vector  
                grid(:,sfcppos(psschsubframes(i)):sfcppos(psschsubframes(i)+1)-1) = lgrid;
                waveform((psschsubframes(i)-1)*nsamplespersf+1:psschsubframes(i)*nsamplespersf) = lteSLSCFDMAModulate(ue,lgrid,windowing);  
            end               
                
            % Create sync subframe sequence (merged into same grid as PSCCH/PSSCH)
            % PSBCH sent on antenna port 1010
            % PSSS/SSSS sent on antenna port 1020
            
            if ~isempty(obj.SyncSubframes)
                
                % Set up processing parameters that are specific to the sync sequence         
                ue.CyclicPrefixSL = obj.Config.syncConfig.syncCP_Len_r12;
                ue.NSLID = obj.Config.syncConfig.slssid_r12;   
                       
                syncsubframes = obj.SyncSubframes + 1;  % Make subframe indices 1-based for container indexing    
                              
                % Calculate the transport/physical channel content that is sent in all sync subframes
                % The content doesn't change across these subframes except for the PSBCH since 
                % the NSubframe/NFrame fields of the SL-MIB was change according to the subframe use
                % Therefore we can pre-create a complete sync subframe and write the PSBCH as we go  
                
                % PSSS symbols and indices
                pssssymbols = sqrt(72/62)*ltePSSS(ue);
                psssindices = ltePSSSIndices(ue);

                % SSSS symbols and indices
                sssssymbols = sqrt(72/62)*lteSSSS(ue);
                ssssindices = lteSSSSIndices(ue);

                % PSBCH indices and PSBCH DRS symbols and indices
                psbchindices = ltePSBCHIndices(ue); 
                psbchdrssymbols = ltePSBCHDRS(ue);
                psbchdrsindices = ltePSBCHDRSIndices(ue);
           
                % Write all fixed signals into an empty subframe grid
                lgrid = lteSLResourceGrid(ue);
                lgrid(psssindices) = pssssymbols;
                lgrid(ssssindices) = sssssymbols;
                lgrid(psbchdrsindices) = psbchdrssymbols;
                
                % Processing each sync subframe instance in the period by updating the subframe/frame number
                % in the MIB-SL message and recoding through the SL-BCH and onto the PSBCH
                for i=1:length(syncsubframes)

                    % Create MIB-SL (letting unset fields default to 0)
                    nsf = obj.NSubframeBegin + syncsubframes(i)-1; % Current absolute subframe number
                    ue.NSubframe = mod(nsf,10);         % 0...9
                    ue.NFrame = floor(nsf/10);          % Will be 0...1024-1 in the message since it's 10 bits in the SL-MIB
                    ue.InCoverage = 1;                  % Indicates UE is in E-UTRAN coverage
                    mibsl = lteSLMIB(ue);               % Turn parameters into the 40 bit MIB-SL

                    % Encode MIB-SL on the SL-BCH and write PSBCH into subframe grid
                    codeword = lteSLBCH(ue,mibsl);
                    psbchsymbols = ltePSBCH(ue,codeword);
                    lgrid(psbchindices) = psbchsymbols;
 
                    % Write subframe into overall period grid
                    % Modulate subframe and write the waveform section into the overall waveform vector  
                    grid(:,sfcppos(syncsubframes(i)):sfcppos(syncsubframes(i)+1)-1) = lgrid;
                    waveform((syncsubframes(i)-1)*nsamplespersf+1:syncsubframes(i)*nsamplespersf) = lteSLSCFDMAModulate(ue,lgrid,windowing);  
                end
            
            end
            
        end
        
        %% displayPeriod
        function displayPeriod(obj,dci)
        %displayPeriod Display image representing the PSCCH period resource pools 
        %   displayPeriod(obj) displays an image representing the locations of 
        %   the subframe pools and physical resource block pools within the PSCCH
        %   period grid. The PSCCH resource pool is displayed in lighter blue and
        %   the PSSCH resource pool is given in yellow. Any active synchronization 
        %   resource blocks/subframes are overlaid on top of the pools.
        %
        %   displayPeriod(obj,dci) is the same as above except it also displays
        %   the specific PSCCH and PSSCH resources that are selected from the 
        %   respective pools for a scheduled data transmission. The used resources 
        %   are marked in green. The input structure, dci, contains the scheduling
        %   parameters that are found in a DCI format 5/SCI format 0 message. 
        %   This should include the following fields: 
        %   PSCCHResource       - PSCCH resource index ("Resource for PSCCH"). Optional.
        %   TimeResourcePattern - Time resource pattern index (I_TRP). Optional.
        %   FreqHopping         - Frequency hopping flag
        %                         (0 for non-hopping, 1 for hopping)
        %   Allocation          - Resource allocation parameter substructure
        %   The Allocation substructure should contain the following fields:
        %       Allocation.HoppingBits  - Hopping bits (0...2 bits). Only required for frequency hopping.
        %       Allocation.RIV          - Resource indication value (5...13 bits)
        %
        %   These parameters are describe further in <a href="matlab: doc('lteDCI')">lteDCI</a> and <a href="matlab: doc('lteSCI')">lteSCI</a>.
        %
        %   Examples: 
        %   % Display an image representing the PSCCH and PSSCH resource pools 
        %   % within a PSCCH period.
        %
        %   period = PSCCHPeriod;
        %   displayPeriod(period);
        % 
        %   % Display an image representing the resource blocks used by a set of
        %   % PSCCH and PSSCH transmissions selected from the resource pools within
        %   % a PSCCH period.
        % 
        %   dci.PSCCHResource = 0;
        %   dci.TimeResourcePattern = 106;
        %   dci.FreqHopping = 1;
        %   dci.Allocation.HoppingBits = 3;
        %   dci.Allocation.RIV = 0;
        %   displayPeriod(period,dci);
        %
        %   See also generateWaveform, lteSCIResourceAllocation, 
        %   lteDCIResourceAllocation.

            % Input argument check 
            narginchk(1,2);
  
            % Clear current figure for use by function
            clf;
 
            % Create 2-D arrays representing all the RB per subframe grid across PSCCH period
            % Create separate grids for the PSCCH/PSSCH and sync subframes, then these grids
            % will be combined later to display a final image 
            commgrid = zeros(obj.Config.NSLRB,obj.Config.sc_Period_r12);
            syncgrid = commgrid;
                       
            % Mark the sync RB/subframes in the PSCCH period grid
            mp = fix(obj.Config.NSLRB/2);
            srb = mp-2:mp+3+mod(obj.Config.NSLRB,2);        % The range of RB affected by sync channels
            syncgrid(srb,obj.SyncSubframes+1) = 20;             % Mark with lighter blue - Sync                 
                             
            % Mark PSCCH/PSSCH RB/subframe pool in the PSCCH period grid 
            commgrid(obj.PSCCHResourceBlockPool+1,obj.PSCCHSubframePool+1) = 10;   % Mark with blue - Control pool
            commgrid(obj.PSSCHResourceBlockPool+1,obj.PSSCHSubframePool+1) = 60;   % Mark with orange/yellow - Shared Pool

            % Depending on whether a DCI grant was input, mark the grid with the resource allocation 
            % associated with the grant and create allocation specific annotation text, 
            % otherwise create annotation text describing the pools only
            if nargin == 1
                            
                % Create cell array of text (each cell is a new line) describing the resource pools
                % First, include some text about the transmission mode and create the PSCCH resource
                % pool description
                str1 = {...
                       sprintf('PSSCH %s',obj.TxMode),...
                       'PSCCH Pool Config:',...
                       sprintf('  Offset: %d',obj.Config.sc_TF_ResourceConfig_r12.offsetIndicator_r12),...
                       sprintf('  Bitmap: %s',num2str(logicalBitmap(obj.Config.sc_TF_ResourceConfig_r12.subframeBitmap_r12),'%d')),...      
                       sprintf('  PRB Start,End,Num: [%d,%d,%d]',obj.Config.sc_TF_ResourceConfig_r12.prb_Start_r12,...
                                                                 obj.Config.sc_TF_ResourceConfig_r12.prb_End_r12,...
                                                                 obj.Config.sc_TF_ResourceConfig_r12.prb_Num_r12)...
                };
                % If transmission mode 2 then add additional information about the PSSCH resource pool 
                if strcmpi(obj.TxMode,'Mode1')
                    str2 = '  Mode1';
                else
                  str2 = {...              
                       'PSSCH Pool Config:',...
                       sprintf('  Offset: %d',obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.offsetIndicator_r12),...
                       sprintf('  Bitmap: %s',num2str(logicalBitmap(obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.subframeBitmap_r12),'%d')),... 
                       sprintf('  PRB Start,End,Num: [%d,%d,%d]',obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_Start_r12,...
                                                                 obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_End_r12,...
                                                                 obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_Num_r12)... 
                       };                
                end

                % Combine both arrays of text (PSCCH and PSSCH) to give final annotation text 
                annostr = [str1 str2];
                fontsize = 7;
                tboxdim = [0.70 0.6 0.3 0.3];
                
            else
                % If a DCI grant was input then derive the associated RB/subframes
                % and mark them on the PSCCH/PSSCH grid
                % In this case the annotation text will only describe the grant and not the pools
                
                % Get the PSCCH resources (RB used per subframe) and mark them on the grid
                [subframesC,rbsC] = obj.getPSCCHResources(dci);          
                commgrid(sub2ind(size(commgrid),rbsC+1,subframesC+1)) = 40; % Mark with green

                % Get the PSSCH resources and mark them on the grid
                [subframesS,rbsS] = obj.getPSSCHResources(dci);
                commgrid(sub2ind(size(commgrid),rbsS+1,kron(subframesS+1,ones(size(rbsS,1),1)))) = 40;  % Mark with green

                % Create cell array of text describing the DCI grant 
                if isfield(dci.Allocation,'HoppingBits')
                    hbstr =  num2str(dci.Allocation.HoppingBits);
                else
                    hbstr = 'N/A';
                end       
                annostr = {...
                       sprintf('nPSCCH: %d',dci.PSCCHResource),...
                       sprintf('ITRP: %d',dci.TimeResourcePattern),...
                       sprintf('Freq Hopping: %d',dci.FreqHopping),...
                       sprintf('Hopping Bits Value: %s',hbstr),...
                       sprintf('RIV: %d',dci.Allocation.RIV)...
                       };
                fontsize = 8; 
                tboxdim = [0.75 0.6 0.3 0.3];
                
            end
                 
            % Create image of the additively combined the two grids 
            % and giving the figure a title and labels
            image(commgrid+syncgrid);axis xy;
            tdstr = [];
            if strcmpi(obj.Config.DuplexMode,'TDD')
                tdstr = sprintf(' Config %d',obj.Config.TDDConfig);
            end
            ttext = sprintf('PSCCH period #%d (%s%s)',obj.Config.NPSCCHPeriod,obj.Config.DuplexMode,tdstr);
            title(ttext); xlabel('Subframes'); ylabel('Resource blocks');
            
            % Add annotation text box to image figure
            % For reference, default axes position is [0.1300 0.1100 0.7750 0.8150], [lowerleftx lowerlefty width height] 
            t = annotation('textbox',tboxdim,'String',annostr,'FitBoxToText','on','BackgroundColor','white');
            t.FontSize = fontsize;
            t.FaceAlpha = 0.85;
            
            % Make figure visible
            figure(gcf);       
            
        end
             
        function [riv,range] = getAllowedRIV(obj,hoppingbits)
        % getAllowedRIV Get allowed RIV values
        %   [riv,allocrange] = getAllowedRIV(obj,alloc) returns the set of allowed
        %   RIV values given information about the frequency hopping behaviour of
        %   the PSSCH resource allocation. These RIV values define contiguous 
        %   allocations, accounting for any gaps on the PSCCH resource pool.
        %   The allocrange output is a two column matrix which defines the VRB
        %   range associated with the RIV values. Each two element row of the
        %   matrix is the range pair, [Lcrb RBstart], for the associated RIV value.
        % 
        %   The alloc input parameter should be a structure including the 
        %   following fields: 
        %   FreqHopping         - Frequency hopping flag
        %                         (0 for non-hopping, 1 for hopping)
        %   Allocation          - Resource allocation parameter substructure
        %   The Allocation substructure should contain the following fields:
        %       Allocation.HoppingBits  - Hopping bits (0...2 bits). Only required for frequency hopping.
        %   
        %   Example:
        %   % Get the set of valid contiguous PSSCH allocations that can be used
        %   % from the configured PSSCH resource block pool.
        %
        %   period = PSCCHPeriod;
        %   sci.FreqHopping = 0;
        %   [riv,allocrange] = getAllowedRIV(period,sci);

            if isstruct(hoppingbits)
                % Extract the frequency hopping flag and hopping bits value
                freqhopping = hoppingbits.FreqHopping;
                if freqhopping
                    hoppingbits = hoppingbits.Allocation.HoppingBits;
                end
            else
                % If the parameter was not a structure then interpret it as the hopping bits value
                % and set frequency hopping on
                freqhopping = true;
            end        
            
            % Output variables
            riv = [];
            range = [];
            
            nslrb = obj.Config.NSLRB;
            % Use the terminology from TS 36.213 for the effective tx bandwidth for calculations
            nulrb = length(obj.PSSCHResourceBlockPool);
            
            % If non-hopping allocation
            if ~freqhopping  
                lcrbmax = nulrb;
                for crb = 1:lcrbmax                
                    % Set allocations within the subbands so that they don't cross them...                       
                    if nulrb == nslrb
                        nt_s1 = (0:nulrb-crb)';
                    else
                        % Since the pool is not the full bandwidth it may contain an intermediate gap
                        % that has to be accounted for in the RIV calculations
                        prbstart = obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_Start_r12;
                        prbnum = obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_Num_r12;
                        prbend = obj.Config.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12.prb_End_r12;
                        nt_s1 = [prbstart:prbstart+prbnum-crb, prbend-prbnum+1:prbend-crb+1]';
                        nt_s1 = unique(nt_s1);
                    end
                    rbstart = nt_s1;
                    % Concatenate this set of same length allocations to the output variables 
                    riv = [riv; calculateRIV(nslrb,crb,rbstart)];               %#ok<AGROW>
                    range = [range; [crb*ones(size(rbstart)) rbstart]];         %#ok<AGROW>       
                end
                return;  % RETURN... finished non-freq hopping 
            end
       
            % Hopping allocation
            n_ho = obj.Config.dataHoppingConfig_r12.rb_Offset_r12;
            nt_ho = n_ho + mod(n_ho,2);

            nhoppingbits = 1+(nslrb>=50);       % Note that we ould use also lteDCI/lteSCI too for these lengths formula
            
            y = ceil(log2(nslrb*(nslrb+1)/2)) - nhoppingbits;    
            lcrbmax = fix((2^y)/nslrb);   % The maximum length that can be signalled in the SCI

            % Establish whether it is type 1 or type 2 hopping
            type2 = (nhoppingbits==1 && hoppingbits==1) || hoppingbits == 3;

            if type2
                
                % Type 2 hopping
                nsb = obj.Config.dataHoppingConfig_r12.numSubbands_r12;
                hoppingbw = nulrb-nt_ho*(nsb>1);    % Hopping bandwidth 
                nsbrb = fix(hoppingbw/nsb);         % Number of RB in a subband
                lcrbmax = min(lcrbmax,nsbrb);       % Maximum number of RB in a contiguous allocation
                % Run through all the allocation lengths
                for crb = 1:lcrbmax
                    
                    % Calculate all possible starting RB indices for this allocation length 
                    % Set allocations within the subbands so that they don't cross them under hopping
                    %
                    % Add (column of) RB starting indices within a subband 
                    % to (a row of) the subband offsets, to give, 
                    % [0    nsbrb     2*nsbrb    ...]
                    % [1    nsbrb+1   2*nsbrb+1     ]
                    % [2    nsbrb+2   2*nsbrb+2     ]
                    % [...                          ]
                    nt_s1 = (0:nsbrb-crb)' + nsbrb*(0:nsb-1);
                    rbstart = nt_s1(:) + (nt_ho/2)*(nsb>1);
                    % Concatenate this set of same length allocations to the output variables
                    riv = [riv; calculateRIV(nslrb,crb,rbstart)];           %#ok<AGROW>
                    range = [range; [crb*ones(size(rbstart)) rbstart]];     %#ok<AGROW>
                end  
            else
                
                % Type 1 hopping           
                hoppingbw = nulrb-nt_ho-mod(nulrb,2);  % Hopping BW will be even
                if nhoppingbits == 1 || hoppingbits == 2
                    nsb = 2;    % 0 (1 bit) or 2 (2 bits)
                else
                    nsb = 4;    % 0 (2 bits) or 1 (2 bits)
                end       
                roffset = fix(hoppingbw/nsb);  % Length of implicit each hopping subband
                
                % Find the length of the contiguous part of the PSSCH resource block pool
                blklen = find(diff(obj.PSSCHResourceBlockPool)~=1);
                if isempty(blklen)
                    blklen = length(obj.PSSCHResourceBlockPool);
                end
                
                % Calculate hopping offset
                if bitand(hoppingbits,1)==0     % 0 or 2 - positive shift
                    r1 = roffset;               % Shift toward top of hopping band
                else                            % 1 - negative shift
                    r1 = -roffset;              % Shift back one subband 
                end
                
                % Run through all the allocation lengths       
                for crb = 1:min(blklen,lcrbmax)
                    
                    if blklen <  length(obj.PSSCHResourceBlockPool)
                        % Contiguous hopping bandwidth 'block' pair - 0 based VRB indices
                        blk1 = [nt_ho/2, blklen];
                        blk2 = [blklen, nt_ho/2+hoppingbw];
                    else
                        % Single hopping bandwidth 'block'
                        blk1 = [nt_ho/2, nt_ho/2+hoppingbw];
                        blk2 = [0, 0];
                    end

                    % Column of the starting indices of valid 'first slot' allocations 
                    rbstartA = [blk1(1):blk1(2)-crb, blk2(1):blk2(2)-crb]';

                    % Create two columns of [start,end] indices of valid 'first slot' allocations 
                    range1 = rbstartA + [0 crb];

                    % Create the hopped version of the above ranges
                    range2 = mod(range1 - nt_ho/2 + r1,hoppingbw) + nt_ho/2;

                    % Check whether the ranges are valid after hopping (create a logical vector associated with above ranges)
                    inc = (range2(:,2) > range2(:,1)) & ((range2(:,1) >= blk1(1) & range2(:,2) < blk1(2)) | (range2(:,1) >= blk2(1) & range2(:,2) < blk2(2)));
                    % Select valid 'first slot' starting indices
                    rbs = rbstartA(inc,1);                  
                   
                    % Concatenate this set of same length allocations to the output variables 
                    riv = [riv; calculateRIV(nslrb,crb,rbs)];           %#ok<AGROW>  
                    range = [range; [crb*ones(size(rbs)) rbs]];         %#ok<AGROW>                                             
                end
            end
                 
        end

        %% reset configuration
        function resetConfig(obj)
        % resetConfig Set configuration parameters to the default
        %   resetConfig(obj) sets the Config property to the default, which follows   
        %   TS 36.101 Section A.7.2, pool configuration #5-FDD (Table A.7.2.1-5).
        %
        %   Example:
        %   % This function is equivalent to calling,
        %   period.Config = PSCCHPeriod.defaultConfig(5,'5MHz');
        %   
        %   See also defaultConfig.

            obj.Config = obj.defaultConfig();

        end

        %% getPSCCHResources
        function [subframes,prb,n_PSCCH] = getPSCCHResources(obj,resparam)
        % getPSCCHResources Get PSCCH transmission resources for an nPSCCH value
        %   [SIND,PIND,NP] = getPSCCHResources(obj,npscch) returns indices for the
        %   pair of single PRB PSCCH transmissions derived from the scalar input,
        %   npscch ("Resource for PSCCH"). The pair of subframe indices, [a1 a2],
        %   are contained in the 1-by-2 vector SIND, and the pair of resource block
        %   indices, [b1 b2], are given in the 1-by-2 vector PIND. The ordering 
        %   of the indices follows that of TS 36.213 Section 14.2.1.1 and, as 
        %   such, may not be monotonically increasing. The indices are 0-based and
        %   for the subframes they are relative to the start of the PSCCH period.
        %   If the npscch input is empty then a resource value is chosen
        %   uniformly at random from the allowed value set. This follows the
        %   behaviour required for PSCCH resource selection in transmission mode 2.
        %   The value used is always returned in NP output.
        %   
        %   [SIND,PIND,NP] = getPSCCHResources(obj,dci) is the same as the above
        %   except the PSCCH resource value is given by the following field in 
        %   the dci input structure:
        %   PSCCHResource - Resource for PSCCH
        %    
        %   As with the first signature, if this field value is empty then an 
        %   allowed resource value is selected uniformly at random. In the case
        %   of transmission mode 1, these parameters will have been sent to the UE
        %   via a SCI format 5 message.
        %
        %   Example: 
        %   % Get the PSCCH resources associated with the last allowed resource for
        %   % PSCCH index value.
        %   period = PSCCHPeriod;
        %   period.PSCCHSubframePool
        %   period.PSCCHResourceBlockPool
        %   dci.PSCCHResource = period.NumPSCCHResource - 1;
        %   [sfindices,prbindices,npssch] = period.getPSCCHResources(dci);
        %   
        %   See also getPSSCHResources, displayPeriod.

            % Prepare parameters
            if isstruct(resparam)
                n_PSCCH = resparam.PSCCHResource;
            else
                n_PSCCH = resparam;
            end
                
            % If "resource for PSCCH" wasn't presented then randomly select an allowed index with a random distribution
            % otherwise ensure that the presented index is valid
            % TS 36.321 Section 5.14.1.1 - SL Grant reception and SCI transmission
            if isempty(n_PSCCH)
                n_PSCCH = randi(obj.NumPSCCHResource)-1;
            else
                n_PSCCH = mod(n_PSCCH,obj.NumPSCCHResource);
            end
            
            % Calculate subframe and resource pool indices for the two transmission RBs
            a1 = fix(n_PSCCH/obj.L_PSCCH);                             % RB index A
            b1 = mod(n_PSCCH,obj.L_PSCCH);                             % SF index A
            a2 = a1 + fix(obj.M_RB_PSCCH_RP/2);                        % RB index B
            b2 = mod( n_PSCCH+1+mod(a1,obj.L_PSCCH-1), obj.L_PSCCH);   % SF index B

            % Look-up the actual subframe and resource indices from the pools for the two transmission RBs
            m1 = obj.PSCCHResourceBlockPool(a1+1);
            l1 = obj.PSCCHSubframePool(b1+1);
            m2 = obj.PSCCHResourceBlockPool(a2+1);
            l2 = obj.PSCCHSubframePool(b2+1);
            
            % Output the PRB and subframe indices for the pair of PSCCH
            prb = [m1 m2];
            subframes  = [l1 l2];
            
        end            
        
        %% getPSSCHResources
        function [subframes,prb,poolindices] = getPSSCHResources(obj,resparam)
        % getPSSCHResources Get resource indices for PSSCH transmission
        %   [SIND,PIND,NPSF] = getPSSCHResources(obj,SCISTR) returns the resource 
        %   indices for the sequence of possible PSSCH transmissions associated 
        %   with the resource parts of an SCI format 0 message. The subframe 
        %   indices are returned in the 1-by-NPSSCH vector SIND. The physical 
        %   resource block indices are given in the NPRB-by-NPSSCH matrix PIND. 
        %   Each column of the matrix contains the PRB allocation indices for the
        %   associated subframe. These indices are 0-based and for the subframes
        %   they are relative to the start of the PSCCH period. The 1-by-NPSSCH 
        %   vector NPSF contains the 0-based indices of the used subframes relative
        %   to their position in the subframe pool. The number of subframes NPSCCH
        %   will always be a multiple of four (in order that four RV transmissions
        %   can be made for each available transport block), however it could be
        %   zero to satisfy this constraint. 
        %   
        %   The input SCISTR must be a structure containing the fields:
        %   TimeResourcePattern - Time resource pattern index (I_TRP)
        %   FreqHopping         - Frequency hopping flag
        %                         (0 for non-hopping, 1 for hopping)
        %   Allocation          - Resource allocation parameter substructure
        %   The Allocation substructure should contain the following fields:
        %       Allocation.HoppingBits  - Hopping bits (0...3, depending on NSLRB)
        %                                 Only required for frequency hopping
        %       Allocation.RIV          - Resource indication value
        %
        %   These parameters represent the resource related fields which are 
        %   part of the associated SCI format 0 carried on the PSCCH. In the case
        %   of transmission mode 1, these parameters will have been sent to the UE
        %   via an SCI format 5 message. The TimeResourcePattern index parameter
        %   controls the subframes selected from the PSSCH pool, while the other 
        %   parameters relate to the PRB selected from the PSSCH physical resource 
        %   block pool.
        %   
        %   Example:
        %   % Get the PSSCH resource information for the specified transmission 
        %   % sequence (single PRB, type 1 frequency hopping)
        %
        %   period = PSCCHPeriod;
        %   sci.TimeResourcePattern = 106;
        %   sci.FreqHopping = 1;
        %   sci.Allocation.HoppingBits = 0;
        %   sci.Allocation.RIV = 0;
        %   [sfi,prbi,sfpi] = getPSSCHResources(period,sci)
        %
        %   See also AllowedITRP, getAllowedRIV, lteSCIResourceAllocation.

            % Prepare parameters
            if isstruct(resparam)
                sci = resparam;
                itrp = sci.TimeResourcePattern;
            else
                sci = [];
                itrp = resparam;
            end
                     
            % Get extended bitmap associated with I_TRP
            bitmap = determineSubframeIndicatorBitmap(obj.Config,itrp,length(obj.PSSCHSubframePool));
            
            % Use extended bitmap to select from the PSSCH subframe pool 
            subframes = obj.PSSCHSubframePool(bitmap);
            poolindices = find(bitmap)-1;
            
            % Fix the number of subframes to be a multiple of 4
            tl = 4*fix(length(subframes)/4);    
            subframes = subframes(1:tl);
            poolindices = poolindices(1:tl);
                
            if isempty(sci)
                % If the SCI parameter was not explicitly defined then use the entire RB pool
                prb = repmat(obj.PSSCHResourceBlockPool,1,length(subframes));
            else
                % Set up the parameter set required by lteSCIResourceAllocation
                ue.NSLRB = obj.Config.NSLRB;
                ue.PSSCHHoppingParameter = obj.Config.dataHoppingConfig_r12.hoppingParameter_r12;
                ue.NSubbands             = obj.Config.dataHoppingConfig_r12.numSubbands_r12;
                ue.PSSCHHoppingOffset    = obj.Config.dataHoppingConfig_r12.rb_Offset_r12;
                ue.PRBPool               = obj.PSSCHResourceBlockPool;
                % Iterate through the active subframes and concatenate the allocated PRB as columns
                % in the output matrix
                prb = [];
                sci.SCIFormat = 'Format0';
                for i = poolindices
                    ue.NSubframePSSCH = i;
                    prb = [prb lteSCIResourceAllocation(ue,sci)]; %#ok<AGROW>
                end
            
            end
        end
    
    end
    
    methods (Access = private)
        
        %% Object initialization function
        function obj = initialize(obj)
                             
            % Refer to the entire Configuration parameter structure
            params = obj.Config;
                               
            % Calculate first subframe in this PSCCH period
            [jbegin,jend] = determinejbegin(params);
            obj.NSubframeBegin = jbegin;
            obj.PeriodLength = jend - jbegin + 1;
                        
            % Determine the sync subframes (if any) 
            obj.SyncSubframes = determineSyncSubframes(params);
            
            % Determine the subframe and resource block pools for PSCCH       
            configC = params.sc_TF_ResourceConfig_r12;
            [obj.PSCCHSubframePool,obj.PSCCHResourceBlockPool] = determinePSCCHSubframeAndResourceAllocation(params,configC);
            calcPSCCHResourceInfo(obj);  % Call after the PSCCH pools have been set
            
            % Set the transmission mode property
            modenames = {'Mode1','Mode2'};
            mode2 = isfield(params,'ue_SelectedResourceConfig_r12') && (~isfield(params,'UESelected') || strcmpi(params.UESelected,'On'));
            obj.TxMode = modenames{mode2+1};
             
            % Determine the subframe and resource block pools for PSSCH
            configS = [];
            trptsubset = [];          
            if mode2
                configS = params.ue_SelectedResourceConfig_r12.data_TF_ResourceConfig_r12;
                trptsubset = params.ue_SelectedResourceConfig_r12.trpt_Subset_r12;
            end
            [obj.PSSCHSubframePool,obj.PSSCHResourceBlockPool] = determinePSSCHSubframeAndResourceAllocation(params,configS);
            obj.AllowedITRP = determineAllowedITPRMode2(params,trptsubset);
                   
            % Make all subframe pool indices relative to PSCCH period start
            obj.SyncSubframes = obj.SyncSubframes-jbegin;
            obj.PSCCHSubframePool = obj.PSCCHSubframePool-jbegin;
            obj.PSSCHSubframePool = obj.PSSCHSubframePool-jbegin;
            
            %  Validate pools
            if max(obj.PSCCHResourceBlockPool) >= params.NSLRB
                error('The PSCCH resource block pool indices exceed the number of sidelink resource blocks.');
            end
            if max(obj.PSSCHResourceBlockPool) >= params.NSLRB
                error('The PSSCH resource block pool indices exceed the number of sidelink resource blocks.');
            end
            if max(obj.PSCCHSubframePool) >= params.sc_Period_r12
                error('The PSCCH subframe pool indices exceed the PSCCH period length.');
            end
            if max(obj.PSSCHSubframePool) >= params.sc_Period_r12
                error('The PSSCH subframe pool indices exceed the PSCCH period length.');
            end 
            
            % There are further restrictions placed on the RRC message parameters
            % which are not enforced, for example,
            %
            % TS 36.331 Section 6.3.8 - SL-CommResourcePool field descriptions
            %
            % PSCCH period length restrictions...
            % If FDD or TDD and configuration 1...5 then {40,80,160,320} 
            % If TDD and configuration 0 then {70,140,280}
            % If TDD and configuration 6 then {60,120,240}
            %
            % SL-TF-ResourceConfig subframeBitmap restrictions...
            % Indicates the subframe bitmap indicating resources used for sidelink. E-UTRAN configures,
            % value bs40 for FDD,
            % and the following values for TDD:
            % value bs42 for configuration0,
            % value bs16 for configuration1, 
            % value bs8 for configuration2,
            % value bs12 for configuration3,
            % value bs8 for configuration4,
            % value bs4 for configuration5,
            % value bs30 for configuration6.
            
        end
        
        %% Initialize PSCCH resource related object properties
        function calcPSCCHResourceInfo(obj)
            
            % Intermediate variables
            obj.L_PSCCH = length(obj.PSCCHSubframePool);                % Number of subframes in resource pool
            obj.M_RB_PSCCH_RP = length(obj.PSCCHResourceBlockPool);     % Number of RB in resource pool

            % Suppose the UE is selecting the resources (mode 2) then we need to known the range of nPSCCH ('Resource for PSCCH')
            % The allowed values of resource for PSCCH are in the range between [0 NumPSCCHResource)
            obj.NumPSCCHResource =  fix(obj.M_RB_PSCCH_RP/2)*obj.L_PSCCH; % Exclusive upper bound
            
        end 
        
    end
      
end


%% TS 34.213 Section 8.1.1 RIV calculation (Uplink resource allocation type 0)
function riv = calculateRIV(nulrb,lcrb,rbstart)

    if lcrb-1 <= fix(nulrb/2)
        riv = nulrb*(lcrb-1) + rbstart;
    else
        riv = nulrb*(nulrb-lcrb+1) + (nulrb-1-rbstart);
    end

end        

%% Subframe range indices (relative to SFN/DFN) of current PSCCH period 
function [jbegin,jend] = determinejbegin(params)

    jbegin = params.sc_TF_ResourceConfig_r12.offsetIndicator_r12 + params.NPSCCHPeriod*params.sc_Period_r12;
    jend = jbegin + params.sc_Period_r12 - 1;

end

%% PSSS/SSSS/PSBCH: Sidelink SYNC
% Procedures for determining sync subframes
function subframes = determineSyncSubframes(params)
    
    % First SFN of period 
    [jbegin,jend] = determinejbegin(params);

    % Offset of first sync subframe in PSCCH period
    syncoffset = mod(params.syncConfig.syncOffsetIndicator_r12 - jbegin,40);

    % Calculate the active sync subframes in the period
    if isfield(params,'syncConfig') && (~isfield(params,'SyncEnable') || strcmpi(params.SyncEnable,'On'))
             subframes =  jbegin+syncoffset : 40 : jend;
        else
             subframes = [];
    end

end 
       
%% Determine the uplink subframe indices starting at 'first', where the 
% 'range' parameter is optionally the number of indices to find or the 
% last indices in a range of interest
function uplinksfindices = determineUplinkSubframes(config,first,range,opts)

    % inclusive first, exclusive last 
    if nargin < 4
        opts = [];
    end
    
    % Calculate uplink subframe pattern from 10 subframes from 'first' 
    % subframe (inclusive). Ensure that UL CP is set for the duplex info
    % code to work (the value doesn't matter)    
    config.CyclicPrefixUL = 'Normal';   % Present to prioritize the uplink 
    config.CyclicPrefix = 'Normal';     % Always optional parameter for lteDuplexInfo so add it to remove any warnings
    config.SSC = 0;
    uplinksf = arrayfun(@(x)strcmp(getfield(lteDuplexingInfo(setfield(config,'NSubframe',x)),'SubframeType'),'Uplink'),first:first+9);    %#ok<SFLD> % Get the uplink/downlink pattern of the first frame of the period
    
    if any(strcmpi(opts,'first'))
       % If 'first' then identify first 'range' uplink subframes after 'first' subframe
       nuplinksf = sum(uplinksf);      % Number of uplink subframes in a frame
       uplinksf = repmat(uplinksf,1,ceil(range/nuplinksf));
       uplinksfindices = find(uplinksf,range);
    else
       % Otherwise identify uplink subframes between 'first' and 'range' subframes
       nsf = max(range-first,0);
       nuplinksf = length(uplinksf);   % Number of subframes in a frame
       uplinksf = repmat(uplinksf,1,ceil(nsf/nuplinksf));
       uplinksfindices = find(uplinksf(1:nsf));
    end
    % Adjust relative to the 'first' index (and '1based' if required, '0based' is the default)
    offset = ~any(strcmpi(opts,'1based'));
    uplinksfindices = uplinksfindices+first-offset; 
    
end

%% Determine RB pool from the RB transport format parameters (PSCCH and PSSCH mode 2)
function rbpool = determineResourcePool(config)
    rbpool = ...
        [ config.prb_Start_r12:(config.prb_Start_r12+config.prb_Num_r12-1), ...
          config.prb_End_r12-config.prb_Num_r12+1:config.prb_End_r12 ]';
    rbpool = unique(rbpool);
end

%% Procedures for determining PSSCH subframe pool and resource block pool
function [subframepool,rbpool] = determinePSSCHSubframeAndResourceAllocation(params,configS)

    % PSCCH TF part
    configC = params.sc_TF_ResourceConfig_r12;
    
    [jbegin,jend] = determinejbegin(params);
        
    if nargin == 1 || isempty(configS)

        % mode 1 - E-UTRAN directed (not constrained)
        % We only need the last (uplink) subframe of the PSCCH subframe pool, 
        % and the PSSCH subframe pool will be all the uplink subframes from 
        % the end of the PSCCH pool to the end of the PSCCH period

        subframepoolPSCCH = determinePSCCHSubframeAndResourceAllocation(params,configC);

        subframepool = determineUplinkSubframes(params,...
                                                subframepoolPSCCH(end)+1,...  % Start, for mode 1 the PSSCH subframe pool starts after the end of the PSCCH pool (the +1)
                                                jend+1,'last');               % End, and ends at the end of the PSCCH period    

        % Resource block pool is the full bandwidth for PSSCH mode 1 
        rbpool = (0:params.NSLRB-1)';
    else

        % mode 2 - UE selected, but heavily constrained by protocol
        % Uses data-TF-ResourceConfig-r12, a separate instance of SL-TF-ResourceConfig-r12

        % RB pool
        rbpool = determineResourcePool(configS);

        % Subframes within the PSSCH part of the PSCCH period    
        uplinksfindices = determineUplinkSubframes(params,...
                                                   jbegin+configS.offsetIndicator_r12,... % For mode 2 the PSSCH subframe pool starts after the PSSCH resource configuration offset parameter added to jbegin
                                                   jend+1,'last');                        % And ends at the end of the PSCCH period    

        % Then use bitmap (extended) to select from the above uplink subframes 
        % to create the actual pool
        bitmap = logicalBitmap(configS.subframeBitmap_r12);
        extendedbitmap = repmat(bitmap,1,ceil(length(uplinksfindices)/length(bitmap)));
        subframepool = uplinksfindices(extendedbitmap(1:length(uplinksfindices)));
    end

end

%% Procedures for determining PSCCH subframe pool and resource block pool
function [subframepool,rbpool] = determinePSCCHSubframeAndResourceAllocation(params,configC)

    % The subframe period is aligned relative to the SFN/DFN 
    % so the offset just shifts the start of the PSCCH period relative to SFN/DFN #0.

    % Calculate for a single PSCCH period (the ith period) 
    jbegin = determinejbegin(params);               % jbegin is the subframe index (relative to SFN/DFN #0) at the beginning of the ith PSCCH period
    Nd = length(configC.subframeBitmap_r12);        % N' is the length of the bitmap

    % Need the first N' uplink subframes
    uplinksfindices = determineUplinkSubframes(params,jbegin,Nd,{'first','0based'});
    
    % Use the bitmap to index and select from the set of N' uplink indices
    subframepool = uplinksfindices(logicalBitmap(configC.subframeBitmap_r12));

    % RB pool from SL-TF-ResourceConfig-r12
    rbpool = determineResourcePool(configC);
  
end

%% PSCCH: Sidelink CONTROL
% Procedures for determining PSCCH subframe and resource allocation
function output = determineSubframeIndicatorBitmap(config,itrp,bitmaplength)
    
    % One-time construction of the tables
    persistent table_NTRP8;
    persistent table_NTRP7;
    persistent table_NTRP6;
    
    if isempty(table_NTRP8)  % One time initialization of cell tables

    % ----------------------------
    % Table 14.1.1.1.1-1: Time Resource pattern Index mapping for N_TRP = 8
    table_NTRP8 = ...
    {    
    0	1	[1,0,0,0,0,0,0,0];
    1	1	[0,1,0,0,0,0,0,0];
    2	1	[0,0,1,0,0,0,0,0];
    3	1	[0,0,0,1,0,0,0,0];
    4	1	[0,0,0,0,1,0,0,0];
    5	1	[0,0,0,0,0,1,0,0];
    6	1	[0,0,0,0,0,0,1,0];
    7	1	[0,0,0,0,0,0,0,1];
    8	2	[1,1,0,0,0,0,0,0];
    9	2	[1,0,1,0,0,0,0,0];
    10	2	[0,1,1,0,0,0,0,0];
    11	2	[1,0,0,1,0,0,0,0];
    12	2	[0,1,0,1,0,0,0,0];
    13	2	[0,0,1,1,0,0,0,0];
    14	2	[1,0,0,0,1,0,0,0];
    15	2	[0,1,0,0,1,0,0,0];
    16	2	[0,0,1,0,1,0,0,0];
    17	2	[0,0,0,1,1,0,0,0];
    18	2	[1,0,0,0,0,1,0,0];
    19	2	[0,1,0,0,0,1,0,0];
    20	2	[0,0,1,0,0,1,0,0];
    21	2	[0,0,0,1,0,1,0,0];
    22	2	[0,0,0,0,1,1,0,0];
    23	2	[1,0,0,0,0,0,1,0];
    24	2	[0,1,0,0,0,0,1,0];
    25	2	[0,0,1,0,0,0,1,0];
    26	2	[0,0,0,1,0,0,1,0];
    27	2	[0,0,0,0,1,0,1,0];
    28	2	[0,0,0,0,0,1,1,0];
    29	2	[1,0,0,0,0,0,0,1];
    30	2	[0,1,0,0,0,0,0,1];
    31	2	[0,0,1,0,0,0,0,1];
    32	2	[0,0,0,1,0,0,0,1];
    33	2	[0,0,0,0,1,0,0,1];
    34	2	[0,0,0,0,0,1,0,1];
    35	2	[0,0,0,0,0,0,1,1];
    36	4	[1,1,1,1,0,0,0,0];
    37	4	[1,1,1,0,1,0,0,0];
    38	4	[1,1,0,1,1,0,0,0];
    39	4	[1,0,1,1,1,0,0,0];
    40	4	[0,1,1,1,1,0,0,0];
    41	4	[1,1,1,0,0,1,0,0];
    42	4	[1,1,0,1,0,1,0,0];
    43	4	[1,0,1,1,0,1,0,0];
    44	4	[0,1,1,1,0,1,0,0];
    45	4	[1,1,0,0,1,1,0,0];
    46	4	[1,0,1,0,1,1,0,0];
    47	4	[0,1,1,0,1,1,0,0];
    48	4	[1,0,0,1,1,1,0,0];
    49	4	[0,1,0,1,1,1,0,0];
    50	4	[0,0,1,1,1,1,0,0];
    51	4	[1,1,1,0,0,0,1,0];
    52	4	[1,1,0,1,0,0,1,0];
    53	4	[1,0,1,1,0,0,1,0];
    54	4	[0,1,1,1,0,0,1,0];
    55	4	[1,1,0,0,1,0,1,0];
    56	4	[1,0,1,0,1,0,1,0];
    57	4	[0,1,1,0,1,0,1,0];
    58	4	[1,0,0,1,1,0,1,0];
    59	4	[0,1,0,1,1,0,1,0];
    60	4	[0,0,1,1,1,0,1,0];
    61	4	[1,1,0,0,0,1,1,0];
    62	4	[1,0,1,0,0,1,1,0];
    63	4	[0,1,1,0,0,1,1,0];
    64	4	[1,0,0,1,0,1,1,0];
    65	4	[0,1,0,1,0,1,1,0];
    66	4	[0,0,1,1,0,1,1,0];
    67	4	[1,0,0,0,1,1,1,0];
    68	4	[0,1,0,0,1,1,1,0];
    69	4	[0,0,1,0,1,1,1,0];
    70	4	[0,0,0,1,1,1,1,0];
    71	4	[1,1,1,0,0,0,0,1];
    72	4	[1,1,0,1,0,0,0,1];
    73	4	[1,0,1,1,0,0,0,1];
    74	4	[0,1,1,1,0,0,0,1];
    75	4	[1,1,0,0,1,0,0,1];
    76	4	[1,0,1,0,1,0,0,1];
    77	4	[0,1,1,0,1,0,0,1];
    78	4	[1,0,0,1,1,0,0,1];
    79	4	[0,1,0,1,1,0,0,1];
    80	4	[0,0,1,1,1,0,0,1];
    81	4	[1,1,0,0,0,1,0,1];
    82	4	[1,0,1,0,0,1,0,1];
    83	4	[0,1,1,0,0,1,0,1];
    84	4	[1,0,0,1,0,1,0,1];
    85	4	[0,1,0,1,0,1,0,1];
    86	4	[0,0,1,1,0,1,0,1];
    87	4	[1,0,0,0,1,1,0,1];
    88	4	[0,1,0,0,1,1,0,1];
    89	4	[0,0,1,0,1,1,0,1];
    90	4	[0,0,0,1,1,1,0,1];
    91	4	[1,1,0,0,0,0,1,1];
    92	4	[1,0,1,0,0,0,1,1];
    93	4	[0,1,1,0,0,0,1,1];
    94	4	[1,0,0,1,0,0,1,1];
    95	4	[0,1,0,1,0,0,1,1];
    96	4	[0,0,1,1,0,0,1,1];
    97	4	[1,0,0,0,1,0,1,1];
    98	4	[0,1,0,0,1,0,1,1];
    99	4	[0,0,1,0,1,0,1,1];
    100	4	[0,0,0,1,1,0,1,1];
    101	4	[1,0,0,0,0,1,1,1];
    102	4	[0,1,0,0,0,1,1,1];
    103	4	[0,0,1,0,0,1,1,1];
    104	4	[0,0,0,1,0,1,1,1];
    105	4	[0,0,0,0,1,1,1,1];
    106	8	[1,1,1,1,1,1,1,1];
    % 107-127	 reserved	reserved
    };

    % ----------------------------
    % Table 14.1.1.1.1-2: Time Resource pattern Index mapping for N_TRP=7
    table_NTRP7 = ...
    {  	 	 
    % 0	reserved	reserved
    0	0	[0,0,0,0,0,0,0];    % Reserved
    1	1	[1,0,0,0,0,0,0];
    2	1	[0,1,0,0,0,0,0];
    3	2	[1,1,0,0,0,0,0];
    4	1	[0,0,1,0,0,0,0];
    5	2	[1,0,1,0,0,0,0];
    6	2	[0,1,1,0,0,0,0];
    7	3	[1,1,1,0,0,0,0];
    8	1	[0,0,0,1,0,0,0];
    9	2	[1,0,0,1,0,0,0];
    10	2	[0,1,0,1,0,0,0];
    11	3	[1,1,0,1,0,0,0];
    12	2	[0,0,1,1,0,0,0];
    13	3	[1,0,1,1,0,0,0];
    14	3	[0,1,1,1,0,0,0];
    15	4	[1,1,1,1,0,0,0];
    16	1	[0,0,0,0,1,0,0];
    17	2	[1,0,0,0,1,0,0];
    18	2	[0,1,0,0,1,0,0];
    19	3	[1,1,0,0,1,0,0];
    20	2	[0,0,1,0,1,0,0];
    21	3	[1,0,1,0,1,0,0];
    22	3	[0,1,1,0,1,0,0];
    23	4	[1,1,1,0,1,0,0];
    24	2	[0,0,0,1,1,0,0];
    25	3	[1,0,0,1,1,0,0];
    26	3	[0,1,0,1,1,0,0];
    27	4	[1,1,0,1,1,0,0];
    28	3	[0,0,1,1,1,0,0];
    29	4	[1,0,1,1,1,0,0];
    30	4	[0,1,1,1,1,0,0];
    31	5	[1,1,1,1,1,0,0];
    32	1	[0,0,0,0,0,1,0];
    33	2	[1,0,0,0,0,1,0];
    34	2	[0,1,0,0,0,1,0];
    35	3	[1,1,0,0,0,1,0];
    36	2	[0,0,1,0,0,1,0];
    37	3	[1,0,1,0,0,1,0];
    38	3	[0,1,1,0,0,1,0];
    39	4	[1,1,1,0,0,1,0];
    40	2	[0,0,0,1,0,1,0];
    41	3	[1,0,0,1,0,1,0];
    42	3	[0,1,0,1,0,1,0];
    43	4	[1,1,0,1,0,1,0];
    44	3	[0,0,1,1,0,1,0];
    45	4	[1,0,1,1,0,1,0];
    46	4	[0,1,1,1,0,1,0];
    47	5	[1,1,1,1,0,1,0];
    48	2	[0,0,0,0,1,1,0];
    49	3	[1,0,0,0,1,1,0];
    50	3	[0,1,0,0,1,1,0];
    51	4	[1,1,0,0,1,1,0];
    52	3	[0,0,1,0,1,1,0];
    53	4	[1,0,1,0,1,1,0];
    54	4	[0,1,1,0,1,1,0];
    55	5	[1,1,1,0,1,1,0];
    56	3	[0,0,0,1,1,1,0];
    57	4	[1,0,0,1,1,1,0];
    58	4	[0,1,0,1,1,1,0];
    59	5	[1,1,0,1,1,1,0];
    60	4	[0,0,1,1,1,1,0];
    61	5	[1,0,1,1,1,1,0];
    62	5	[0,1,1,1,1,1,0];
    63	6	[1,1,1,1,1,1,0];
    64	1	[0,0,0,0,0,0,1];
    65	2	[1,0,0,0,0,0,1];
    66	2	[0,1,0,0,0,0,1];
    67	3	[1,1,0,0,0,0,1];
    68	2	[0,0,1,0,0,0,1];
    69	3	[1,0,1,0,0,0,1];
    70	3	[0,1,1,0,0,0,1];
    71	4	[1,1,1,0,0,0,1];
    72	2	[0,0,0,1,0,0,1];
    73	3	[1,0,0,1,0,0,1];
    74	3	[0,1,0,1,0,0,1];
    75	4	[1,1,0,1,0,0,1];
    76	3	[0,0,1,1,0,0,1];
    77	4	[1,0,1,1,0,0,1];
    78	4	[0,1,1,1,0,0,1];
    79	5	[1,1,1,1,0,0,1];
    80	2	[0,0,0,0,1,0,1];
    81	3	[1,0,0,0,1,0,1];
    82	3	[0,1,0,0,1,0,1];
    83	4	[1,1,0,0,1,0,1];
    84	3	[0,0,1,0,1,0,1];
    85	4	[1,0,1,0,1,0,1];
    86	4	[0,1,1,0,1,0,1];
    87	5	[1,1,1,0,1,0,1];
    88	3	[0,0,0,1,1,0,1];
    89	4	[1,0,0,1,1,0,1];
    90	4	[0,1,0,1,1,0,1];
    91	5	[1,1,0,1,1,0,1];
    92	4	[0,0,1,1,1,0,1];
    93	5	[1,0,1,1,1,0,1];
    94	5	[0,1,1,1,1,0,1];
    95	6	[1,1,1,1,1,0,1];
    96	2	[0,0,0,0,0,1,1];
    97	3	[1,0,0,0,0,1,1];
    98	3	[0,1,0,0,0,1,1];
    99	4	[1,1,0,0,0,1,1];
    100	3	[0,0,1,0,0,1,1];
    101	4	[1,0,1,0,0,1,1];
    102	4	[0,1,1,0,0,1,1];
    103	5	[1,1,1,0,0,1,1];
    104	3	[0,0,0,1,0,1,1];
    105	4	[1,0,0,1,0,1,1];
    106	4	[0,1,0,1,0,1,1];
    107	5	[1,1,0,1,0,1,1];
    108	4	[0,0,1,1,0,1,1];
    109	5	[1,0,1,1,0,1,1];
    110	5	[0,1,1,1,0,1,1];
    111	6	[1,1,1,1,0,1,1];
    112	3	[0,0,0,0,1,1,1];
    113	4	[1,0,0,0,1,1,1];
    114	4	[0,1,0,0,1,1,1];
    115	5	[1,1,0,0,1,1,1];
    116	4	[0,0,1,0,1,1,1];
    117	5	[1,0,1,0,1,1,1];
    118	5	[0,1,1,0,1,1,1];
    119	6	[1,1,1,0,1,1,1];
    120	4	[0,0,0,1,1,1,1];
    121	5	[1,0,0,1,1,1,1];
    122	5	[0,1,0,1,1,1,1];
    123	6	[1,1,0,1,1,1,1];
    124	5	[0,0,1,1,1,1,1];
    125	6	[1,0,1,1,1,1,1];
    126	6	[0,1,1,1,1,1,1];
    127	7	[1,1,1,1,1,1,1];
    };

    % ---------------------------------------
    % Table 14.1.1.1.1-3: Time Resource pattern Index mapping for N_TRP = 6
    table_NTRP6 = ...
    {  	 	 	 
    % 0	reserved	reserved
    0	0	[0,0,0,0,0,0];      % Reserved
    1	1	[1,0,0,0,0,0];	
    2	1	[0,1,0,0,0,0];
    3	2	[1,1,0,0,0,0];	
    4	1	[0,0,1,0,0,0];	
    5	2	[1,0,1,0,0,0];	
    6	2	[0,1,1,0,0,0];	
    7	3	[1,1,1,0,0,0];	
    8	1	[0,0,0,1,0,0];	
    9	2	[1,0,0,1,0,0];	
    10	2	[0,1,0,1,0,0];	
    11	3	[1,1,0,1,0,0];	
    12	2	[0,0,1,1,0,0];
    13	3	[1,0,1,1,0,0];	
    14	3	[0,1,1,1,0,0];	
    15	4	[1,1,1,1,0,0];	
    16	1	[0,0,0,0,1,0];	
    17	2	[1,0,0,0,1,0];	
    18	2	[0,1,0,0,1,0];	
    19	3	[1,1,0,0,1,0];	
    20	2	[0,0,1,0,1,0];	
    21	3	[1,0,1,0,1,0];	
    22	3	[0,1,1,0,1,0];
    23	4	[1,1,1,0,1,0];
    24	2	[0,0,0,1,1,0];
    25	3	[1,0,0,1,1,0];
    26	3	[0,1,0,1,1,0];
    27	4	[1,1,0,1,1,0];
    28	3	[0,0,1,1,1,0];
    29	4	[1,0,1,1,1,0];
    30	4	[0,1,1,1,1,0];
    31	5	[1,1,1,1,1,0];
    32	1	[0,0,0,0,0,1];
    33	2	[1,0,0,0,0,1];
    34	2	[0,1,0,0,0,1];
    35	3	[1,1,0,0,0,1];
    36	2	[0,0,1,0,0,1];
    37	3	[1,0,1,0,0,1];
    38	3	[0,1,1,0,0,1];
    39	4	[1,1,1,0,0,1];
    40	2	[0,0,0,1,0,1];
    41	3	[1,0,0,1,0,1];
    42	3	[0,1,0,1,0,1];
    43	4	[1,1,0,1,0,1];
    44	3	[0,0,1,1,0,1];
    45	4	[1,0,1,1,0,1];
    46	4	[0,1,1,1,0,1];
    47	5	[1,1,1,1,0,1];
    48	2	[0,0,0,0,1,1];
    49	3	[1,0,0,0,1,1];
    50	3	[0,1,0,0,1,1];
    51	4	[1,1,0,0,1,1];
    52	3	[0,0,1,0,1,1];
    53	4	[1,0,1,0,1,1];
    54	4	[0,1,1,0,1,1];
    55	5	[1,1,1,0,1,1];
    56	3	[0,0,0,1,1,1];
    57	4	[1,0,0,1,1,1];
    58	4	[0,1,0,1,1,1];
    59	5	[1,1,0,1,1,1];
    60	4	[0,0,1,1,1,1];
    61	5	[1,0,1,1,1,1];
    62	5	[0,1,1,1,1,1];
    63	6	[1,1,1,1,1,1];
    % 64-127	 reserved	 reserved
    };
    end

    % Select the active table
    if strcmpi(config.DuplexMode,'FDD') || any(config.TDDConfig==[1 2 4 5])
        table = table_NTRP8;
    elseif config.TDDConfig==0
        table = table_NTRP7;
    else  
        table = table_NTRP6;    % TDDConfig ==[3 6]
    end    
    
    % If no I_TRP was presented then return the entire kTRP column
    if nargin == 1
        output = cell2mat(table(:,2));  % Second column of active table
        return;
    end

    % Look up the bitmap associated with the input index (I_TRP)
    output = logical(table{itrp+1,3});

    % Extend this bitmap if required
    if nargin > 2
       output = repmat(output,1,ceil(bitmaplength/length(output)));
       output = output(1:bitmaplength);
    end

end 

%% PSCCH: Sidelink CONTROL
% Procedures for determining PSCCH subframe and resource allocation
function [indexset,allowedktrpset] = determineAllowedITPRMode2(config,trpt_subset_r12)
    
    % If no trpt_subset_r12 presented then return all valid kTRP
    if nargin==1 || isempty(trpt_subset_r12)
        ktrplist = determineSubframeIndicatorBitmap(config);
        indexset = find(ktrplist~=0)-1;
        allowedktrpset = [];
        return
    end
    
    % Ensure that the trp subset bitmap is a logical array 
    trpt_subset_r12 = logicalBitmap(trpt_subset_r12);
      
    % From the duplex mode and tdd config get the subset of k_trp
    if strcmpi(config.DuplexMode,'FDD') || any(config.TDDConfig==[1 2 4 5])
        allowedktrpset = [1 2 4];
    elseif config.TDDConfig==0
        allowedktrpset = 1:5;         % TDDConfig == 0
    else  
        allowedktrpset = 1:4;         % TDDConfig == [3 6]
    end    

    % trpt_subset_r12 - optional in the RRC messaging
    % This subset parameter is not applicable if TDD and TDDConfig = 5
    if nargin > 1 && ~isempty(trpt_subset_r12) && ~(strcmpi(config.DuplexMode,'TDD') && any(config.TDDConfig==5))
        allowedktrpset = allowedktrpset(logical(trpt_subset_r12(1:min(length(trpt_subset_r12),length(allowedktrpset)))));
    end

    % From set of k_trp, get list of associated i_trp
    % Get ktrplist for this config...
    ktrplist = determineSubframeIndicatorBitmap(config);
    indexset = find(ismember(ktrplist,allowedktrpset))-1;

    % In mode 2 the ue can then select from this list at random

    % Note that other useful dimensions wrt the calculations include:
    % The set of kTRP that are allowed given the duplexing information
    % XTRP (the length of the set) so that the bitmap length can be 
    % established/validated, or at least the number of used bits established
    
end

%% Convert a numerical or 'bitmap' character vector into a logical one   
function bitmap = logicalBitmap(bitmap)

    if ischar(bitmap)
        bitmap = bitmap~='0';
    else
        bitmap = logical(bitmap);
    end
    
end
