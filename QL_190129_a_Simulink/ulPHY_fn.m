%%

function [timeDomainSig, frame, codedTrBlock, puschIndices, txInfo, txChInfo, rmcconfig] = ulPHY_fn(varargin)
    
    global ulChannel

    if(isempty(varargin))
        if(nargout >0 )
            % GUI call if no parameter is provided
            [timeDomainSig,frame,rmcconfig] = saLteRMCULToolGUI(nargout>0);
        else
            saLteRMCULToolGUI(nargout>0);
        end
    else
        if(isstruct(varargin{1}) && (nargin >= 2))
            rmc = varargin{1};
            transportStream = varargin{2};
            
            % Get CQI bits if provided
            if(nargin>2)
                cqiBits = varargin{3};
            else
                cqiBits = [];
            end
            
            % Get RI bits if provided
            if(nargin>3)
                riBits = varargin{4};
            else
                riBits = [];
            end
            
            % Get HARQ-ACK bits if provided
            if(nargin>4)
                ackBits = varargin{5};
            else
                ackBits = [];
            end
            
            if(~isfield(rmc,'RC'))
                error('lte:error','Input structure must contain the RC number');
            end
        elseif(ischar(varargin{1}) || isstring(varargin{1})) && (nargin >= 2)
            rmcNo = varargin{1};
            transportStream = varargin{2};

            % Get the number of subframes if supplied
            if(nargin == 4)
                totSubframes = varargin{4};
            else
                totSubframes = 10;
            end
            
            if(nargin > 2)
                duplexMode = varargin{3};
            else
                duplexMode = 'FDD';
            end
            cqiBits = [];
            riBits = [];
            ackBits = [];
            
            rmc.RC = rmcNo;
            rmc.TotSubframes = totSubframes;
            rmc.DuplexMode = duplexMode;
        else
            error('lte:error','Input argument(s) is not correct. Please see help.');
        end

        if(iscell(transportStream) || isempty(transportStream))
            [rows,cols] = size(transportStream);
            Ncodewords = max(rows,cols);            % Number of codewords
        elseif(isvector(transportStream)|| isempty(transportStream))
            Ncodewords = 1;
        end
        
        if(Ncodewords == 0 || isempty(transportStream))
            Ncodewords = 1;
        end
        
        % Get full RMC configuration for given inputs
        rmc = lteRMCUL(rmc);
 
        ttiPerBundle = 1;     % Set default value (No TTI bundling)
        if strcmpi(rmc.RC,'A11-1')
            ttiPerBundle = 4; % Enable TTI bundling for A11-1 RMC
        end
        
        PUSCH = rmc.PUSCH;
        
        % Input configuration structure is passed at output without any
        % parameter update so that the waveform can be regenerated with the
        % function output configuration
        rmcconfig = rmc;
        
        % Set up data generator source(s).
        if(iscell(transportStream))
            idx = cellfun(@isempty,transportStream);
            if(isempty(transportStream))
                dataSource1=saVectorDataSource([]);
                dataSource2=saVectorDataSource([]);
                PUSCH.TrBlkSizes(:) = 0;
                PUSCH.CodedTrBlkSizes(:) = 0;
            else
                dataSource1=saVectorDataSource(transportStream{1});
                if(numel(transportStream)>1)
                    dataSource2=saVectorDataSource(transportStream{2});
                elseif(idx)
                    PUSCH.TrBlkSizes(:) = 0;
                    PUSCH.CodedTrBlkSizes(:) = 0;
                end
            end
        else
            if(isempty(transportStream))
                PUSCH.TrBlkSizes(:) = 0;
                PUSCH.CodedTrBlkSizes(:) = 0;
            end
            dataSource1=saVectorDataSource(transportStream);
        end
        
        frame = [];
        reGrid = [];
        cBlkSizeIdx = 0;
        
        % Empty subframe construction for given configuration filled with all
        % zeros.
        subframe = lteULResourceGrid(rmc);
        
        totalSubframes = rmc.TotSubframes;
        if(~isnumeric(totalSubframes) || isempty(totalSubframes) || totalSubframes<0)
            error('lte:error','The function call resulted in an error: The total number of subframes must be a positive integer');
        end
        
        % Check if DL cyclic prefix is specified
        if(~isfield(rmc,'CyclicPrefix'))
            rmc.CyclicPrefix='Normal';
        end
        
        % Set default value to SerialCat if its not defined.
        if(~isfield(rmc,'SerialCat'))
            rmc.SerialCat = true;
        end
        serialCat = rmc.SerialCat;
        
        if(~isfield(PUSCH,'RVSeq') || (isfield(PUSCH,'RVSeq') && isempty(PUSCH.RVSeq)))
            if(Ncodewords == 2)
                PUSCH.RVSeq = [0;0];
            elseif(Ncodewords == 1)
                PUSCH.RVSeq = 0;
            end
        end
        
        % HARQ setup initialization
        noHarqProcesses = PUSCH.NHARQProcesses;
        harqTableIdx = 1;
        % If TTI bundling is enabled ('A11-1'), each HARQ process will
        % transmit 'ttiPerBundle' times with different RV values before
        % moving on to the next HARQ process. Define the HARQ table with
        % all HARQ processes repeating 'ttiPerBundle' times
        harqTable = reshape(repmat(1:noHarqProcesses,ttiPerBundle,1),1,noHarqProcesses*ttiPerBundle);
        
        harqprocess.data = struct('blk1',[],'blk2',[]);
        harqprocess.RVIdx = ones(Ncodewords,1);
        harqProcesses(1:max(harqTable)) = harqprocess;
        newData = [1 1];
        
        [nRows,nColms] = size(PUSCH.TrBlkSizes);
        tempBlkSize = PUSCH.TrBlkSizes;
        if(Ncodewords==2)
            if(nColms==2)
                tempBlkSize = tempBlkSize.';
            end
        else
            if(nRows>2)
                tempBlkSize = tempBlkSize.';
            end
        end
        
        % Get the absolute subframe number from NSubframe and NFrame
        NSubframe = rmc.NFrame*10+rmc.NSubframe;
        for subframeIdx=NSubframe:(NSubframe+totalSubframes)-1

            % Update subframe number and clear subframe
            subframe(:)=0;
            rmc.NSubframe = mod(subframeIdx,10);
            rmc.NFrame = floor(subframeIdx/10);

            % If this subframe is a cell-specific SRS subframe, configure
            % the PUCCH for shortened transmission.
            if(isfield(rmc,'SRS') && isstruct(rmc.SRS))
                srsInfo = lteSRSInfo(rmc,rmc.SRS);
                rmc.Shortened = srsInfo.IsSRSSubframe;
            end
            
            info=lteDuplexingInfo(rmc);
            transportBlock = {[]}; % Initialize cw1
            if (info.NSymbolsUL)
                % Get the transport block size for the current subframe
                if(Ncodewords==2)
                    transportBlkSize = tempBlkSize(:,mod(rmc.NSubframe,size(tempBlkSize,2)) + 1);
                    if(size(transportBlkSize,1)==2)
                        transportBlkSize = transportBlkSize.';
                    elseif(size(transportBlkSize,1)==1 && size(transportBlkSize,2)==1)
                        transportBlkSize = [transportBlkSize transportBlkSize]; %#ok<AGROW>
                    end
                    
                    if(isempty(transportBlkSize(1)) || isempty(transportBlkSize(2)))
                        effectiveCodewords = 1;
                        transportBlkSize(1) = transportBlkSize(1)*layersPerCW(PUSCH.NLayers,effectiveCodewords);
                        transportBlkSize(2) = transportBlkSize(2)*layersPerCW2(PUSCH.NLayers,effectiveCodewords);
                    end
                    transportBlock{2} = [];% Initialize cw2
                else
                    transportBlkSize = tempBlkSize(1,mod(rmc.NSubframe,size(tempBlkSize,2)) + 1);
                end
                
                % Generate PUSCH, PUSCH can only be transmitted in UL subframes
                if(strcmpi(info.SubframeType,'Uplink')) && any(transportBlkSize~=0)
                    harqIdx =  harqTable(harqTableIdx);
                    tempRv(1,1) = PUSCH.RVSeq(1,harqProcesses(harqIdx).RVIdx(1));
                    if(Ncodewords==2)
                        tempRv(2,1) = PUSCH.RVSeq(end,harqProcesses(harqIdx).RVIdx(2));
                    end
                    PUSCH.RV = tempRv(:,1);
                    newData(:) = (harqProcesses(harqIdx).RVIdx == 1);

                    % UL-SCH transport block size of configured RMC as in TS 36.104.
                    if(Ncodewords==1 && (newData(1) || isempty(harqProcesses(harqIdx).data.blk1)))
                        harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize);
                    elseif(Ncodewords == 2)
                        transportBlkSize(idx) = 0;
                        if(isempty(PUSCH.RV))
                            PUSCH.RV = 0;
                        end
                        if(newData(1) || isempty(harqProcesses(harqIdx).data.blk1))
                            harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize(1));
                        end
                        if(newData(2) || isempty(harqProcesses(harqIdx).data.blk2))
                            harqProcesses(harqIdx).data.blk2 = dataSource2.getData(transportBlkSize(2));
                        end
                    end
                    if(transportBlkSize(1)~=0)
                        transportBlock = {harqProcesses(harqIdx).data.blk1};
                    end
                    if(Ncodewords == 2)
                        if(transportBlkSize(2)~=0)
                            transportBlock{2} = harqProcesses(harqIdx).data.blk2;
                        end
                    end
                    % If transport block size wasn't 0, update rv and
                    % harq table index
                    harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) = mod(harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) ,size(PUSCH.RVSeq,2))+1;
                    harqTableIdx = mod(harqTableIdx,length(harqTable)) + 1;
                    
                    % Generating UL-SCH coded bits by performing complete channel coding including
                    % CRC calculation, code block segmentation and CRC attachment, turbo
                    % coding, rate matching and code block concatenation.
                    codedTrBlock = lteULSCH(rmc,PUSCH,transportBlock,cqiBits,riBits,ackBits);

                    cBlkSizeIdx = cBlkSizeIdx+1;
                    if(Ncodewords==2)
                        PUSCH.TrBlkSizes(1,cBlkSizeIdx) = transportBlkSize(1);
                        PUSCH.TrBlkSizes(2,cBlkSizeIdx) = transportBlkSize(2);
                        PUSCH.CodedTrBlkSizes(1,cBlkSizeIdx) = length(codedTrBlock{1});
                        PUSCH.CodedTrBlkSizes(2,cBlkSizeIdx) = length(codedTrBlock{2});
                    else
                        PUSCH.TrBlkSizes(1,cBlkSizeIdx) = transportBlkSize;
                        PUSCH.CodedTrBlkSizes(1,cBlkSizeIdx) = length(codedTrBlock{1});
                    end

                    % Complex-valued modulated symbol generation for PUSCH. This involves
                    % scrambling, modulation and precoding processes.
                    if(Ncodewords==2 && isempty(codedTrBlock{1}) && ~isempty(codedTrBlock{2}))
                        % Swapping codeword and corresponding modulation scheme
                        puschSymbols = ltePUSCH(rmc,setfield(PUSCH,'Modulation',{PUSCH.Modulation{2} PUSCH.Modulation{1}}),{codedTrBlock{2} []}); %#ok<SFLD>
                    else
                        puschSymbols = ltePUSCH(rmc,PUSCH,codedTrBlock);
                    end
                    
                    % PUSCH symbols mapping on to the resource grid 
                    puschIndices = ltePUSCHIndices(rmc,PUSCH);
                    subframe(puschIndices) = puschSymbols;  
                    
                    % PUSCH DRS symbol creation and mapping on to resource grid
                    puschDrsSeq = ltePUSCHDRS(rmc,PUSCH);
                    puschDrsSeqIndices = ltePUSCHDRSIndices(rmc,PUSCH);                
                    subframe(puschDrsSeqIndices) = puschDrsSeq;        
                end
                
                if(isfield(rmc,'SRS') && isstruct(rmc.SRS) && srsInfo.IsSRSSubframe)
                    % Transmit SRS (if active under UE-specific SRS and cell-specific configuration)
                    [srsIndices,srsIndicesInfo] = lteSRSIndices(rmc,rmc.SRS);
                    srsSymbols = lteSRS(rmc,rmc.SRS);
                    if (rmc.SRS.NTxAnts==1 && rmc.NTxAnts>1)
                        subframe(offsetIndices(rmc,srsIndices,srsIndicesInfo.Port)) = srsSymbols;
                    else
                        subframe(srsIndices) = srsSymbols;
                    end
                end
            end
            % concatenate subframes to form a complete frame
            frame = cat(2,frame,subframe);
            if(~serialCat)
                reGrid(:,:,:,mod(subframeIdx-NSubframe,NSubframe+totalSubframes)+1) = subframe; %#ok<AGROW>
            end
        end
        
        if(~isfield(rmc,'Windowing'))
            rmc.Windowing = 0;
        end
        
        % Time Domain mapping by performing SC-FDMA modulation for uplink
        % symbols.
        [timeDomainSig,txInfo] = lteSCFDMAModulate(rmc,frame);
        if(~serialCat)
            frame = reGrid;
        end
        rmcconfig.SamplingRate = txInfo.SamplingRate;
        rmcconfig.Nfft = txInfo.Nfft;
        rmcconfig.Windowing = txInfo.Windowing;
        rmcconfig.PUSCH.HARQProcessSequence = getHARQTable(rmcconfig,ttiPerBundle);
        
        % Apply fading channel
        % ------------------------- %
        subframeNo = NSubframe;
        ulChannel.ulchcfg.InitTime = subframeNo/1000;
        ulChannel.ulchcfg.SamplingRate  = txInfo.SamplingRate;

        [timeDomainSig, txChInfo] = ulChannel.lteCh_fn(timeDomainSig, 'Uplink');
%         txChInfo = [];
        % =============================================================== %
    end
end

% This function calculates the number of layers per codeword.
function nLyrs = layersPerCW(layers,nCWs)
    nLyrs = floor(layers/nCWs);
end

function nLyrs = layersPerCW2(layers,nCWs)
    nLyrs = ceil(layers/nCWs);
end

% obtains indices 'out' corresponding to the indices 'in' offset to address
% antenna plane 'p'. 
function out = offsetIndices(ue,in,p)
    griddims=lteULResourceGridSize(ue);
    K=griddims(1);          % number of subcarriers
    L=griddims(2);          % number of OFDM symbols in a subframe    
    out=in+(K*L*p);
end

% Function to calculate the HARQ table assuming all uplink subframes are
% carrying data and have the same transport block size
function harqTable = getHARQTable(rmc,ttiPerBundle)
    noHarqProcesses = rmc.PUSCH.NHARQProcesses;
    % Define the default values for cyclic prefix
    if ~isfield(rmc,'CyclicPrefixUL')
        rmc.CyclicPrefixUL = 'Normal';
    end
    if ~isfield(rmc,'CyclicPrefix')
        rmc.CyclicPrefix = 'Normal';
    end
    info = arrayfun(@(x)lteDuplexingInfo(setfield(rmc,'NSubframe',x)),0:9); 
    activesfs = arrayfun(@(x)strcmpi(x.SubframeType,'Uplink'),info);
    harqTable = ones(10,lcm(noHarqProcesses*ttiPerBundle,sum(activesfs))/sum(activesfs))*-1;
    harqTable(activesfs==0,:) = 0; % Non-transmititng subframes
    harqTable(harqTable==-1) = repmat(repmat(1:noHarqProcesses,ttiPerBundle,1),1,lcm(noHarqProcesses*ttiPerBundle,sum(activesfs))/(noHarqProcesses*ttiPerBundle));
	harqTable = harqTable(:).';
end