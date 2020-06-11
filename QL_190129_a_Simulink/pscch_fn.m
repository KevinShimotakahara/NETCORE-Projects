
%%

function [waveform] = pscch_fn(nodeID, pscchPeriod, dci, sfNumber)

% update subframe number
sfNumber = mod(sfNumber, pscchPeriod.Config.sc_Period_r12)+1;

% Preallocate output for PSCCH period resource grid and baseband waveform
ue.NSLRB = pscchPeriod.Config.NSLRB;
ue.CyclicPrefixSL = pscchPeriod.Config.syncConfig.syncCP_Len_r12;
sfdims = lteSLResourceGridSize(ue);
grid = zeros(sfdims(1), sfdims(2));
waveinfo = lteSLSCFDMAInfo(ue);
nsamplespersf = waveinfo.SamplingRate/1000;
waveform = zeros(nsamplespersf, 1);

% The SC-FDMA modulation will not use windowing due to the variable CP
% and the piece-wise nature of the waveform construction
windowing = 0;

% Create basic set of common parameters used by the low-level PHY functions
ue = struct('NSLRB', pscchPeriod.Config.NSLRB, 'CyclicPrefixSL', pscchPeriod.Config.sc_CP_Len_r12,...
    'DuplexMode', pscchPeriod.Config.DuplexMode, 'TDDConfig', pscchPeriod.Config.TDDConfig);

ue.SidelinkMode = 'D2D';

% Get the resources (PRB/subframes) that will be used to carry
% the pair of PSCCH which in turn carry the coded SCI
[pscchsubframes, pscchprb] = pscchPeriod.getPSCCHResources(dci);
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
if sfNumber ~= (pscchPeriod.SyncSubframes+1)
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

end

