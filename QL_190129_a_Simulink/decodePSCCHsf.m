%%

function [sciMessageRx, sciDecoded] = decodePSCCHsf(rxWaveform, udTemp, pscch, sfNumber)

sciMessageRx = [];
    
% If SCI BLER measurement is configured:
sciDecoded = false;

% If this subframe is a PSCCH, synchronize, extract and decode
% the subframe to get the SCI message.
if(~isempty(sfNumber))
    subframeWaveform = rxWaveform;
    bestCorr = 0.0;
    for prb = udTemp.pscchPeriod.PSCCHResourceBlockPool.'
        pscch.PRBSet = prb;
        [thisoffset, thiscorr] = lteSLFrameOffsetPSCCH(pscch, subframeWaveform);
        % Find the best correlation across all antennas and
        % store the corresponding offset
        maxCorr = max(thiscorr(:));
        if maxCorr > bestCorr
            bestCorr = maxCorr;
            offset = thisoffset;
        end
    end
    % assume perfect synchronization
    offset = 0;     % Fixme
    subframe = lteSLSCFDMADemodulate(pscch, rxWaveform(offset+1:end, :));
    
    % Repeat for each PSCCH resource until the SCI is decoded
    pscchResource = 0;
    while (pscchResource < udTemp.pscchPeriod.NumPSCCHResource && ~sciDecoded)
        
        % Get the PSCCH subframes 'sf' and PRB allocations 'prb'
        % for the current PSCCH resource
        [sf, prb] = udTemp.pscchPeriod.getPSCCHResources(pscchResource);
        nsfpscchPos = find(sfNumber == sf);
        
        % Remove any resources that overlap with synchronization
        % subframes. PSCCH will not be transmitted in these
        % subframes
        [~,syncIndex] = intersect(sf,udTemp.pscchPeriod.SyncSubframes);
        prb(syncIndex) = [];
        
        % Configure the PSCCH receiver for the PRB allocation
        pscch.PRBSet = prb(:, nsfpscchPos);

        % Perform channel estimation, extract the received
        % PSCCH symbols and the corresponding channel estimate,
        % and perform equalization
        cec.PilotAverage = 'UserDefined';
        cec.TimeWindow = 15;
        cec.FreqWindow = 23;
        cec.InterpType = 'linear';

        [hest,nest] = lteSLChannelEstimatePSCCH(pscch, cec, subframe);
        [pscchIndices,pscchIndicesInfo] = ltePSCCHIndices(pscch);
        [pscchRx,pscchHest] = lteExtractResources(pscchIndices,subframe,hest);
        pscchSymbols = lteEqualizeMMSE(pscchRx,pscchHest,nest);

        % If we are receiving the first PSCCH transmission
        % instance, reset the receiver buffer
        codedSciBits = zeros(pscchIndicesInfo.G,1);

        % Demodulate the PSCCH and add the result into the
        % receiver buffer
        codedSciBits = codedSciBits + ltePSCCHDecode(pscchSymbols);

        % Decode the SCI message. If successful (CRC=0), check
        % the NSAID field of the decoded message against the
        % transmitted message
        sciInfo = lteSCIInfo(pscch);
        [sciBits,sciCRC] = lteSCIDecode(sciInfo.Format0,codedSciBits);
        if (sciCRC==0)
            sciMessageRx = lteSCI(pscch,sciBits);
            if (sciMessageRx.NSAID == udTemp.expectedNSAID)
                sciDecoded = true;
            else
                sciMessageRx = [];
            end
        end            
        
        % Increment the PSCCH resource number
        pscchResource = pscchResource + 1;
        
    end
end



end