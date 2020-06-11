function runDevice5App(eventID)
    global UDs
    switch eventID
        case 1
                rlcSDU = [zeros(1,77), randi([0,1],1, 592-77)]; %Set to IEEE C37.188.2-2011 Standard Payload for Synchrophasor Data Transfer
                rlcSDU = configMetaDataHeader(UDs(5),rlcSDU);

                if isempty(UDs(5).pktNumber)
                    UDs(5).pktNumber = 1;
                else
                    UDs(5).pktNumber = UDs(5).pktNumber + 1;
                end

                if isempty(UDs(5).RLCtxEntity.TXbuff)
                    UDs(5).RLCtxEntity.rxRLCSDU(rlcSDU);
                else
                    UDs(5).RLCtxEntity.rxRLCSDU(rlcSDU);
                end
        otherwise
    end
end