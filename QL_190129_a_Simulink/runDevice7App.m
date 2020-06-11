function runDevice7App(eventID)
    global UDs
%     try
        switch eventID

            case 1
                if isempty(UDs(7).pktNumber)
                    UDs(7).pktNumber = 1;
                else
                    UDs(7).pktNumber = UDs(7).pktNumber + 1;
                end
                rlcSDU = [zeros(1,77),randi([0,1],1, 1120-77)];
                rlcSDU = configMetaDataHeader(UDs(7),rlcSDU);

                if isempty(UDs(7).RLCtxEntity.TXbuff)
                    UDs(7).RLCtxEntity.rxRLCSDU(rlcSDU);
                else
                    UDs(7).RLCtxEntity.rxRLCSDU(rlcSDU);
                end

            otherwise
        end
%     catch
%         diversion = 1;
%     end
end