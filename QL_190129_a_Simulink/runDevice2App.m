function runDevice2App(event,varargin)
    global UDs
    global baseComfortD2
    global PriceReadingD2
    global demandReductionD2
    global currentDemandD2
    global currentDemandFunk
    global gClk
    
    switch event
        case 1
            %Send demand reduction availability update
            %instead of random, do power reading
            lengthPfield = 14;
            currentDemandD2 = varargin{1};
            reductionAvailability = currentDemandD2 - demandReductionD2 - baseComfortD2;
            if reductionAvailability < 0
                reductionAvailability = 0;
            end
            ra = num2binArray(round(reductionAvailability),lengthPfield);
            currentDemand = num2binArray(round(currentDemandD2),lengthPfield);
            if isempty(UDs(2).pktNumber)
                UDs(2).pktNumber = 1;
            else
                UDs(2).pktNumber = UDs(2).pktNumber + 1;
            end
            rlcSDU = [zeros(1,77), ra, currentDemand, randi([0,1],1,2000-105)]; %
            rlcSDU = configMetaDataHeader(UDs(2),rlcSDU);

            UDs(2).RLCtxEntity.rxRLCSDU(rlcSDU);

        case 2
            %Recieved message(s)
            for k = 1:length(UDs(2).rxBuffer)
                sdu = UDs(2).rxBuffer{k};
                PriceReadingD2 = binArray2num(sdu(78:91));
                %check if price is high enough to self-regulate demand
                if PriceReadingD2 > 350
                   demandReductionD2 = currentDemandFunk(gClk/1000) - baseComfortD2;
                else
                   demandReductionD2 = 0;
                end
                rxedDemandReduction = binArray2num(sdu(92:105));
                
                if rxedDemandReduction > demandReductionD2
                    if rxedDemandReduction > currentDemandFunk(gClk/1000) - baseComfortD2
                        demandReductionD2 = currentDemandFunk(gClk/1000) - baseComfortD2;
                    else
                        demandReductionD2 = rxedDemandReduction;
                    end
                end
            end
            
        otherwise
    end
end