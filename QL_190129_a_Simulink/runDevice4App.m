function runDevice4App(event,varargin)
    global UDs
    global baseComfortD4
    global PriceReadingD4
    global demandReductionD4
    global currentDemandD4
    global currentDemandFunk
    global gClk
    
    switch event
        case 1
            %Send demand reduction availability update
            %instead of random, do power reading
            lengthPfield = 14;
            currentDemandD4 = varargin{1};
            reductionAvailability = currentDemandD4 - demandReductionD4 - baseComfortD4;
            if reductionAvailability < 0
                reductionAvailability = 0;
            end
            ra = num2binArray(round(reductionAvailability),lengthPfield);
            currentDemand = num2binArray(round(currentDemandD4),lengthPfield);
            if isempty(UDs(4).pktNumber)
                UDs(4).pktNumber = 1;
            else
                UDs(4).pktNumber = UDs(4).pktNumber + 1;
            end
            rlcSDU = [zeros(1,77), ra, currentDemand, randi([0,1],1,2000-105)]; %
            rlcSDU = configMetaDataHeader(UDs(4),rlcSDU);

            UDs(4).RLCtxEntity.rxRLCSDU(rlcSDU);

        case 2
            %Recieved message(s)
            for k = 1:length(UDs(4).rxBuffer)
                sdu = UDs(4).rxBuffer{k};
                PriceReadingD4 = binArray2num(sdu(78:91));
                %check if price is high enough to self-regulate demand
                if PriceReadingD4 >= 300
                   demandReductionD4 = currentDemandFunk(gClk/1000) - baseComfortD4;
                else
                   demandReductionD4 = 0;
                end
                rxedDemandReduction = binArray2num(sdu(120:133));
                
                if rxedDemandReduction > demandReductionD4
                    if rxedDemandReduction > currentDemandFunk(gClk/1000) - baseComfortD4
                        demandReductionD4 = currentDemandFunk(gClk/1000) - baseComfortD4;
                    else
                        demandReductionD4 = rxedDemandReduction;
                    end
                end
            end
            
        otherwise
    end
end